using Test
using ForwardDiff
using SparseArrays
using Pyrolysis
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.Kernels
using Pyrolysis.Kernels: Contributions, contribute!, find_kernel
using Pyrolysis.Assembly
using Pyrolysis.Assembly: residual!
using Pyrolysis.ColoredJacobian
using Pyrolysis.Integration
using Pyrolysis.Integration: setup, integrate, consistent_state!, CoreProblem
using Pyrolysis.Closures: EquidistributionCorrection

const _L = 0.01

# Two solids of equal density and different heat capacity (`c_b` the second
# one's), a smooth front between them, adiabatic and impermeable faces; with
# `law` the material carries a lateral shrinkage law (and so the history pair).
function _composition(n; T0 = 800.0, law = nothing, c_b = 3000.0)
    a = SolidComponent(:a, ρ = 1000.0, c = 1000.0, k = 0.5)
    b = SolidComponent(:b, ρ = 1000.0, c = c_b, k = 0.5)
    mat = Material(
        name = :composition,
        components = (a, b),
        reactions = (),
        lateral_shrinkage_law = law,
    )
    front(z) = (1 + tanh((z - _L / 2) / (_L / 10))) / 2
    return PyrolysisProblem(
        mesh = create_uniform_mesh(_L, n),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
        ),
        T_initial = z -> T0,
        ξ_initial = [z -> 1000.0 * (1 - front(z)), z -> 1000.0 * front(z)],
        tspan = (0.0, 10.0),
    )
end

# The mesh slides over the material: the end nodes stay, the middle moves a
# quarter of the thickness up over the solve.
const _slide = (x, t) -> _L / 40 * sinpi(x / _L)

@testset "van Albada slope" begin
    @test van_albada(2.0, 2.0, 1e-12) ≈ 2.0
    @test van_albada(0.0, 0.0, 1e-12) == 0.0
    @test abs(van_albada(1.0, -1.0, 1e-12)) < 1e-12
    @test 0 < van_albada(1.0, 3.0, 1e-12) < 3.0
    @test all(
        isfinite,
        ForwardDiff.gradient(x -> van_albada(x[1], x[2], 1e-12), [0.0, 0.0]),
    )
end

# The relative-flux kernel without its enthalpy flux: the negative control of
# isothermal advection. Runs after the kernel and removes what it added.
struct _WithoutEnergy{K} <: Pyrolysis.Kernels.Kernel
    kernel::K
end
function Pyrolysis.Kernels.contribute!(
    c::Contributions{NC,T},
    k::_WithoutEnergy,
    ds,
    topo,
    material,
    t,
) where {NC,T}
    scratch = Contributions{T}(NC, 0, n_cells(topo), n_faces(topo), n_nodes(topo))
    contribute!(scratch, k.kernel, ds, topo, material, t)
    c.F_energy .-= scratch.F_energy
    return c
end

@testset "a NaN trial temperature gives a NaN residual through a tabulated enthalpy: $(energy) form" for energy in
                                                                                                         (
    :temperature,
    :enthalpy,
)
    # The second solid's heat capacity has no closed-form integral, so its
    # enthalpy, the relative fluxes' face enthalpies among them, is tabulated.
    c_b = CallableProperty(T -> 3000.0 + 0.2T * sin(T / 300))
    cp = setup(
        _composition(8; c_b);
        energy_state = energy,
        use_ale = true,
        w_override = _slide,
    )
    @test find_kernel(RelativeFluxKernel, cp.def.kernels) !== nothing
    @test cp.def.setup.laws.h[2] isa TabulatedEnthalpy
    u = copy(cp.u0)
    u[energy_index(cp.def.layout, 4)] = NaN
    du = zeros(length(u))
    residual!(du, u, cp.def, 0.0, cp.ws)
    @test any(isnan, du)
end

@testset "isothermal advection: $(energy) form" for energy in (:temperature, :enthalpy)
    cp = setup(
        _composition(20);
        energy_state = energy,
        use_ale = true,
        w_override = _slide,
        abstol = 1e-10,
        reltol = 1e-8,
    )
    sol = integrate(cp; saveat = Float64[])
    @test sol.retcode === :Success
    @test length(sol.t) > 10
    @test maximum(abs, sol.T .- 800.0) <= 1e-9
    # The front moved through the mesh: a cell near it changed composition.
    change = maximum(abs, sol.ξ[1, :, end] .- sol.ξ[1, :, 1])
    @test change > 300
    # Mass is conserved; the end nodes did not move.
    @test sum(sol.m[:, :, end]) ≈ sum(sol.m[:, :, 1]) rtol = 1e-12
    @test sol.z[1, end] == 0.0
    @test sol.z[end, end] ≈ _L rtol = 1e-12

    # Without the enthalpy the mass fluxes carry, the composition change moves
    # the temperature.
    def = cp.def
    relative = find_kernel(RelativeFluxKernel, def.kernels)
    control = CoreDefinition(
        def.layout,
        def.setup,
        (def.kernels..., _WithoutEnergy(relative)),
        def.closure,
        def.sources,
    )
    ws = CoreWorkspace(control, cp.T0)
    cj = ColoredJacobianCache(control, ws, cp.u0, 0.0)
    cpc = CoreProblem(
        cp.problem,
        control,
        ws,
        cp.u0,
        cp.reltol,
        cp.abstol,
        cj,
        cp.T0,
        cp.m0,
        cp.deactivation,
    )
    solc = integrate(cpc; saveat = Float64[])
    @test solc.retcode === :Success
    @test maximum(abs, solc.T .- 800.0) > 100
end

# At 8 cells the top cell starts with 0.158 kg/m³ of the first solid beside
# the front's steep side, and the mesh crosses it inward: its face state
# toward the second cell has no cell beyond it.
@testset "a boundary cell beside a steep front is not drained: $(energy) form" for energy in
                                                                                   (
    :temperature,
    :enthalpy,
)
    cp = setup(
        _composition(8);
        energy_state = energy,
        use_ale = true,
        w_override = _slide,
        abstol = 1e-10,
        reltol = 1e-8,
    )
    sol = integrate(cp; saveat = Float64[])
    @test sol.retcode === :Success
    @test sol.t[end] == 10.0
    @test minimum(sol.ξ) > 0
    @test sum(sol.m[:, :, end]) ≈ sum(sol.m[:, :, 1]) rtol = 1e-12
end

@testset "an inventory's face state lies in [0, 2p] and is zero for an empty cell" begin
    # Cells 1, 2, 3 at unit spacing (X, U, D), the face half a cell from U.
    reconstruction(q) = Kernels._face_value(i -> q[i], 2, 3, 1, 0.5, 1.0, 1.0, 1e-12)
    state(q) = Kernels._inventory_face_value(i -> q[i], 2, 3, 1, 0.5, 1.0, 1.0, 1e-12, 1e-6)
    # Between full cells, the steeper gradient behind the empty cell gives
    # the reconstruction a positive slope; the face sends nothing.
    @test reconstruction((2.0, 0.0, 1.0)) > 0
    @test state((2.0, 0.0, 1.0)) == 0
    @test state((1.0, -1e-9, 3.0)) == 0
    for q in ((2.0, 1e-3, 1.0), (0.0, 0.5, 10.0), (10.0, 0.5, 0.0), (1.0, 2.0, 3.0))
        @test 0 <= state(q) <= 2q[2]
    end
    # A small increment against the inventory keeps the reconstruction.
    @test state((999.0, 1000.0, 1001.0)) ≈ reconstruction((999.0, 1000.0, 1001.0)) rtol =
        1e-9
end

# Four inert cells of two solids, the first one's concentrations 2, 0, 1 and
# 1 kg/m³: the second cell holds none of it between cells that do.
function _interior_gap()
    a = SolidComponent(:a; ρ = 1000.0, c = 1000.0, k = 0.5)
    b = SolidComponent(:b; ρ = 1000.0, c = 1000.0, k = 0.5)
    q(z) = z < 0.25 ? 2.0 : z < 0.5 ? 0.0 : 1.0
    return PyrolysisProblem(
        mesh = create_uniform_mesh(1.0, 4),
        material = Material(name = :interior_gap, components = (a, b), reactions = ()),
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
        ),
        T_initial = z -> 600.0,
        ξ_initial = [q, z -> 1000.0 - q(z)],
        tspan = (0.0, 0.02),
        dt_initial = 1e-4,
        output_times = [0.0, 0.02],
    )
end

# The middle node moves at −0.01 m/s and the others stay, so the mesh crosses
# the empty cell's face toward the cell beyond it.
const _gap_motion = (x, t) -> -0.64 * x * (x - 0.25) * (x - 0.75) * (x - 1.0)

@testset "an empty interior cell is not drained: $(energy) form" for energy in (
    :temperature,
    :enthalpy,
)
    cp = setup(
        _interior_gap();
        energy_state = energy,
        use_ale = true,
        w_override = _gap_motion,
    )
    du = similar(cp.u0)
    residual!(du, cp.u0, cp.def, 0.0, cp.ws)
    k = mass_index(cp.def.layout, 2, 1)
    @test cp.u0[k] == 0
    @test du[k] == 0
    sol = integrate(cp)
    @test sol.retcode === :Success
    @test sol.t[end] == 0.02
    @test minimum(sol.m) >= 0
    @test sum(sol.m[:, :, end]) ≈ sum(sol.m[:, :, 1]) rtol = 1e-12
end

@testset "an empty cell is not drained under the equidistribution correction" begin
    s = SolidComponent(:solid; ρ = 1000.0, c = 1000.0, k = 0.5)
    g = GaseousComponent(:gas; M = 0.032, c = 1000.0, k = 0.03, λ = 0.0)
    q(z) = z < 0.05 ? 2.0 : z < 0.1 ? 0.0 : 1.0
    problem = PyrolysisProblem(
        mesh = create_uniform_mesh(0.2, 4),
        material = Material(name = :gap_correction, components = (s, g), reactions = ()),
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (top = ImpermeableBC(), bottom = ImpermeableBC()),
        ),
        T_initial = z -> 600.0,
        ξ_initial = [z -> 1000.0, q],
        tspan = (0.0, 0.2),
        dt_initial = 1e-4,
        output_times = [0.0, 0.2],
    )
    cp = setup(
        problem;
        use_ale = true,
        correction = EquidistributionCorrection(layers = :all),
    )
    # The third cell compressed to a fifth of its thickness, which activates
    # the correction at its nodes: the mesh crosses the empty second cell's
    # outer face into it.
    layout = cp.def.layout
    u = copy(cp.u0)
    x = [0.0, 0.05, 0.1, 0.11, 0.2]
    for (nd, xn) in enumerate(x)
        u[position_index(layout, nd)] = xn
    end
    for i = 1:4
        V = x[i+1] - x[i]
        u[mass_index(layout, i, 1)] = 1000.0 * V
        u[mass_index(layout, i, 2)] = (2.0, 0.0, 1.0, 1.0)[i] * V
    end
    consistent_state!(u, cp.def, cp.ws, 0.0)
    du = similar(u)
    residual!(du, u, cp.def, 0.0, cp.ws)
    @test cp.ws.derived.c[3] < 0
    k = mass_index(layout, 2, 2)
    @test u[k] == 0
    @test du[k] == 0
    sol = integrate(cp; u0 = u)
    @test sol.retcode === :Success
    @test sol.t[end] == 0.2
    @test minimum(sol.m) >= 0
end

# The front at rest under the slide: each cell holds the front's integral over
# its own nodes at every time, so the cell masses have an exact reference.
@testset "the relative fluxes converge at second order on a front at rest" begin
    w = _L / 10
    I(z) = 1000.0 * (z / 2 + (w / 2) * log(cosh((z - _L / 2) / w)))
    errors = map((32, 64, 128)) do n
        cp = setup(
            _composition(n);
            use_ale = true,
            w_override = _slide,
            abstol = 1e-12,
            reltol = 1e-10,
        )
        sol = integrate(cp; saveat = Float64[])
        @test sol.retcode === :Success
        z = sol.z[:, end]
        b = [I(z[i+1]) - I(z[i]) for i = 1:n]
        a = 1000.0 .* diff(z) .- b
        return (sum(abs, sol.m[1, :, end] .- a) + sum(abs, sol.m[2, :, end] .- b)) /
               (2000.0 * _L)
    end
    @test all(>(1.6), log2.(errors[1:(end-1)] ./ errors[2:end]))
end

# One solid at rest at a uniform concentration, compressed by a prescribed
# uniform stretch of the mesh.
function _uniform(n; coordinates = Planar())
    s = SolidComponent(:s, ρ = 800.0, c = 1000.0, k = 0.3)
    return PyrolysisProblem(
        mesh = create_uniform_mesh(_L, n; coordinates),
        material = Material(name = :uniform, components = (s,), reactions = ()),
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
        ),
        T_initial = z -> 500.0,
        ξ_initial = [z -> 600.0],
        tspan = (0.0, 10.0),
        output_times = [0.0, 10.0],
    )
end
_stretch(t) = 1 - 0.15 * (1 - cos(2π * t / 20))
const _compress = (x, t) -> (-0.15 * (2π / 20) * sin(2π * t / 20) / _stretch(t)) * x

function _free_stream_error(coords, dt; energy = :temperature)
    sol = Integration.solve(
        _uniform(12; coordinates = coords);
        energy_state = energy,
        use_ale = true,
        w_override = _compress,
        dt,
        adaptive = false,
        abstol = 1e-12,
        reltol = 1e-12,
    )
    @test sol.retcode === :Success
    @test sol.z[end, end] ≈ _L * _stretch(10.0) rtol = 1e-6
    return maximum(abs, sol.ξ[1, :, :] ./ 600 .- 1), maximum(abs, sol.T .- 500)
end

@testset "free stream: exact on the planar mesh, at the integrator's order on radial ones" begin
    for energy in (:temperature, :enthalpy), dt in (0.2, 0.1, 0.05)
        ξ_err, T_err = _free_stream_error(Planar(), dt; energy)
        @test ξ_err <= 1e-13
        @test T_err <= 1e-10
    end
    for coords in (Cylindrical(), Spherical())
        errors = [_free_stream_error(coords, dt)[1] for dt in (0.8, 0.4, 0.2)]
        orders = log2.(errors[1:2] ./ errors[2:3])
        @info "free-stream error, $(nameof(typeof(coords)))" errors orders
        @test errors[3] > 1e-11
        @test all(p -> p > 3.5, orders)
    end
end

@testset "history transport under a prescribed mesh velocity" begin
    cp = setup(
        _composition(16; law = χ -> 1 - 0.3χ);
        use_ale = true,
        w_override = _slide,
        abstol = 1e-10,
        reltol = 1e-8,
    )
    layout = cp.def.layout
    u0 = copy(cp.u0)
    for i = 1:n_cells(layout)
        u0[history_index(layout, i, 1)] =
            0.4 * i / n_cells(layout) * u0[history_index(layout, i, 2)]
    end
    consistent_state!(u0, cp.def, cp.ws, 0.0)
    sums(u) = (
        sum(u[history_index(layout, i, 1)] for i = 1:n_cells(layout)),
        sum(u[history_index(layout, i, 2)] for i = 1:n_cells(layout)),
    )
    G0, M0 = sums(u0)
    sol = integrate(cp; u0, saveat = Float64[])
    @test sol.retcode === :Success
    for u in sol.raw_solution.u
        G, M = sums(u)
        @test G ≈ G0 rtol = 1e-12
        @test M ≈ M0 rtol = 1e-12
        @test u[area_index(layout)] ≈ u0[area_index(layout)] rtol = 1e-12
    end
    # The history moved with the solid relative to the mesh.
    first_G = [u0[history_index(layout, i, 1)] for i = 1:n_cells(layout)]
    last_G =
        [sol.raw_solution.u[end][history_index(layout, i, 1)] for i = 1:n_cells(layout)]
    @test maximum(abs, last_G .- first_G) > 1e-3 * maximum(first_G)
end

@testset "relative fluxes: containment, $(energy) form, $(n) cells, $(name)" for energy in
                                                                                 (
        :temperature,
        :enthalpy,
    ),
    n in (6, 8),
    name in (:plain, :history)

    law = name === :history ? (χ -> 1 - 0.3χ) : nothing
    cp = setup(
        _composition(n; law);
        energy_state = energy,
        use_ale = true,
        w_override = _slide,
    )
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    layout = def.layout
    u0 = copy(cp.u0)
    if has_chi(layout)
        for i = 1:n
            u0[history_index(layout, i, 1)] = 0.1i * u0[history_index(layout, i, 2)]
        end
    end
    du = zeros(length(u0))
    residual!(du, u0, def, 0.0, ws)
    @test (@allocated residual!(du, u0, def, 0.0, ws)) == 0
    hot = copy(u0)
    for i = 1:n
        hot[energy_index(layout, i)] =
            energy_state(layout) isa TemperatureState ? 700.0 + 20i :
            u0[energy_index(layout, i)] * 1.5
    end
    # A trial mesh velocity that changes sign along the column, so both
    # upwind sides of the face states are exercised.
    for i = 2:n
        hot[velocity_index(layout, i)] = 1e-4 * (-1)^i
    end
    perturbed = u0 .* (1 .+ 0.2 .* sin.(1:length(u0)))
    J = jac_prototype(cj)
    for (state, u) in
        pairs((; u0, hot, perturbed, zero_state = zeros(length(u0)), negative = -abs.(u0)))
        @testset "$state" begin
            residual!(du, u, def, 1.0, ws)
            @test all(isfinite, du)
            Jd = dense_jacobian(def, ws, u, 1.0)
            @test all(isfinite, Jd)
            @test isempty(missing_entries(Jd, cj.pattern))
            jacobian!(J, cj, u, 1.0)
            @test compare_entries(Matrix(J), Jd; rtol = 1e-12).passed
        end
    end
end
