using Test
using Pyrolysis
using Pyrolysis.Materials: DepletionLimiter, RADIATIVE_TRANSFER, SURFACE_ABSORPTION
using Pyrolysis.Materials: STEFAN_BOLTZMANN
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.Kernels:
    RadiativeTransferKernel,
    MirrorBack,
    OpenBack,
    consistent_radiation!,
    incoming_radiation,
    exposed_face,
    find_kernel,
    reset!,
    contribute!,
    DIFFUSIVITY
using Pyrolysis.Assembly
using Pyrolysis.Assembly: residual!, buffers
using Pyrolysis.ColoredJacobian
using Pyrolysis.Integration
using Pyrolysis.Integration: setup

include(joinpath(@__DIR__, "..", "verification", "exact", "two_stream.jl"))

# A semi-transparent decomposing slab: the virgin solid's emissivity and
# density depend on temperature, its absorption coefficient differs from the
# char's, so the radiation rows read the temperature and the composition of
# their cells. No component sets an absorptivity: diffuse radiation enters at
# the emissivity.
function _slab(
    n;
    top,
    bottom = AdiabaticBC(),
    λ = 0.0,
    reacting = true,
    coordinates = Planar(),
    mesh = stretched_mesh(0.01, n, 1.15; coordinates),
)
    virgin = SolidComponent(
        :virgin;
        ρ_intrinsic = 1600.0,
        λ = λ,
        ρ = T -> 1100.0 - 0.05 * (T - 300.0),
        c = (900.0, 1.5),
        k = (0.2, 2e-4),
        ε = (0.95, -1.5e-4),
        α = 1.5,
    )
    char = SolidComponent(
        :char;
        ρ_intrinsic = 1800.0,
        λ = λ,
        ρ = 300.0,
        c = 1100.0,
        k = 0.1,
        ε = 0.85,
        α = 20.0,
    )
    gas = GaseousComponent(:gas, M = 0.03, c = 1100.0, k = 0.03)
    mat = Material(
        name = :translucent,
        components = (virgin, char, gas),
        reactions = reacting ?
                    (
            Reaction(:pyro, 1 => (2, 3); A = 5e9, E = 1.4e5, h = 8e5, yields = (0.3, 0.7)),
        ) : (),
        depletion_limiter = DepletionLimiter(threshold = 1.0),
    )
    return PyrolysisProblem(
        mesh = mesh,
        material = mat,
        bc_set = BoundaryConditionSet(thermal = (top = top, bottom = bottom)),
        T_initial = z -> 400.0 + 4e4 * z,
        ξ_initial = [z -> 1100.0 * (1 - 30z), z -> 300.0 * 30z, z -> 0.5],
        tspan = (0.0, 10.0),
    )
end

_flux(t) = 3e4 * (1 + 0.5 * sin(t))
const _tops = (
    heated = RadiativeFluxBC(_flux) +
             RadiativeBC(T_inf = 300.0) +
             ConvectiveBC(h = 10.0, T_inf = 300.0),
    surroundings = RadiativeBC(T_inf = t -> 900.0 + 20t, F = 0.8) +
                   ConvectiveBC(h = 10.0, T_inf = 300.0),
)
const _backs = (
    mirror = AdiabaticBC(),
    open = RadiativeBC(T_inf = 350.0) +
           RadiativeFluxBC(2e3) +
           ConvectiveBC(h = 5.0, T_inf = 300.0),
)

function _states(cp)
    u0 = cp.u0
    layout = cp.def.layout
    hot = copy(u0)
    for i = 1:n_cells(layout)
        hot[energy_index(layout, i)] =
            energy_state(layout) isa TemperatureState ? 700.0 + 20i :
            u0[energy_index(layout, i)] * 3
    end
    perturbed = u0 .* (1 .+ 0.2 .* sin.(1:length(u0)))
    return (; u0, hot, perturbed, zero_state = zeros(length(u0)), negative = -abs.(u0))
end

# The radiation rows of every face, and the tolerance scale of the rows.
_radiation_residuals(du, layout) = [
    du[radiation_index(layout, f, k)] for f = 1:n_nodes(layout) for
    k = 1:radiation_rows(layout)
]
_q_scale(cp) = cp.abstol[radiation_index(cp.def.layout, 1, 1)] / cp.reltol

# The faces from the exposed face into the depth, and the cells' thicknesses
# on a fixed line.
_faces(kernel) = [kernel.up[first(kernel.order)]; [kernel.down[i] for i in kernel.order]]
function _thickness(ds, i)
    x = ds.geometry.node_position
    return x[i+1][1] - x[i][1]
end

@testset "Radiative transfer: $(energy) form, $(n) cells, $(top) top, $(back) back, gas λ = $(λ)" for energy in
                                                                                                      (
        :temperature,
        :enthalpy,
    ),
    n in (6, 8),
    top in keys(_tops),
    back in keys(_backs),
    λ in (0.0, 1e-5)

    cp = setup(
        _slab(n; top = _tops[top], bottom = _backs[back], λ);
        energy_state = energy,
        radiation_model = RADIATIVE_TRANSFER,
    )
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    layout = def.layout
    kernel = find_kernel(RadiativeTransferKernel, def.kernels)
    @test kernel.back isa (back === :mirror ? MirrorBack : OpenBack)
    @test has_radiation(layout) && radiation_rows(layout) == 2
    @test count(!, is_differential(layout)) == 2 * (n + 1)
    pattern = cj.pattern
    @test all(pattern[r, r] for r in axes(pattern, 1))

    du = zeros(length(cp.u0))
    residual!(du, cp.u0, def, 0.0, ws)
    @test all(isfinite, du)
    # The initial fluxes satisfy their rows.
    @test maximum(abs, _radiation_residuals(du, layout)) <= 1e-12 * _q_scale(cp)
    residual!(du, cp.u0, def, 0.0, ws)
    @test (@allocated residual!(du, cp.u0, def, 0.0, ws)) == 0

    J = jac_prototype(cj)
    for (state, u) in pairs(_states(cp))
        @testset "$state" begin
            residual!(du, u, def, 1.0, ws)
            @test all(isfinite, du)
            Jd = dense_jacobian(def, ws, u, 1.0)
            @test all(isfinite, Jd)
            @test isempty(missing_entries(Jd, pattern))
            jacobian!(J, cj, u, 1.0)
            @test compare_entries(Matrix(J), Jd; rtol = 1e-12).passed
        end
    end
end

# The rows solve the two-stream equations of the frozen cell properties: the
# consistent fluxes are the exact solution of `exact/two_stream.jl`, the rows
# at arbitrary trial fluxes are its residuals, and the cell sources its
# sources, with the first live cell at the exposed boundary or below it.
@testset "Radiative transfer: rows against the exact two-stream solution, $(back) back, first live cell $(s)" for back in
                                                                                                                  keys(
        _backs,
    ),
    s in (1, 3)

    n = 7
    top, bottom = _tops.heated, _backs[back]
    cp =
        setup(_slab(n; top, bottom, reacting = false); radiation_model = RADIATIVE_TRANSFER)
    def, ws = cp.def, cp.ws
    layout, topo, material = def.layout, def.setup.topology, def.setup.material
    kernel = find_kernel(RadiativeTransferKernel, def.kernels)
    u = copy(cp.u0)
    for i = 1:n
        u[energy_index(layout, i)] = 500.0 + 40i
    end
    ds, c = buffers(ws, def, Float64)
    derive!(ds, u, layout, def.setup)
    for k = 1:(s-1)
        ds.active[kernel.order[k]] = false
    end
    t = 1.3
    cells, faces = kernel.order, _faces(kernel)
    @test exposed_face(kernel, ds) == faces[s]
    τ = [ds.α[i] * _thickness(ds, i) for i in cells]
    B = [STEFAN_BOLTZMANN * ds.T_trial[i]^4 for i in cells]
    inputs = (;
        D = DIFFUSIVITY,
        ε = ds.ε[cells[s]],
        Q = incoming_radiation(top, t),
        ε_b = back === :open ? ds.ε[last(cells)] : 0.0,
        Q_b = back === :open ? incoming_radiation(bottom, t) : 0.0,
        first = s,
    )
    @test 0.1 < minimum(τ[s:end]) && maximum(τ) < 5
    exact = two_stream_solve(τ, B; inputs...)
    scale = max(maximum(exact.qm), maximum(exact.qp))

    q = fill(NaN, 2, n + 1)
    consistent_radiation!(q, kernel, ds, topo, t)
    @test maximum(abs, q[1, faces] .- exact.qm) <= 1e-12 * scale
    @test maximum(abs, q[2, faces] .- exact.qp) <= 1e-12 * scale
    @test all(iszero, q[:, faces[1:(s-1)]])

    # Arbitrary trial fluxes: the rows are the exact solution's residuals.
    qm = scale .* (0.5 .+ 0.4 .* sin.(1:(n+1)))
    qp = scale .* (0.5 .+ 0.4 .* cos.(1:(n+1)))
    ds.q[1, faces] .= qm
    ds.q[2, faces] .= qp
    reset!(c)
    contribute!(c, kernel, ds, topo, material, t)
    none = zeros(n + 1, 0)
    rows = two_stream_residual(qm, qp, none, none, τ, B; inputs...)
    live = faces[s:end]
    @test maximum(abs, c.r_radiation[1, live] .- rows.thermal[1, s:end]) <= 1e-12 * scale
    @test maximum(abs, c.r_radiation[2, live] .- rows.thermal[2, s:end]) <= 1e-12 * scale
    S = two_stream_sources(qm, qp, none, none; first = s)
    A = ds.geometry.face_area
    @test maximum(abs, [c.S_energy[i] / A[1] for i in cells] .- S) <= 1e-12 * scale
end

# An isothermal medium facing diffuse radiation at its own temperature on
# both faces, or on its exposed face with a mirror behind it, is in
# equilibrium at every optical thickness and emissivity: q⁻ = q⁺ = B and no
# cell gains or loses energy.
@testset "Radiative transfer: an isothermal enclosure is in equilibrium, $(back) back" for back in
                                                                                           (
    :mirror,
    :open,
)

    n, T_h = 12, 700.0
    bottom = back === :mirror ? AdiabaticBC() : RadiativeBC(T_inf = T_h)
    problem = _slab(n; top = RadiativeBC(T_inf = T_h), bottom, reacting = false)
    cp = setup(problem; radiation_model = RADIATIVE_TRANSFER)
    def, ws = cp.def, cp.ws
    layout, topo, material = def.layout, def.setup.topology, def.setup.material
    kernel = find_kernel(RadiativeTransferKernel, def.kernels)
    u = copy(cp.u0)
    for i = 1:n
        u[energy_index(layout, i)] = T_h
    end
    ds, c = buffers(ws, def, Float64)
    derive!(ds, u, layout, def.setup)
    @test all(==(T_h), ds.T_trial)
    B = STEFAN_BOLTZMANN * T_h^4
    Δ = [_thickness(ds, i) for i = 1:n]
    q = zeros(2, n + 1)
    worst_S, worst_q = 0.0, 0.0
    for ε in (0.3, 0.96, 1.0), κΔ in (1e-6, 1e-3, 0.1, 1.0, 10.0, 50.0)
        ds.ε .= ε
        ds.α .= κΔ / (sum(Δ) / n)
        consistent_radiation!(q, kernel, ds, topo, 0.0)
        ds.q .= q
        reset!(c)
        contribute!(c, kernel, ds, topo, material, 0.0)
        worst_S = max(worst_S, maximum(abs, c.S_energy) / ds.geometry.face_area[1])
        worst_q = max(worst_q, maximum(abs, q .- B))
        @test maximum(abs, c.r_radiation) <= 1e-12 * B
    end
    @test worst_S <= 1e-12 * B
    @test worst_q <= 1e-12 * B
end

# The cell sources telescope to the net flux through the exposed face less
# the net flux through the back, at any trial fluxes, on a moving mesh too.
@testset "Radiative transfer: the cell sources telescope, $(back) back, moving = $(moving)" for back in
                                                                                                keys(
        _backs,
    ),
    moving in (false, true)

    n = 8
    cp = setup(
        _slab(n; top = _tops.heated, bottom = _backs[back]);
        use_ale = moving,
        radiation_model = RADIATIVE_TRANSFER,
    )
    def, ws = cp.def, cp.ws
    layout, topo, material = def.layout, def.setup.topology, def.setup.material
    kernel = find_kernel(RadiativeTransferKernel, def.kernels)
    ds, c = buffers(ws, def, Float64)
    for s in (1, 4)
        derive!(ds, cp.u0, layout, def.setup)
        for k = 1:(s-1)
            ds.active[kernel.order[k]] = false
        end
        ds.q .= 2e4 .* (1 .+ 0.7 .* sin.(reshape(1:length(ds.q), size(ds.q))))
        reset!(c)
        contribute!(c, kernel, ds, topo, material, 2.0)
        A = ds.geometry.face_area
        e, b = exposed_face(kernel, ds), kernel.down[last(kernel.order)]
        net = A[e] * (ds.q[1, e] - ds.q[2, e]) - A[b] * (ds.q[1, b] - ds.q[2, b])
        @test sum(c.S_energy) ≈ net rtol = 1e-12
        @test all(iszero, c.S_energy[kernel.order[1:(s-1)]])
    end
end

# The radiation rows read the first live cell's emissivity, its temperature
# dependence included, and every live cell's temperature and absorption
# coefficient: the colored Jacobian's columns of the exposed cell's
# temperature and virgin mass, and on a moving mesh of the exposed node's
# position, against central differences of the residual, an oracle
# independent of automatic differentiation.
@testset "Radiative transfer: exposed-cell columns against central differences, moving = $(moving)" for moving in
                                                                                                        (
    false,
    true,
)

    n = 7
    cp = setup(
        _slab(n; top = _tops.heated, bottom = _backs.open);
        use_ale = moving,
        radiation_model = RADIATIVE_TRANSFER,
    )
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    layout = def.layout
    J = jac_prototype(cj)
    jacobian!(J, cj, cp.u0, 1.0)
    row_floor = vec(maximum(abs, J; dims = 2))
    columns = [(energy_index(layout, n), 1e-3), (mass_index(layout, n, 1), 1e-6)]
    moving && push!(columns, (position_index(layout, n + 1), 1e-9))
    for (k, h) in columns
        up, down = copy(cp.u0), copy(cp.u0)
        up[k] += h
        down[k] -= h
        r_up, r_down = zeros(length(up)), zeros(length(up))
        residual!(r_up, up, def, 1.0, ws)
        residual!(r_down, down, def, 1.0, ws)
        fd = (r_up .- r_down) ./ 2h
        col = Vector(J[:, k])
        @test all(abs.(col .- fd) .<= 1e-6 .* max.(abs.(fd), row_floor))
    end
    # The exposed row reads the exposed cell's emissivity, which depends on
    # its composition and temperature; the q⁺ row of the exposed face reads
    # the cell's emission.
    e_down, e_up = radiation_index(layout, n + 1, 1), radiation_index(layout, n + 1, 2)
    T, m = energy_index(layout, n), mass_index(layout, n, 1)
    @test J[e_down, T] != 0 && J[e_down, m] != 0
    @test J[e_up, T] != 0 && J[e_up, m] != 0
end

@testset "Radiative transfer on a moving mesh: $(energy) form, $(top) top" for energy in (
        :temperature,
        :enthalpy,
    ),
    top in keys(_tops)

    cp = setup(
        _slab(6; top = _tops[top], bottom = _backs.open);
        energy_state = energy,
        use_ale = true,
        radiation_model = RADIATIVE_TRANSFER,
    )
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    layout = def.layout
    @test is_moving(layout) && radiation_rows(layout) == 2
    du = zeros(length(cp.u0))
    residual!(du, cp.u0, def, 0.0, ws)
    @test maximum(abs, _radiation_residuals(du, layout)) <= 1e-12 * _q_scale(cp)
    J = jac_prototype(cj)
    for (state, u) in pairs(_states(cp))
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

# The rows carry fluxes per unit face area and the sources read the face
# areas, which agree only on planar lines.
@testset "Radiative transfer is refused off planar coordinates" begin
    for coordinates in (Cylindrical(), Spherical())
        problem = _slab(6; top = _tops.heated, coordinates)
        @test_throws ArgumentError setup(problem; radiation_model = RADIATIVE_TRANSFER)
        cp = setup(problem; radiation_model = SURFACE_ABSORPTION)
        @test !has_radiation(cp.def.layout)
    end
end
