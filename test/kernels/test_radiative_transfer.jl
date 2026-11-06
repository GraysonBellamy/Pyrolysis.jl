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

include(joinpath(@__DIR__, "..", "verification", "exact", "gauss_legendre.jl"))
include(joinpath(@__DIR__, "..", "verification", "exact", "two_stream.jl"))
include(joinpath(@__DIR__, "..", "verification", "exact", "gray_slab.jl"))

# A semi-transparent decomposing slab: the virgin solid's emissivity and
# density depend on temperature, its absorption coefficient differs from the
# char's, so the radiation rows read the temperature and the composition of
# their cells. Unless `absorbing`, no component sets an absorptivity and
# every radiation enters at the emissivity; with it the virgin solid absorbs
# beams at 0.7, which no diffuse flux may accompany.
function _slab(
    n;
    top,
    bottom = AdiabaticBC(),
    λ = 0.0,
    reacting = true,
    absorbing = false,
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
        absorptivity = absorbing ? 0.7 : nothing,
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
_beam_absorptivity(t) = 0.85 + 0.05 * sin(t)
_side_flux(t) = 5e3 * (1 + 0.1t)
const _cone = Cone(deg2rad(20), deg2rad(73))
# Diffuse incidence, surroundings only, a collimated heater entering at its own
# absorptivity, and a cone beside a collimated beam entering at the exposed
# cell's (their beam tuple mixes a constant and a time-varying flux).
const _tops = (
    heated = RadiativeFluxBC(_flux) +
             RadiativeBC(T_inf = 300.0) +
             ConvectiveBC(h = 10.0, T_inf = 300.0),
    surroundings = RadiativeBC(T_inf = t -> 900.0 + 20t, F = 0.8) +
                   ConvectiveBC(h = 10.0, T_inf = 300.0),
    collimated = RadiativeFluxBC(
                     _flux;
                     absorptivity = _beam_absorptivity,
                     distribution = Collimated(),
                 ) +
                 RadiativeBC(T_inf = 300.0) +
                 ConvectiveBC(h = 10.0, T_inf = 300.0),
    cone = RadiativeFluxBC(2.5e4; distribution = _cone) +
           RadiativeFluxBC(_side_flux; distribution = Collimated()) +
           RadiativeBC(T_inf = 300.0) +
           ConvectiveBC(h = 10.0, T_inf = 300.0),
)
# The radiation rows per face: q⁻, q⁺ and one row per beam.
const _rows = (heated = 2, surroundings = 2, collimated = 3, cone = 6)
# The beam tops run on the slab whose virgin solid sets its absorptivity
# where no diffuse flux arrives: behind a mirror, not the open back.
_absorbing(top, back) = top in (:collimated, :cone) && back === :mirror
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
        _slab(
            n;
            top = _tops[top],
            bottom = _backs[back],
            λ,
            absorbing = _absorbing(top, back),
        );
        energy_state = energy,
        radiation_model = RADIATIVE_TRANSFER,
    )
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    layout = def.layout
    kernel = find_kernel(RadiativeTransferKernel, def.kernels)
    @test kernel.back isa (back === :mirror ? MirrorBack : OpenBack)
    @test has_radiation(layout) && radiation_rows(layout) == _rows[top]
    @test count(!, is_differential(layout)) == _rows[top] * (n + 1)
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
# sources, with the first live cell at the exposed boundary or below it. The
# beams enter at the condition's absorptivity or, without one, at the first
# live cell's, which differs from its emissivity on the absorbing slab; the
# beams reaching the back join q⁺ there.
@testset "Radiative transfer: rows against the exact two-stream solution, $(top) top, $(back) back, first live cell $(s)" for top in
                                                                                                                              (
        :heated,
        :collimated,
        :cone,
    ),
    back in keys(_backs),
    s in (1, 3)

    n = 7
    bc, bottom = _tops[top], _backs[back]
    cp = setup(
        _slab(n; top = bc, bottom, reacting = false, absorbing = _absorbing(top, back));
        radiation_model = RADIATIVE_TRANSFER,
    )
    def, ws = cp.def, cp.ws
    layout, topo, material = def.layout, def.setup.topology, def.setup.material
    kernel = find_kernel(RadiativeTransferKernel, def.kernels)
    R = radiation_rows(layout)
    @test R == _rows[top]
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
    a = ds.a[cells[s]]
    beams =
        top === :collimated ? [TwoStreamBeam(1.0, _beam_absorptivity(t) * _flux(t))] :
        top === :cone ?
        [
            two_stream_cone_beams(2.5e4, deg2rad(20), deg2rad(73); a)
            TwoStreamBeam(1.0, a * _side_flux(t))
        ] : TwoStreamBeam[]
    _absorbing(top, back) && @test abs(a - ds.ε[cells[s]]) > 0.05
    inputs = (;
        D = DIFFUSIVITY,
        ε = ds.ε[cells[s]],
        Q = incoming_radiation(bc, t),
        ε_b = back === :open ? ds.ε[last(cells)] : 0.0,
        Q_b = back === :open ? incoming_radiation(bottom, t) : 0.0,
        beams,
        first = s,
    )
    @test 0.1 < minimum(τ[s:end]) && maximum(τ) < 5
    exact = two_stream_solve(τ, B; inputs...)
    scale = max(maximum(exact.qm), maximum(exact.qp), maximum(exact.Id; init = 0.0))
    nb = length(beams)
    @test R == 2 + nb

    q = fill(NaN, R, n + 1)
    consistent_radiation!(q, kernel, ds, topo, t)
    @test maximum(abs, q[1, faces] .- exact.qm) <= 1e-12 * scale
    @test maximum(abs, q[2, faces] .- exact.qp) <= 1e-12 * scale
    for h = 1:nb
        @test maximum(abs, q[2+h, faces] .- exact.Id[:, h]) <= 1e-12 * scale
    end
    @test all(iszero, q[:, faces[1:(s-1)]])

    # Arbitrary trial fluxes: the rows are the exact solution's residuals.
    qm = scale .* (0.5 .+ 0.4 .* sin.(1:(n+1)))
    qp = scale .* (0.5 .+ 0.4 .* cos.(1:(n+1)))
    Id = scale .* (0.3 .+ 0.2 .* sin.((1:(n+1)) .+ (1:nb)'))
    ds.q[1, faces] .= qm
    ds.q[2, faces] .= qp
    for h = 1:nb
        ds.q[2+h, faces] .= Id[:, h]
    end
    reset!(c)
    contribute!(c, kernel, ds, topo, material, t)
    Iu = zeros(n + 1, nb)
    rows = two_stream_residual(qm, qp, Id, Iu, τ, B; inputs...)
    live = faces[s:end]
    @test maximum(abs, c.r_radiation[1, live] .- rows.thermal[1, s:end]) <= 1e-12 * scale
    @test maximum(abs, c.r_radiation[2, live] .- rows.thermal[2, s:end]) <= 1e-12 * scale
    for h = 1:nb
        @test maximum(abs, c.r_radiation[2+h, live] .- rows.down[s:end, h]) <= 1e-12 * scale
    end
    S = two_stream_sources(qm, qp, Id, Iu; first = s)
    A = ds.geometry.face_area
    @test maximum(abs, [c.S_energy[i] / A[1] for i in cells] .- S) <= 1e-12 * scale
end

# A cone heater travels as three beams, the Gauss–Legendre nodes of its
# direction cosine. Against the exact transmission of the cone through a cold
# slab, the flux the beams carry to each face and the power each cell absorbs
# from them are within 2×10⁻⁴ of the flux entering.
@testset "Radiative transfer: a cone heater against the exact transmission" begin
    n = 12
    cp = setup(
        _slab(
            n;
            top = RadiativeFluxBC(2.5e4; distribution = _cone) + RadiativeBC(T_inf = 300.0),
            reacting = false,
        );
        radiation_model = RADIATIVE_TRANSFER,
    )
    def, ws = cp.def, cp.ws
    layout, topo = def.layout, def.setup.topology
    kernel = find_kernel(RadiativeTransferKernel, def.kernels)
    @test radiation_rows(layout) == 5
    ds, _ = buffers(ws, def, Float64)
    derive!(ds, cp.u0, layout, def.setup)
    cells, faces = kernel.order, _faces(kernel)
    q = zeros(5, n + 1)
    consistent_radiation!(q, kernel, ds, topo, 0.0)
    entering = ds.a[cells[1]] * 2.5e4
    depth = [0.0; cumsum([ds.α[i] * _thickness(ds, i) for i in cells])]
    @test depth[end] > 3
    exact =
        [gray_slab_cone_transmission(entering, d, deg2rad(20), deg2rad(73)) for d in depth]
    carried = vec(sum(q[3:5, faces]; dims = 1))
    @test carried[1] ≈ entering rtol = 1e-14
    @test maximum(abs, carried .- exact) <= 2e-4 * entering
    @test maximum(abs, diff(carried) .- diff(exact)) <= 2e-4 * entering
end

# The beams of a cone carry shares of its flux proportional to their weight
# times their direction cosine, which lie inside the cone: the nodes and
# shares of `exact/two_stream.jl`'s cone. A collimated heater is one beam
# along the normal.
@testset "Radiative transfer: the directions of a heater's beams" begin
    for (θ₁, θ₂) in ((0.0, π / 2), (deg2rad(20), deg2rad(73)), (0.3, 0.31)), d = 1:5
        cone = Cone(θ₁, θ₂; directions = d)
        cp = setup(
            _slab(4; top = RadiativeFluxBC(1e4; distribution = cone), reacting = false);
            radiation_model = RADIATIVE_TRANSFER,
        )
        source = only(find_kernel(RadiativeTransferKernel, cp.def.kernels).beams)
        reference = two_stream_cone_beams(1.0, θ₁, θ₂; directions = d)
        @test radiation_rows(cp.def.layout) == 2 + d
        @test source.row == 3
        @test maximum(abs, source.μ .- [b.μ for b in reference]) <= 1e-14
        @test maximum(abs, source.w .- [b.q0 for b in reference]) <= 1e-14
        @test sum(source.w) ≈ 1 rtol = 1e-14
        @test all(cos(θ₂) .< source.μ .< cos(θ₁))
    end
    cp = setup(
        _slab(4; top = RadiativeFluxBC(1e4; distribution = Collimated()), reacting = false);
        radiation_model = RADIATIVE_TRANSFER,
    )
    source = only(find_kernel(RadiativeTransferKernel, cp.def.kernels).beams)
    @test source.μ == [1.0] && source.w == [1.0] && source.row == 3
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
# the net flux through the back, beams included, at any trial fluxes, on a
# moving mesh too.
@testset "Radiative transfer: the cell sources telescope, $(top) top, $(back) back, moving = $(moving)" for top in
                                                                                                            (
        :heated,
        :cone,
    ),
    back in keys(_backs),
    moving in (false, true)

    n = 8
    cp = setup(
        _slab(
            n;
            top = _tops[top],
            bottom = _backs[back],
            absorbing = _absorbing(top, back),
        );
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
        beams(f) = sum(ds.q[3:end, f]; init = 0.0)
        net =
            A[e] * (ds.q[1, e] - ds.q[2, e] + beams(e)) -
            A[b] * (ds.q[1, b] + beams(b) - ds.q[2, b])
        @test size(ds.q, 1) == _rows[top]
        @test sum(c.S_energy) ≈ net rtol = 1e-12
        @test all(iszero, c.S_energy[kernel.order[1:(s-1)]])
    end
end

# The radiation rows read the first live cell's emissivity, its temperature
# dependence included, the beams its absorptivity (behind a mirror on the
# slab that sets one), and every live cell's
# temperature and absorption coefficient: the colored Jacobian's columns of
# the exposed cell's temperature and virgin mass, and on a moving mesh of the
# exposed node's position, against central differences of the residual, an
# oracle independent of automatic differentiation.
@testset "Radiative transfer: exposed-cell columns against central differences, $(top) top, moving = $(moving)" for top in
                                                                                                                    (
        :heated,
        :cone,
    ),
    moving in (false, true)

    n = 7
    back = top === :cone ? :mirror : :open
    cp = setup(
        _slab(
            n;
            top = _tops[top],
            bottom = _backs[back],
            absorbing = _absorbing(top, back),
        );
        use_ale = moving,
        radiation_model = RADIATIVE_TRANSFER,
    )
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    layout = def.layout
    J = jac_prototype(cj)
    jacobian!(J, cj, cp.u0, 1.0)
    row_floor = vec(maximum(abs, J; dims = 2))
    columns = [(energy_index(layout, n), 1e-3), (mass_index(layout, n, 1), 1e-6)]
    moving && push!(columns, (position_index(layout, n + 1), 1e-8))
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
    # A beam without its own absorptivity enters at the exposed cell's, which
    # depends on its composition.
    if top === :cone
        @test all(J[radiation_index(layout, n + 1, k), m] != 0 for k = 3:6)
    end
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
    @test is_moving(layout) && radiation_rows(layout) == _rows[top]
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
