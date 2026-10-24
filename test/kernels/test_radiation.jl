using Test
using ForwardDiff
using SparseArrays
using Pyrolysis
using Pyrolysis.Materials: DepletionLimiter, BEER_LAMBERT, SURFACE_ABSORPTION
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.States: StateLayout, state_length
using Pyrolysis.Kernels: BeerLambertKernel, incident_intensity, transmitted_face
using Pyrolysis.Assembly
using Pyrolysis.Assembly: residual!, buffers
using Pyrolysis.ColoredJacobian
using Pyrolysis.Integration
using Pyrolysis.Integration: setup

# The surface radiative properties of the slab's solids. With
# `component_absorptivity` the virgin solid's emissivity depends on
# temperature but it sets its own absorptivity, and the char absorbs at its
# emissivity: the exposed cell's absorptivity follows its composition alone.
# With `absorptivity_of_temperature` the virgin solid absorbs at its
# temperature-dependent emissivity (Kirchhoff's law) and its density depends
# on temperature, while the char sets its own absorptivity: the exposed cell's
# absorptivity also depends on its temperature.
const _absorbers = (
    component_absorptivity = (
        virgin = (ρ = 1100.0, ε = (0.8, 5e-5), absorptivity = 0.7),
        char = (ε = 0.95,),
    ),
    absorptivity_of_temperature = (
        virgin = (ρ = T -> 1100.0 - 0.05 * (T - 300.0), ε = (0.95, -1.5e-4)),
        char = (ε = 0.85, absorptivity = 0.6),
    ),
)

# A semi-transparent decomposing slab: temperature-dependent emissivities,
# a component absorptivity different from its emissivity, mass absorption
# coefficients that change as the virgin solid turns to char, and a
# time-dependent incident flux.
function _translucent_slab(
    n;
    top,
    gas_λ = 0.0,
    bottom = AdiabaticBC(),
    reacting = true,
    coordinates = Planar(),
    absorber = :component_absorptivity,
)
    virgin = SolidComponent(
        :virgin;
        c = (900.0, 1.5),
        k = (0.2, 2e-4),
        α = 1.5,
        _absorbers[absorber].virgin...,
    )
    char = SolidComponent(
        :char;
        ρ = 300.0,
        c = 1100.0,
        k = 0.1,
        α = 20.0,
        _absorbers[absorber].char...,
    )
    gas = GaseousComponent(:gas, M = 0.03, c = 1100.0, k = 0.03, λ = gas_λ)
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
        mesh = create_stretched_mesh(0.01, n, 1.15; coordinates),
        material = mat,
        bc_set = BoundaryConditionSet(thermal = (top = top, bottom = bottom)),
        T_initial = z -> 400.0 + 4e4 * z,
        ξ_initial = [z -> 1100.0 * (1 - 30z), z -> 300.0 * 30z, z -> 0.5],
        tspan = (0.0, 10.0),
    )
end

_flux(t) = 3e4 * (1 + 0.5 * sin(t))
const _radiation_tops = (
    kirchhoff = RadiativeFluxBC(flux = _flux) +
                RadiativeBC(T_inf = 300.0) +
                ConvectiveBC(h = 10.0, T_inf = 300.0),
    explicit = RadiativeFluxBC(flux = 4e4, absorptivity = t -> 0.9 - 0.01t) +
               ConvectiveBC(h = 10.0, T_inf = 300.0),
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

@testset "Beer–Lambert: $(energy) form, $(n) cells, $(name), $(absorber), gas λ = $(λ)" for energy in
                                                                                            (
        :temperature,
        :enthalpy,
    ),
    n in (6, 8),
    name in keys(_radiation_tops),
    absorber in keys(_absorbers),
    λ in (0.0, 1e-5)

    cp = setup(
        _translucent_slab(n; top = _radiation_tops[name], gas_λ = λ, absorber);
        energy_state = energy,
        radiation_model = BEER_LAMBERT,
    )
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    layout = def.layout
    @test has_intensities(layout)
    @test count(!, is_differential(layout)) == n + 1
    pattern = cj.pattern
    @test all(pattern[r, r] for r in axes(pattern, 1))
    @test cj.n_colors <= 3 * (n_components(layout) + 1) + 2

    du = zeros(length(cp.u0))
    residual!(du, cp.u0, def, 0.0, ws)
    @test all(isfinite, du)
    # The initial intensities satisfy their rows.
    I_rows = [du[intensity_index(layout, f)] for f = 1:(n+1)]
    @test maximum(abs, I_rows) <= 1e-12 * 4e4
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

@testset "Beer–Lambert: the beam's energy telescopes, Kirchhoff and explicit absorptivity" begin
    for name in keys(_radiation_tops)
        cp = setup(
            _translucent_slab(7; top = _radiation_tops[name], reacting = false);
            radiation_model = BEER_LAMBERT,
        )
        def, ws = cp.def, cp.ws
        kernel = only(k for k in def.kernels if k isa BeerLambertKernel)
        du = zeros(length(cp.u0))
        t = 2.0
        Integration._consistent_intensities!(cp.u0, def, ws, t)
        residual!(du, cp.u0, def, t, ws)
        ds, c = buffers(ws, def, Float64)
        a = name === :kirchhoff ? ds.a[7] : 0.9 - 0.01t
        q = name === :kirchhoff ? _flux(t) : 4e4
        I_in = incident_intensity(kernel, ds, t)
        @test I_in ≈ a * q rtol = 1e-14
        @test ds.I[kernel.exposed] ≈ I_in rtol = 1e-14
        # Exact cell averages: each cell absorbs I_up (1 − e^{−τ}), and the
        # sum is the entering intensity less the transmitted one.
        I_up = I_in
        for i = 7:-1:1
            g = ds.geometry
            Δ = g.node_position[i+1][1] - g.node_position[i][1]
            τ = ds.α[i] * Δ
            @test c.S_energy[i] ≈ I_up * (1 - exp(-τ)) rtol = 1e-12
            I_up *= exp(-τ)
        end
        @test ds.I[transmitted_face(kernel)] ≈ I_up rtol = 1e-12
        absorbed = sum(c.S_energy[i] for i = 1:7)
        @test absorbed + ds.I[1] ≈ I_in rtol = 1e-12
    end
end

# The Kirchhoff and ε(T) chains reach the energy rows: the exposed cell's
# temperature and virgin-mass columns of the colored Jacobian against central
# differences of the residual, an oracle independent of automatic
# differentiation. The exposed cell's absorptivity follows its composition,
# so under Kirchhoff's law the mass column reaches the exposed face's
# intensity row, and under an explicit absorptivity it does not. The
# temperature column reaches that row only when the cell's absorptivity
# depends on its temperature and the condition sets no absorptivity of its
# own.
@testset "Beer–Lambert: exposed-cell columns against central differences" begin
    for absorber in keys(_absorbers), name in keys(_radiation_tops)
        cp = setup(
            _translucent_slab(7; top = _radiation_tops[name], absorber);
            radiation_model = BEER_LAMBERT,
        )
        def, ws, cj = cp.def, cp.ws, cp.jacobian
        layout = def.layout
        J = jac_prototype(cj)
        jacobian!(J, cj, cp.u0, 1.0)
        # Per entry, with the row's largest Jacobian entry as its floor.
        row_floor = vec(maximum(abs, J; dims = 2))
        for (k, h) in ((energy_index(layout, 7), 1e-3), (mass_index(layout, 7, 1), 1e-6))
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
        e = intensity_index(layout, 8)
        m = mass_index(layout, 7, 1)
        name === :kirchhoff ? (@test J[e, m] != 0) : (@test J[e, m] == 0)
        c = energy_index(layout, 7)
        if name === :kirchhoff && absorber === :absorptivity_of_temperature
            @test J[e, c] != 0
        else
            @test abs(J[e, c]) <= 1e-12 * maximum(abs, J[e, :])
        end
    end
end

# On a moving mesh the exposed cell's absorptivity, its temperature
# dependence included, reaches the exposed face's intensity row as on a fixed
# mesh, and the pattern contains the Jacobian with the position columns.
@testset "Beer–Lambert on a moving mesh: $(energy) form, $(name), absorptivity of temperature" for energy in
                                                                                                   (
        :temperature,
        :enthalpy,
    ),
    name in keys(_radiation_tops)

    cp = setup(
        _translucent_slab(
            6;
            top = _radiation_tops[name],
            absorber = :absorptivity_of_temperature,
        );
        energy_state = energy,
        use_ale = true,
        radiation_model = BEER_LAMBERT,
    )
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    layout = def.layout
    @test is_moving(layout) && has_intensities(layout)
    du = zeros(length(cp.u0))
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
            if state === :u0
                e, c = intensity_index(layout, 7), energy_index(layout, 6)
                if name === :kirchhoff
                    @test Jd[e, c] != 0
                else
                    @test abs(Jd[e, c]) <= 1e-12 * maximum(abs, Jd[e, :])
                end
            end
        end
    end
end

# On a cylinder or a sphere the intensity rows attenuate the intensity while
# the cell sources read the face areas, so a transparent cell would absorb
# `I (A_up − A_down)`: setup refuses the combination.
@testset "Beer–Lambert is refused off planar coordinates" begin
    for coordinates in (Cylindrical(), Spherical())
        problem = _translucent_slab(6; top = _radiation_tops.kirchhoff, coordinates)
        @test_throws ArgumentError setup(problem; radiation_model = BEER_LAMBERT)
        cp = setup(problem; radiation_model = SURFACE_ABSORPTION)
        @test !has_intensities(cp.def.layout)
    end
end

@testset "Beer–Lambert requires a single incident source" begin
    first_source = RadiativeFluxBC(5e4; absorptivity = 0.9)
    second_source = RadiativeFluxBC(2e4; absorptivity = 0.9)
    for top in (
        first_source + second_source,
        CombinedThermalBC((
            ConvectiveBC(h = 10.0, T_inf = 300.0),
            CombinedThermalBC((first_source, second_source)),
        )),
    )
        problem = _translucent_slab(6; top, reacting = false)
        @test_throws r"supports one RadiativeFluxBC" setup(
            problem;
            radiation_model = BEER_LAMBERT,
        )
        cp = setup(problem; radiation_model = SURFACE_ABSORPTION)
        @test !has_intensities(cp.def.layout)
    end
    problem = _translucent_slab(
        6;
        top = RadiativeFluxBC(t -> 5e4 + 2e4; absorptivity = 0.9),
        reacting = false,
    )
    cp = setup(problem; radiation_model = BEER_LAMBERT)
    kernel = only(k for k in cp.def.kernels if k isa BeerLambertKernel)
    @test incident_intensity(kernel, first(buffers(cp.ws, cp.def, Float64)), 0.0) ≈ 6.3e4
end
