using Test
using ForwardDiff
using SparseArrays
using Pyrolysis
using Pyrolysis.Materials: DepletionLimiter, BEER_LAMBERT, SURFACE_ABSORPTION
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.States: StateLayout, state_length
using Pyrolysis.Assembly
using Pyrolysis.Assembly: residual!
using Pyrolysis.ColoredJacobian
using Pyrolysis.Integration: setup

# ============================================================================
# Materials and columns
# ============================================================================

# The PMMA-like and wood-like materials of `test/perf/perf_helpers.jl`.
function _pmma(; gas_c = 1500.0)
    virgin = SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86)
    gas = GaseousComponent(:MMA, M = 0.100, c = gas_c, k = 0.02, λ = 1e-5)
    rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
    return Material(name = :PMMA, components = (virgin, gas), reactions = (rxn,))
end

function _wood(; α = 1000.0, κ = 1e-14, law = nothing)
    virgin = SolidComponent(
        :wood,
        ρ = 500.0,
        c = 1500.0,
        k = 0.15,
        ε = 0.9,
        α = α,
        κ = κ,
        λ = 1e-5,
    )
    gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5)
    rxn = Reaction(:pyrolysis, 1 => 2, A = 1e10, E = 150e3, h = -300e3)
    return Material(
        name = :wood,
        components = (virgin, gas),
        reactions = (rxn,),
        lateral_shrinkage_law = law,
    )
end

# A 1 cm column of `n` uniform cells, adiabatic at the back.
_column(mat, n; top, bottom = AdiabaticBC(), mass = nothing) = PyrolysisProblem(
    mesh = create_uniform_mesh(0.01, n, Val(length(mat.components))),
    material = mat,
    bc_set = mass === nothing ?
             BoundaryConditionSet(thermal = (top = top, bottom = bottom)) :
             BoundaryConditionSet(thermal = (top = top, bottom = bottom), mass = mass),
    tspan = (0.0, 100.0),
)

# The configurations of `test/perf/perf_helpers.jl` as problems for the
# finite-volume core, with the incident flux of the old core's
# `radiation_intensity` keyword as a radiative flux condition on the top face.
const _incident = RadiativeFluxBC(flux = 50e3)
const _perf_configurations = (
    baseline_eulerian = (
        problem = n -> _column(_pmma(), n; top = HeatFluxBC(50e3)),
        options = (;),
        wood = false,
    ),
    eulerian_surface_absorption = (
        problem = n -> _column(_pmma(), n; top = _incident),
        options = (radiation_model = SURFACE_ABSORPTION,),
        wood = false,
    ),
    baseline_eulerian_darcy = (
        problem = n -> _column(_wood(), n; top = HeatFluxBC(50e3)),
        options = (;),
        wood = true,
    ),
    ale_no_radiation = (
        problem = n -> _column(_pmma(), n; top = HeatFluxBC(50e3)),
        options = (use_ale = true,),
        wood = false,
    ),
    ale_beer_lambert_darcy = (
        problem = n -> _column(_wood(), n; top = _incident),
        options = (use_ale = true, radiation_model = BEER_LAMBERT),
        wood = true,
    ),
    ale_beer_lambert_darcy_moderate = (
        problem = n -> _column(_wood(α = 0.5), n; top = _incident),
        options = (use_ale = true, radiation_model = BEER_LAMBERT),
        wood = true,
    ),
    ale_beer_lambert_fickian_moderate = (
        problem = n -> _column(_wood(α = 0.5, κ = 0.0), n; top = _incident),
        options = (use_ale = true, radiation_model = BEER_LAMBERT),
        wood = true,
    ),
)

# ============================================================================
# States
# ============================================================================

# The state of `cp` with cell temperatures `T(i)` and concentrations `ξ(i)`
# [kg/m³]; every other unknown as at the start.
function _cell_state(cp, T, ξ)
    u = copy(cp.u0)
    layout = cp.def.layout
    V = cp.def.setup.reference.volume
    for i = 1:n_cells(layout)
        m = Tuple(ξ(i) .* V[i])
        for j in eachindex(m)
            u[mass_index(layout, i, j)] = m[j]
        end
        u[energy_index(layout, i)] =
            energy_state_value(energy_state(layout), cp.def.setup, T(i), m, i)
    end
    return u
end

# The states of the perf configurations and of `test/jacobian/test_regressions.jl`
# (uniform hot, hot and gas-heavy, partially decomposed), the perf suite's
# mixed reacting profile, and three trial states.
function _states(cp, wood::Bool)
    n = n_cells(cp.def.layout)
    u0 = cp.u0
    solid = wood ? 500.0 : 1190.0
    uniform(T, ξ) = _cell_state(cp, i -> T, i -> ξ)
    return (;
        u0,
        mixed = _cell_state(
            cp,
            i -> 450.0 + 200.0 * i / n,
            i -> (solid * (0.65 + 0.3 * i / n), 0.02 + 0.08 * i / n),
        ),
        hot = uniform(1000.0, wood ? (400.0, 0.0) : (1100.0, 0.0)),
        gas_heavy = wood ? uniform(1000.0, (400.0, 5.0)) : uniform(1200.0, (400.0, 8.0)),
        partial = uniform(800.0, wood ? (200.0, 2.0) : (700.0, 2.0)),
        perturbed = u0 .* (1 .+ 0.2 .* sin.(1:length(u0))),
        zero_state = zeros(length(u0)),
        negative = -abs.(u0),
    )
end

# The pattern contains the dense Jacobian's nonzeros and the colored Jacobian
# equals it per entry, at a state where the residual is finite. Returns the
# dense Jacobian.
function _contained(cp, u; t = 1.0)
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    du = zeros(length(u))
    residual!(du, u, def, t, ws)
    @test all(isfinite, du)
    Jd = dense_jacobian(def, ws, u, t)
    @test all(isfinite, Jd)
    @test isempty(missing_entries(Jd, cj.pattern))
    J = jacobian!(jac_prototype(cj), cj, u, t)
    @test same_structure(J, cj.pattern)
    @test compare_entries(Matrix(J), Jd; rtol = 1e-12).passed
    return Jd
end

# ============================================================================
# The perf configurations
# ============================================================================

@testset "perf configuration $(name): $(energy) form, $(n) cells" for name in keys(
        _perf_configurations,
    ),
    energy in (:temperature, :enthalpy),
    n in (6, 8)

    config = _perf_configurations[name]
    cp = setup(config.problem(n); energy_state = energy, config.options...)
    for (state, u) in pairs(_states(cp, config.wood))
        @testset "$state" begin
            Jd = _contained(cp, u)
            layout = cp.def.layout
            if state === :hot && is_moving(layout)
                # The velocity rows read the cells' state through the strain
                # rate: the block the old core's prefix sum carried.
                v = velocity_index(layout, n_nodes(layout))
                @test Jd[v, energy_index(layout, n)] != 0
                @test Jd[v, mass_index(layout, n, 1)] != 0
            end
        end
    end
end

# ============================================================================
# Chain-rule regimes
# ============================================================================

# A reaction confined to a temperature window, its smooth gates evaluated
# within a kelvin of both edges.
function _window(n)
    virgin = SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86)
    gas = GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)
    rxn = Reaction(
        :decomp,
        1 => 2,
        A = 1.0e10,
        E = 100e3,
        h = -300e3,
        T_min = 700.0,
        T_max = 1000.0,
        gate_width = 1.0,
    )
    mat = Material(name = :window, components = (virgin, gas), reactions = (rxn,))
    return _column(mat, n; top = HeatFluxBC(50e3))
end

# Two permeable solids of different permeability and intrinsic density
# (the permeability mixing), and their pyrolysis gas.
function _two_solids(n)
    virgin = SolidComponent(
        :wood,
        ρ = 500.0,
        c = 1500.0,
        k = 0.15,
        ε = 0.9,
        α = 1000.0,
        κ = 1e-14,
        λ = 1e-5,
        ρ_intrinsic = 1500.0,
    )
    char = SolidComponent(
        :char,
        ρ = 200.0,
        c = 1100.0,
        k = 0.10,
        ε = 0.95,
        α = 5000.0,
        κ = 1e-12,
        λ = 1e-5,
        ρ_intrinsic = 1800.0,
    )
    gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5)
    rxn = Reaction(
        :pyrolysis,
        1 => (2, 3);
        A = 1e10,
        E = 150e3,
        h = -300e3,
        yields = (0.3, 0.7),
    )
    mat = Material(name = :two_solids, components = (virgin, char, gas), reactions = (rxn,))
    return _column(mat, n; top = HeatFluxBC(50e3))
end

const _regimes = (
    window_lower_edge = (
        problem = () -> _window(8),
        options = (;),
        state = cp -> _cell_state(cp, i -> 700.0 + 0.1i, i -> (1100.0, 0.5)),
    ),
    window_upper_edge = (
        problem = () -> _window(8),
        options = (;),
        state = cp -> _cell_state(cp, i -> 1000.0 - 0.1i, i -> (1100.0, 0.5)),
    ),
    gas_heat_capacity_of_temperature = (
        problem = () -> _column(
            _pmma(gas_c = LinearProperty(1100.0, 0.8)),
            8;
            top = HeatFluxBC(50e3),
        ),
        options = (;),
        state = cp -> _cell_state(cp, i -> 900.0 + 10i, i -> (700.0 - 20i, 1.0 + 0.3i)),
    ),
    darcy_viscosity_hot_gas = (
        problem = () -> _column(_wood(), 8; top = _incident),
        options = (use_ale = true, radiation_model = BEER_LAMBERT),
        state = cp -> _cell_state(cp, i -> 900.0 + 20i, i -> (400.0 - 10i, 2.0 + 0.5i)),
    ),
    single_cell = (
        problem = () -> _column(_pmma(), 1; top = HeatFluxBC(50e3)),
        options = (;),
        state = cp -> _cell_state(cp, i -> 1000.0, i -> (1100.0, 0.0)),
    ),
    single_moving_cell = (
        problem = () -> _column(_pmma(), 1; top = HeatFluxBC(50e3)),
        options = (use_ale = true,),
        state = cp -> _cell_state(cp, i -> 1000.0, i -> (1100.0, 0.0)),
    ),
    two_hundred_cells = (
        problem = () -> _column(_pmma(), 200; top = HeatFluxBC(50e3)),
        options = (;),
        state = cp -> _cell_state(cp, i -> 1000.0, i -> (1100.0, 0.5)),
    ),
    cold_below_window = (
        problem = () -> _column(
            Material(
                name = :dormant,
                components = _pmma().components,
                reactions = (
                    Reaction(
                        :decomp,
                        1 => 2,
                        A = 8.5e12,
                        E = 188e3,
                        h = -870e3,
                        T_min = 800.0,
                    ),
                ),
            ),
            8;
            top = HeatFluxBC(50e3),
        ),
        options = (;),
        state = cp -> _cell_state(cp, i -> 400.0, i -> (1100.0, 1.0)),
    ),
    permeability_mixing = (
        problem = () -> _two_solids(8),
        options = (;),
        state = cp ->
            _cell_state(cp, i -> 900.0 + 10i, i -> (300.0, 60.0 + 5i, 0.5 + 0.1i)),
    ),
)

@testset "chain-rule regime $(name): $(energy) form" for name in keys(_regimes),
    energy in (:temperature, :enthalpy)

    regime = _regimes[name]
    cp = setup(regime.problem(); energy_state = energy, regime.options...)
    u = regime.state(cp)
    Jd = _contained(cp, u)
    layout = cp.def.layout
    if name in (:window_lower_edge, :window_upper_edge)
        # The gate's slope reaches the temperature column of the mass rows.
        @test Jd[mass_index(layout, 4, 1), energy_index(layout, 4)] != 0
    end
end

# ============================================================================
# The variable cross section
# ============================================================================

@testset "the cross section's back-reaction: $(energy) form, $(n) cells" for energy in (
        :temperature,
        :enthalpy,
    ),
    n in (6, 8)

    cp = setup(
        _column(_wood(law = χ -> 1.0 - 0.2χ), n; top = HeatFluxBC(50e3));
        energy_state = energy,
        use_ale = true,
    )
    layout = cp.def.layout
    @test has_chi(layout)
    hot = _cell_state(cp, i -> 900.0 + 100.0 * (i - 1) / (n - 1), i -> (400.0, 2.0))
    for i = 1:n
        hot[history_index(layout, i, 1)] = (0.3 + 0.05i) * hot[history_index(layout, i, 2)]
    end
    # With no dry solid or no conversion the conversion fraction is zero, and
    # the lateral law and its derivatives would not be exercised.
    ΣG = sum(hot[history_index(layout, i, 1)] for i = 1:n)
    ΣM = sum(hot[history_index(layout, i, 2)] for i = 1:n)
    @test ΣM > 0
    @test ΣG / ΣM > 0.3
    for (state, u) in pairs((;
        hot,
        perturbed = hot .* (1 .+ 0.2 .* sin.(1:length(hot))),
        zero_state = zeros(length(hot)),
        negative = -abs.(hot),
    ))
        @testset "$state" begin
            Jd = _contained(cp, u)
            if state === :hot
                θ_A = lateral_rate_index(layout)
                v = velocity_index(layout, n_nodes(layout))
                @test Jd[v, θ_A] != 0
                @test Jd[θ_A, energy_index(layout, 1)] != 0
                @test Jd[area_index(layout), history_index(layout, 1, 1)] != 0
            end
        end
    end
end

# ============================================================================
# Boundary mass fluxes on the energy row
# ============================================================================

const _surface_balance =
    HeatFluxBC(50e3) + ConvectiveBC(h = 25.0, T_inf = 300.0) + RadiativeBC(T_inf = 300.0)
const _boundary_flux_cases = (
    fick_diffusive_top = (darcy = false, top = HeatFluxBC(50e3), bottom = ImpermeableBC()),
    darcy_diffusive_top = (darcy = true, top = HeatFluxBC(50e3), bottom = ImpermeableBC()),
    fick_diffusive_both = (
        darcy = false,
        top = HeatFluxBC(50e3),
        bottom = DiffusiveMassBC(),
    ),
    fick_surface_balance = (
        darcy = false,
        top = _surface_balance,
        bottom = ImpermeableBC(),
    ),
    darcy_surface_balance = (
        darcy = true,
        top = _surface_balance,
        bottom = ImpermeableBC(),
    ),
)

function _boundary_flux_column(case)
    virgin =
        case.darcy ?
        SolidComponent(
            :wood,
            ρ = 500.0,
            c = 1500.0,
            k = 0.15,
            ε = 0.9,
            κ = 1e-14,
            λ = 1e-5,
        ) : SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86, λ = 1e-5)
    gas = GaseousComponent(:gas, M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)
    rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
    mat = Material(name = :m, components = (virgin, gas), reactions = (rxn,))
    return _column(
        mat,
        8;
        top = case.top,
        mass = (top = DiffusiveMassBC(), bottom = case.bottom),
    )
end

@testset "boundary mass flux on the energy row: $(name), $(energy) form" for name in keys(
        _boundary_flux_cases,
    ),
    energy in (:temperature, :enthalpy)

    case = _boundary_flux_cases[name]
    cp = setup(_boundary_flux_column(case); energy_state = energy)
    layout = cp.def.layout
    n = n_cells(layout)
    # A temperature ramp, so both boundary faces conduct, and gas everywhere:
    # the boundary mass flux then dominates the boundary cells' energy rows.
    u = _cell_state(
        cp,
        i -> 600.0 + 400.0 * (i - 1) / (n - 1),
        i -> (case.darcy ? 400.0 : 1100.0, 8.0),
    )
    Jd = _contained(cp, u)
    top = energy_index(layout, n)
    @test Jd[top, mass_index(layout, n, 2)] != 0
    @test Jd[top, energy_index(layout, n - 1)] != 0
    if case.bottom isa DiffusiveMassBC
        @test Jd[energy_index(layout, 1), mass_index(layout, 1, 2)] != 0
    end
end

# ============================================================================
# The surface balance under surface absorption: ε(T) and an absorptivity
# ============================================================================

@testset "surface absorption with ε(T) and a component absorptivity: $(energy) form" for energy in
                                                                                         (
    :temperature,
    :enthalpy,
)
    virgin = SolidComponent(
        :virgin,
        ρ = 1100.0,
        c = (900.0, 1.5),
        k = (0.2, 2e-4),
        ε = (0.8, 5e-5),
        absorptivity = 0.7,
    )
    char = SolidComponent(:char, ρ = 300.0, c = 1100.0, k = 0.1, ε = 0.95)
    gas = GaseousComponent(:gas, M = 0.03, c = 1100.0, k = 0.03, λ = 0.0)
    mat = Material(
        name = :surface,
        components = (virgin, char, gas),
        reactions = (
            Reaction(:pyro, 1 => (2, 3); A = 5e9, E = 1.4e5, h = 8e5, yields = (0.3, 0.7)),
        ),
        depletion_limiter = DepletionLimiter(threshold = 1.0),
    )
    top = _incident + RadiativeBC(T_inf = 300.0) + ConvectiveBC(h = 10.0, T_inf = 300.0)
    cp = setup(
        _column(mat, 6; top);
        energy_state = energy,
        radiation_model = SURFACE_ABSORPTION,
    )
    layout = cp.def.layout
    n = n_cells(layout)
    u = _cell_state(cp, i -> 600.0 + 40i, i -> (800.0 - 50i, 50.0 + 20i, 0.4))
    Jd = _contained(cp, u)
    # The exposed cell's energy row against central differences in its own
    # energy and mass columns: the surface temperature's implicit derivative.
    du = zeros(length(u))
    row = energy_index(layout, n)
    for col in (energy_index(layout, n), mass_index(layout, n, 1), mass_index(layout, n, 2))
        h = 1e-6 * max(abs(u[col]), 1.0)
        up, um = copy(u), copy(u)
        up[col] += h
        um[col] -= h
        residual!(du, up, cp.def, 1.0, cp.ws)
        rp = du[row]
        residual!(du, um, cp.def, 1.0, cp.ws)
        rm = du[row]
        @test Jd[row, col] ≈ (rp - rm) / 2h rtol = 1e-5
    end
end

# ============================================================================
# The smooth face mean
# ============================================================================

@testset "the smooth face mean is finite with a finite derivative through its guard" begin
    guarded = Pyrolysis.Properties.guarded_harmonic_mean
    η = Pyrolysis.Properties.FACE_MEAN_GUARD
    d_a, d_b = 1e-3, 2e-3
    for a in (0.0, 1e-30, 1e-25, 0.5η, η, 2η, 1e-15, 0.2), b in (0.0, η, 0.45)
        g = ForwardDiff.gradient(x -> guarded(x[1], x[2], d_a, d_b), [a, b])
        @test isfinite(guarded(a, b, d_a, d_b))
        @test all(isfinite, g)
        @test guarded(a, b, d_a, d_b) >= 0
    end
    # Continuous across the scale of the guard, where the old face mean cut
    # off: no jump between just below and just above it.
    below, above =
        guarded(η * (1 - 1e-9), 0.2, d_a, d_b), guarded(η * (1 + 1e-9), 0.2, d_a, d_b)
    @test abs(above - below) <= 1e-8 * above
    # Far above the guard it is the harmonic mean.
    @test guarded(0.15, 0.45, d_a, d_b) ≈
          Pyrolysis.Properties.harmonic_mean_weighted(0.15, 0.45, d_a, d_b) rtol = 1e-14
end
