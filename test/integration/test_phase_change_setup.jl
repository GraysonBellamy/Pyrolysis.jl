"""What `setup` says and refuses for a material with a phase change: the fill
warning and fixed-mesh notice, the two warnings of an open reservoir that does
not name the vapour and of a start off equilibrium, and the two refusals a
reaction would also meet."""

using Test
using Logging
using Pyrolysis

const _ADIABATIC =
    BoundaryConditionSet(thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()))

const _MOISTURE =
    PhaseChange(:moisture, :water => :vapour; activity = CAST_SPRUCE, sorbents = (:wood,))

# The `_moist` fixture of test_condensed_volume_checks.jl with the phase change
# in place of its drying reaction: wood of oven-dry density 500 kg/m³ at its
# initial concentration `ξ_wood`, holding water at dry-basis moisture 0.1, the
# pore vapour in equilibrium with it unless `vapour` is given.
function _moist(ξ_wood; vapour = nothing, T = 300.0)
    material = Material(
        name = :moist,
        components = (
            SolidComponent(:wood; ρ = 500.0, c = 1500.0, k = 0.12, φ = 0.6, λ = 0.0),
            LiquidComponent(:water; ρ = 1000.0, c = 4180.0, k = 0.6),
            GaseousComponent(:vapour; M = 0.018, c = 2000.0, k = 0.025),
        ),
        phase_changes = (_MOISTURE,),
    )
    ξ_eq = equilibrium_vapour_concentration(material, (ξ_wood, 0.1ξ_wood, 0.0), T)
    ξ_v = something(vapour, ξ_eq)
    return PyrolysisProblem(
        mesh = uniform_mesh(0.01, 6),
        material = material,
        bc_set = _ADIABATIC,
        T_initial = z -> T,
        ξ_initial = [z -> ξ_wood, z -> 0.1 * ξ_wood, z -> ξ_v],
        tspan = (0.0, 5.0),
    )
end

# Wood that fills its cells with its water, in air, its top face open to the
# reservoir mole fractions `fractions` (a NamedTuple or a function of time),
# the pore gas in equilibrium.
function _open(fractions; T = 300.0)
    material = Material(
        name = :open,
        components = (
            SolidComponent(:wood; ρ = 500.0, c = 1500.0, k = 0.12, φ = 0.6, λ = 5e-7),
            LiquidComponent(:water; ρ = 1000.0, c = 4180.0, k = 0.6),
            GaseousComponent(:vapour; M = 0.018, c = 2000.0, k = 0.025),
            GaseousComponent(:air; M = 0.02897, c = 1006.0, k = 0.026),
        ),
        phase_changes = (_MOISTURE,),
    )
    ξ_wood = 500.0 / 1.05
    ξc = (ξ_wood, 0.1ξ_wood, 0.0, 0.0)
    ξ_v = equilibrium_vapour_concentration(material, ξc, T)
    p_v = equilibrium_vapour_pressure(material, ξc, T)
    ξ_air = background_gas_concentration(material, :air, ξc, T; P = 101325.0 - p_v)
    return PyrolysisProblem(
        mesh = uniform_mesh(0.01, 6),
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (
                top = AmbientGasBC(
                    GasReservoir(P = 101325.0, T = T, mole_fractions = fractions),
                ),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> T,
        ξ_initial = [z -> ξ_wood, z -> 0.1ξ_wood, z -> ξ_v, z -> ξ_air],
        tspan = (0.0, 5.0),
    )
end

_refusal(f) =
    try
        f()
        nothing
    catch err
        err
    end

@testset "a phase change strains a moving mesh and is noticed on a fixed one" begin
    # Wood at its dry density with the water added on top over-fills every
    # cell by the water's volume, 0.1 × 500/1000.
    @test_logs (:warn, r"do not fill 6 of 6 cells .* is 1\.05 there") setup(
        _moist(500.0);
        use_ale = true,
    )
    # Wood and water sharing the cell fill it.
    @test_logs setup(_moist(500.0 / (1 + 0.1 * 500.0 / 1000.0)); use_ale = true)
    # On a fixed mesh, once a session, naming the phase change.
    notice =
        r"phase change :moisture of material moist turns its liquid into vapour, which holds none of the liquid's volume, but the mesh is fixed"
    problem = _moist(500.0)
    @test_logs (:info, notice) (setup(problem); setup(problem))
    # a moving mesh warns of the fill instead
    @test_logs (:warn, r"do not fill 6 of 6 cells") setup(problem; use_ale = true)
end

@testset "an open reservoir that does not name the vapour" begin
    dry_sink =
        r"the reservoir of the :top mass condition does not name :vapour, the vapour of phase change :moisture, so the face is a perfectly dry sink"
    # the moving mesh keeps the fixed-mesh notice out of the logs
    @test_logs (:warn, dry_sink) setup(_open((air = 1.0,)); use_ale = true)
    @test_logs (:warn, dry_sink) setup(_open(t -> (air = 1.0,)); use_ale = true)
    # named, at any fraction including zero
    @test_logs setup(_open((air = 1.0, vapour = 0.0)); use_ale = true)
    @test_logs setup(_open((air = 0.98, vapour = 0.02)); use_ale = true)
    @test_logs setup(_open(t -> (air = 1.0, vapour = 0.0)); use_ale = true)
    # a vacuum is a perfectly dry sink as well
    vacuum = _open((air = 1.0,))
    vacuum = PyrolysisProblem(
        mesh = vacuum.mesh,
        material = vacuum.material,
        bc_set = BoundaryConditionSet(
            thermal = vacuum.bc_set.thermal,
            mass = (
                top = AmbientGasBC(GasReservoir(P = 0.0, T = 300.0)),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = vacuum.T_initial,
        ξ_initial = vacuum.ξ_initial,
        tspan = vacuum.tspan,
    )
    @test_logs (:warn, dry_sink) setup(vacuum; use_ale = true)
    # no phase change, no warning
    dry = _open((air = 1.0,))
    dry_material = Material(name = :dry, components = dry.material.components)
    @test_logs setup(
        PyrolysisProblem(
            mesh = dry.mesh,
            material = dry_material,
            bc_set = dry.bc_set,
            T_initial = dry.T_initial,
            ξ_initial = dry.ξ_initial,
            tspan = dry.tspan,
        );
        use_ale = true,
    )
end

@testset "a start off equilibrium" begin
    filled = 500.0 / (1 + 0.1 * 500.0 / 1000.0)
    ξ_eq = equilibrium_vapour_concentration(
        _moist(filled).material,
        (filled, 0.1filled, 0.0),
        300.0,
    )
    # dry pores: the warning gives the farthest cell, the two values and the helper
    start =
        r"the pore vapour :vapour starts off equilibrium with the liquid of phase change :moisture, farthest in cell 1: 0\.0 kg/m³ against .* kg/m³ in equilibrium.*equilibrium_vapour_concentration gives the start in equilibrium"
    @test_logs (:warn, start) setup(_moist(filled; vapour = 0.0); use_ale = true)
    # once per setup, not per cell
    logs, _ =
        Test.collect_test_logs(() -> setup(_moist(filled; vapour = 0.0); use_ale = true))
    @test count(l -> occursin("starts off equilibrium", l.message), logs) == 1
    # in equilibrium, and within one per cent of it
    @test_logs setup(_moist(filled); use_ale = true)
    @test_logs setup(_moist(filled; vapour = 1.009ξ_eq); use_ale = true)
    @test_logs (:warn, r"starts off equilibrium") setup(
        _moist(filled; vapour = 1.011ξ_eq);
        use_ale = true,
    )
    # the farthest cell is named: vapour above equilibrium in cell 4 only
    profile = _moist(filled)
    farthest = PyrolysisProblem(
        mesh = profile.mesh,
        material = profile.material,
        bc_set = profile.bc_set,
        T_initial = profile.T_initial,
        ξ_initial = [
            profile.ξ_initial[1:2]...,
            z -> 0.005 < z < 0.0067 ? 1.5ξ_eq : 1.02ξ_eq,
        ],
        tspan = profile.tspan,
    )
    @test_logs (:warn, r"farthest in cell 4") setup(farthest; use_ale = true)
    # a cell with no liquid is not tested
    dry_cells = PyrolysisProblem(
        mesh = profile.mesh,
        material = profile.material,
        bc_set = profile.bc_set,
        T_initial = profile.T_initial,
        ξ_initial = [
            z -> filled,
            z -> z < 0.005 ? 0.1filled : 0.0,
            z -> z < 0.005 ? ξ_eq : 0.0,
        ],
        tspan = profile.tspan,
    )
    @test_logs (:warn, r"do not fill 3 of 6 cells") setup(dry_cells; use_ale = true)
end

@testset "the refusals a reaction meets" begin
    # a prescribed mesh velocity moves the mesh over a solid at rest
    err = _refusal(
        () -> setup(
            _moist(500.0 / 1.05);
            use_ale = true,
            correction = PrescribedCorrection((x, t) -> -1e-3 * x),
        ),
    )
    @test err isa ArgumentError
    @test occursin(
        "a prescribed mesh velocity moves the mesh over a solid at rest, so the material must not change phase (it has 1 phase changes)",
        err.msg,
    )
    # a series λ channel with a liquid of zero coefficient: condensation can
    # bring the liquid into a cell, or take it out, so its support changes
    function series(phase_changes; water = z -> 20.0)
        material = Material(
            name = :series,
            components = (
                SolidComponent(:wood; ρ = 500.0, c = 1500.0, k = 0.12, φ = 0.6, λ = 5e-7),
                LiquidComponent(:water; ρ = 1000.0, c = 4180.0, k = 0.6, λ = 0.0),
                GaseousComponent(:vapour; M = 0.018, c = 2000.0, k = 0.025),
            ),
            phase_changes = phase_changes,
            mixing = MaterialMixing(λ = MixingSpec(SERIES)),
        )
        return PyrolysisProblem(
            mesh = uniform_mesh(0.01, 4),
            material = material,
            bc_set = _ADIABATIC,
            T_initial = z -> 300.0,
            ξ_initial = [z -> 450.0, water, z -> 0.01],
            tspan = (0.0, 1.0),
        )
    end
    support =
        r"active series λ channel can change support of zero-coefficient component\(s\) \[:water\]"
    for water in (z -> 20.0, z -> z < 0.005 ? 20.0 : 0.0)
        err = _refusal(() -> setup(series((_MOISTURE,); water)))
        @test err isa ArgumentError
        @test occursin(support, err.msg)
    end
    # without the phase change the water's support never changes
    @test setup(series(()); use_ale = false) isa PreparedProblem
    @test setup(series((); water = z -> z < 0.005 ? 20.0 : 0.0)) isa PreparedProblem
end
