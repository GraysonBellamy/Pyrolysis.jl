using Test
using Pyrolysis
using Pyrolysis.Integration: _condensed_volume_ratio
using Pyrolysis.Assembly: VolumetricSources

const _ADIABATIC =
    BoundaryConditionSet(thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()))

# A charring solid on a line of six cells, `A → 0.3 C + 0.7 G`, whose char
# holds `0.3 ρ_A/ρ_C` of its reactant's volume: 0.6 by default, and all of it at
# `ρ_C = 0.3 ρ_A`.
function _charring(; ρ_A = 500.0, ρ_C = 250.0, ξ_A = z -> ρ_A)
    material = Material(
        name = :charring,
        components = (
            SolidComponent(:A; ρ = ρ_A, c = 1500.0, k = 0.2, φ = 0.6, λ = 0.0),
            SolidComponent(:C; ρ = ρ_C, c = 1000.0, k = 0.1, φ = 0.6, λ = 0.0),
            GaseousComponent(:G; M = 0.03, c = 1100.0, k = 0.03),
        ),
        reactions = (
            Reaction(:pyro, 1 => (2, 3); A = 1e8, E = 1.2e5, h = 0.0, yields = (0.3, 0.7)),
        ),
    )
    return PyrolysisProblem(
        mesh = uniform_mesh(0.01, 6),
        material = material,
        bc_set = _ADIABATIC,
        T_initial = z -> 600.0,
        ξ_initial = [ξ_A, z -> 0.0, z -> 0.0],
        tspan = (0.0, 5.0),
    )
end

# Wood of oven-dry density 500 kg/m³ at its initial concentration `ξ_wood`,
# holding bound water at dry-basis moisture 0.1, which evaporates; `ρ_water`
# the water's bulk density.
function _moist(ξ_wood; ρ_water = 1000.0)
    material = Material(
        name = :moist,
        components = (
            SolidComponent(:wood; ρ = 500.0, c = 1500.0, k = 0.12, φ = 0.6, λ = 0.0),
            LiquidComponent(:water; ρ = ρ_water, c = 4180.0, k = 0.6),
            GaseousComponent(:vapour; M = 0.018, c = 2000.0, k = 0.025),
        ),
        reactions = (Reaction(:drying, 2 => 3; A = 1e5, E = 5e4, h = 2.3e6),),
    )
    return PyrolysisProblem(
        mesh = uniform_mesh(0.01, 6),
        material = material,
        bc_set = _ADIABATIC,
        T_initial = z -> 300.0,
        ξ_initial = [z -> ξ_wood, z -> 0.1 * ξ_wood, z -> 0.0],
        tspan = (0.0, 5.0),
    )
end

# An inert solid whose bulk density is the law `ρ`.
function _inert(ρ; mesh = uniform_mesh(0.01, 4), bc_set = _ADIABATIC)
    material = Material(
        name = :inert,
        components = (SolidComponent(:melt; ρ = ρ, c = 1500.0, k = 0.2),),
        reactions = (),
    )
    return PyrolysisProblem(
        mesh = mesh,
        material = material,
        bc_set = bc_set,
        T_initial = x -> 300.0,
        ξ_initial = [x -> 1000.0],
        tspan = (0.0, 1.0),
    )
end

_refusal(f) =
    try
        f()
        nothing
    catch err
        err
    end

const _MELT = T -> 1000.0 - 0.1 * (T - 300.0)

@testset "a moving mesh refuses a condensed density law" begin
    err = _refusal(() -> setup(_inert(_MELT); use_ale = true))
    @test err isa ArgumentError
    @test occursin("component :melt of material inert", err.msg)
    # A liquid's too, naming the liquid.
    err = _refusal(() -> setup(_moist(476.0; ρ_water = T -> 1000.0); use_ale = true))
    @test err isa ArgumentError
    @test occursin("component :water", err.msg)
    # Under the potential closure of a box as well.
    adiabatic = AdiabaticBC()
    box = _inert(
        _MELT;
        mesh = box_mesh([0.0, 0.01], [0.0, 0.005, 0.01]),
        bc_set = BoundaryConditionSet(
            thermal = (
                top = adiabatic,
                bottom = adiabatic,
                left = adiabatic,
                right = adiabatic,
            ),
        ),
    )
    @test _refusal(() -> setup(box; use_ale = true)) isa ArgumentError
    # A fixed mesh takes the law, and a moving one the constant.
    @test (@test_logs setup(_inert(_MELT))) isa PreparedProblem
    @test (@test_logs setup(_inert(1000.0); use_ale = true)) isa PreparedProblem
end

@testset "a moving mesh warns when the condensed components do not fill a cell" begin
    # Wood at its dry density with the water added on top over-fills every
    # cell by the water's volume, 0.1 × 500/1000.
    @test_logs (:warn, r"do not fill 6 of 6 cells .* is 1\.05 there") setup(
        _moist(500.0);
        use_ale = true,
    )
    # Wood and water sharing the cell fill it.
    @test_logs setup(_moist(500.0 / (1 + 0.1 * 500.0 / 1000.0)); use_ale = true)
    # Half the cells under-filled.
    @test_logs (
        :warn,
        r"do not fill 3 of 6 cells at the start, the first cell 4: .* is 0\.9 there",
    ) setup(_charring(ξ_A = z -> z < 0.005 ? 500.0 : 450.0); use_ale = true)
    # Within the tolerance, and outside it.
    @test_logs setup(_charring(ξ_A = z -> 500.0 * (1 + 5e-4)); use_ale = true)
    @test_logs (:warn,) setup(_charring(ξ_A = z -> 500.0 * (1 + 2e-3)); use_ale = true)
    # A cell with no condensed mass is a gap, not a cell left part empty: a
    # solid that densifies to half its volume, with a gap in cells 3 and 4.
    densifying = Material(
        name = :densifying,
        components = (
            SolidComponent(:A; ρ = 500.0, c = 1500.0, k = 0.2),
            SolidComponent(:B; ρ = 1000.0, c = 1000.0, k = 0.3),
        ),
        reactions = (Reaction(:densify, 1 => 2; A = 1e8, E = 1.2e5, h = 0.0),),
    )
    gap = PyrolysisProblem(
        mesh = uniform_mesh(0.01, 6),
        material = densifying,
        bc_set = _ADIABATIC,
        T_initial = z -> 600.0,
        ξ_initial = [z -> 0.004 < z < 0.006 ? 0.0 : 500.0, z -> 0.0],
        tspan = (0.0, 5.0),
    )
    @test_logs setup(gap; use_ale = true)
    # An external mass source moves the occupancy whatever the start.
    sources = VolumetricSources((z, t) -> (1.0, 0.0, 0.0), (z, t) -> 0.0)
    @test_logs setup(_charring(ξ_A = z -> 450.0); use_ale = true, sources = sources)
    # When no reaction changes the condensed volume the mesh never moves from
    # the initial state, so a cell may start part empty.
    @test_logs setup(_charring(ρ_C = 150.0, ξ_A = z -> 400.0); use_ale = true)
end

@testset "a fixed mesh says once that the skeleton keeps its volume" begin
    notice =
        r"reaction :pyro of material charring leaves condensed products holding 0\.6 of its reactant's volume, but the mesh is fixed"
    problem = _charring()
    @test_logs (:info, notice) (setup(problem); setup(problem))
    @test_logs setup(problem; use_ale = true)
    # Liquid that evaporates leaves no condensed volume.
    @test_logs (
        :info,
        r"reaction :drying of material moist leaves condensed products holding 0\.0",
    ) setup(_moist(500.0))
    # A char that keeps its reactant's volume strains nothing.
    @test_logs setup(_charring(ρ_C = 150.0))
    # A solve prints it once.
    logs, sol = Test.collect_test_logs(() -> solve(problem))
    @test sol.retcode === :Success
    @test count(l -> occursin("the mesh is fixed", l.message), logs) == 1
end

@testset "the condensed volume a reaction leaves" begin
    reaction = first(_charring().material.reactions)
    @test _condensed_volume_ratio(_charring().material, reaction, 300.0) ≈ 0.6
    @test _condensed_volume_ratio(_charring(ρ_C = 150.0).material, reaction, 300.0) ≈ 1.0
    # A gas that deposits a solid adds condensed volume from none.
    deposit = Material(
        name = :deposit,
        components = (
            SolidComponent(:soot; ρ = 1800.0, c = 700.0, k = 1.0, φ = 0.5, λ = 0.0),
            GaseousComponent(:G; M = 0.03, c = 1100.0, k = 0.03),
        ),
        reactions = (Reaction(:deposit, 2 => 1; A = 1.0, E = 0.0, h = 0.0),),
    )
    @test _condensed_volume_ratio(deposit, first(deposit.reactions), 300.0) == Inf
end
