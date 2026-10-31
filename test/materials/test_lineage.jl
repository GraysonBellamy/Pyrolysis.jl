# The reference of a reaction of order other than one: the graph rules of
# `lineage_graph` and the per-cell rules of `lineage_references`.

using Test
using Pyrolysis
using Pyrolysis.Materials: lineage_graph, needs_lineage
using Pyrolysis.Integration: lineage_references, setup
using Pyrolysis.Kernels: ReactionKernel
using Pyrolysis.Closures: EquidistributionCorrection, LayerCorrection
using Pyrolysis.Meshes: box_mesh
using Pyrolysis.Assembly: NoSources, VolumetricSources

_sol(name) = SolidComponent(name, ρ = 500.0, c = 1000.0, k = 0.2, φ = 0.3, λ = 1e-6)
_gas(name = :gas) = GaseousComponent(name, M = 0.03, c = 1000.0, k = 0.02)

_material(components, reactions) =
    Material(name = :lineage, components = components, reactions = reactions)

# One channel of `a`: at most one condensed product at its yield, and the
# remainder to the gas.
function _step(name, a, products; n = 2.0, basis = NaN)
    isempty(products) &&
        return Reaction(name, a => :gas; A = 1.0, E = 0.0, h = 0.0, n, basis)
    p, y = only(products)
    return Reaction(
        name,
        a => (p, :gas);
        A = 1.0,
        E = 0.0,
        h = 0.0,
        yields = (y, 1 - y),
        n,
        basis,
    )
end

# The references of `material` for the initial masses `m0` (`NC × n`).
function _references(material, m0; sources = NoSources())
    condensed = [!is_gas(c) for c in material.components]
    condensed0 = vec(sum(m0[condensed, :]; dims = 1))
    return lineage_references(lineage_graph(material), m0, condensed0, sources)
end

# The reach of reaction `k` onto its reactant from the component `v`.
_reach(material, k, v) = lineage_graph(material).reach[k][v]

@testset "needs_lineage" begin
    first_order = _material((_sol(:a), _gas()), (_step(:r, :a, (); n = 1.0),))
    second_order = _material((_sol(:a), _gas()), (_step(:r, :a, ()),))
    @test !needs_lineage(first_order)
    @test needs_lineage(second_order)
    @test !needs_lineage(_material((_sol(:a), _gas()), ()))
end

@testset "Douglas fir chain" begin
    # Rounded inputs: 94.3 mg wood and 5.7 mg water in 100 mg, yields 0.95,
    # 0.80, 0.37, 0.79, 0.74 along the chain.
    names = (:wood, :int1, :int2, :int3, :int4, :char, :water, :gas)
    components = (
        _sol(:wood),
        _sol(:int1),
        _sol(:int2),
        _sol(:int3),
        _sol(:int4),
        _sol(:char),
        LiquidComponent(:water, ρ = 1000.0, c = 4180.0, k = 0.6, λ = 1e-6),
        _gas(),
    )
    reactions = (
        _step(:r1, :wood, (:int1 => 0.95,)),
        _step(:r2, :int1, (:int2 => 0.80,)),
        _step(:r3, :int2, (:int3 => 0.37,)),
        _step(:r4, :int3, (:int4 => 0.79,)),
        _step(:r5, :int4, (:char => 0.74,)),
        _step(:r6, :char, ()),
        _step(:r7, :water, ()),
    )
    material = _material(components, reactions)
    m0 = zeros(8, 2)
    m0[1, :] .= 94.3
    m0[7, :] .= 5.7
    r = _references(material, m0)
    expected = (94.3, 89.585, 71.668, 26.517, 20.949, 15.502, 5.7)
    @test all(k -> all(≈(expected[k]; atol = 5e-4), r[k]), 1:7)
    @test r[6][1] ≈ 94.3 * 0.95 * 0.80 * 0.37 * 0.79 * 0.74 rtol = 1e-13
end

@testset "Sandia parallel and series" begin
    # Parallel: two lineages side by side, the second 0.91 of the mass.
    parallel =
        _material((_sol(:a), _sol(:b), _gas()), (_step(:ra, :a, ()), _step(:rb, :b, ())))
    r = _references(parallel, reshape([0.09, 0.91, 0.0], 3, 1))
    @test only(r[1]) ≈ 0.09
    @test only(r[2]) ≈ 0.91
    # Series: the intermediate's reference is 0.9 of the mass.
    series = _material(
        (_sol(:a), _sol(:b), _gas()),
        (_step(:ra, :a, (:b => 0.9,)), _step(:rb, :b, ())),
    )
    r = _references(series, reshape([1.0, 0.0, 0.0], 3, 1))
    @test only(r[1]) ≈ 1.0
    @test only(r[2]) ≈ 0.9
end

@testset "reach over competing channels" begin
    three = (_sol(:a), _sol(:b), _gas())
    four = (_sol(:a), _sol(:b), _sol(:c), _sol(:d), _gas())
    target(s) = _step(:z, s, ())
    # Two equal channels count once.
    m = _material(
        three,
        (_step(:p, :a, (:b => 0.5,)), _step(:q, :a, (:b => 0.5,)), target(:b)),
    )
    @test _reach(m, 3, 1) == 0.5
    # Two routes to the same node do not add.
    m = _material(
        four,
        (
            _step(:ab, :a, (:b => 1.0,)),
            _step(:ac, :a, (:c => 1.0,)),
            _step(:bd, :b, (:d => 1.0,)),
            _step(:cd, :c, (:d => 1.0,)),
            target(:d),
        ),
    )
    @test _reach(m, 5, 1) == 1.0
    # A channel producing nothing toward the target is ignored.
    m = _material(
        (_sol(:a), _sol(:b), _sol(:c), _gas()),
        (_step(:ab, :a, (:b => 0.5,)), _step(:ac, :a, (:c => 1.0,)), target(:b)),
    )
    @test _reach(m, 3, 1) == 0.5
    # Reaching channels that disagree depend on the rates.
    @test_throws "competing channels" _material(
        three,
        (_step(:p, :a, (:b => 0.5,)), _step(:q, :a, (:b => 0.2,)), target(:b)),
    )
    # One channel splitting and rejoining.
    m = _material(
        four,
        (
            Reaction(
                :split,
                :a => (:b, :c);
                A = 1.0,
                E = 0.0,
                h = 0.0,
                yields = (0.4, 0.6),
            ),
            _step(:bd, :b, (:d => 1.0,)),
            _step(:cd, :c, (:d => 1.0,)),
            target(:d),
        ),
    )
    @test _reach(m, 4, 1) ≈ 1.0
    # The self-return renormalises.
    self = Reaction(:self, :a => (:a, :b); A = 1.0, E = 0.0, h = 0.0, yields = (0.4, 0.6))
    m = _material((_sol(:a), _sol(:b), _gas()), (self, target(:b)))
    @test _reach(m, 2, 1) == 1.0
end

@testset "a reaction returning its whole reactant has no edge" begin
    m = _material(
        (_sol(:a), _gas()),
        (Reaction(:heat, 1 => 1; A = 1.0, E = 0.0, h = 1e5, n = 0.0),),
    )
    r = _references(m, reshape([2.0, 0.0], 2, 1))
    @test only(r[1]) == 2.0
end

@testset "mixed origin, per cell" begin
    m = _material(
        (_sol(:a), _sol(:b), _sol(:c), _gas()),
        (
            _step(:ab, :a, (:b => 1.0,); n = 1.0),
            _step(:bc, :b, (:c => 0.5,); n = 1.0),
            _step(:cg, :c, ()),
        ),
    )
    @test_throws "give `basis`" _references(m, reshape([1.0, 0.0, 0.2, 0.0], 4, 1))
    @test only(_references(m, reshape([0.0, 0.0, 0.2, 0.0], 4, 1))[3]) ≈ 0.2
    @test only(_references(m, reshape([1.0, 0.0, 0.0, 0.0], 4, 1))[3]) ≈ 0.5
    # The refusal is per cell: the cells are told apart.
    @test_throws "cell 2" _references(m, [0.0 1.0; 0.0 0.0; 0.2 0.2; 0.0 0.0])
    # `basis` is the remedy.
    with_basis = _material(
        (_sol(:a), _sol(:b), _sol(:c), _gas()),
        (
            _step(:ab, :a, (:b => 1.0,); n = 1.0),
            _step(:bc, :b, (:c => 0.5,); n = 1.0),
            _step(:cg, :c, (); basis = 0.5),
        ),
    )
    @test only(_references(with_basis, reshape([1.0, 0.0, 0.2, 0.0], 4, 1))[3]) ≈ 0.6
end

@testset "basis is masked to the cells its lineage reaches" begin
    m = _material((_sol(:inert), _sol(:a), _gas()), (_step(:r, :a, (); basis = 0.25),))
    # Cell 1 is pure inert, cell 2 pure reactant, cell 3 a blend.
    m0 = [4.0 0.0 1.0; 0.0 2.0 1.0; 0.0 0.0 0.0]
    @test _references(m, m0)[1] ≈ [0.0, 0.5, 0.5]
end

@testset "external feed" begin
    # Gas forms the reactant: refused with or without a basis, and even when
    # the reactant is absent, since the gas can still form it.
    for basis in (NaN, 1.0)
        @test_throws "gas-consuming" _material(
            (_sol(:a), _gas()),
            (_step(:r, :a, (); basis), Reaction(:cond, 2 => 1; A = 1.0, E = 0.0, h = 0.0)),
        )
    end
    # A node downstream of a gas-fed node is fed too.
    @test_throws "gas-consuming" _material(
        (_sol(:a), _sol(:b), _gas()),
        (
            Reaction(:cond, 3 => 1; A = 1.0, E = 0.0, h = 0.0),
            _step(:ab, :a, (:b => 1.0,); n = 1.0),
            _step(:r, :b, ()),
        ),
    )
    # A zero-yield product of the gas reaction feeds nothing.
    m = _material(
        (_sol(:a), _sol(:b), _gas()),
        (
            Reaction(:cond, 3 => (1, 2); A = 1.0, E = 0.0, h = 0.0, yields = (0.0, 1.0)),
            _step(:r, :a, ()),
        ),
    )
    @test only(_references(m, reshape([1.0, 0.0, 0.0], 3, 1))[2]) == 1.0
    # A gas reactant of order other than one has no condensed inventory.
    @test_throws "gas reactant" _material(
        (_sol(:a), _gas()),
        (Reaction(:r, 2 => 1; A = 1.0, E = 0.0, h = 0.0, n = 2.0),),
    )
    # External mass sources feed every condensed node, whatever its basis.
    for basis in (NaN, 1.0)
        m = _material((_sol(:a), _gas()), (_step(:r, :a, (); basis),))
        sources = VolumetricSources((x, t) -> (0.0, 0.0), nothing)
        @test_throws "sources add mass" _references(m, reshape([1.0, 0.0], 2, 1); sources)
    end
    m = _material((_sol(:a), _gas()), (_step(:r, :a, ()),))
    @test only(
        _references(
            m,
            reshape([1.0, 0.0], 2, 1);
            sources = VolumetricSources(nothing, (x, t)->0.0),
        )[1],
    ) == 1.0
end

@testset "an empty cell has exact zeros" begin
    m = _material(
        (_sol(:a), _sol(:b), _gas()),
        (_step(:ab, :a, (:b => 0.5,)), _step(:bg, :b, ())),
    )
    r = _references(m, [1.0 0.0; 0.0 0.0; 0.0 3.0])
    @test r[1] == [1.0, 0.0]
    @test r[2] == [0.5, 0.0]
end

@testset "a zero-yield product is not an edge" begin
    # B is listed with yield zero in A's reaction, so A is no ancestor of B.
    m = _material(
        (_sol(:a), _sol(:b), _gas()),
        (
            Reaction(:r, 1 => (2, 3); A = 1.0, E = 0.0, h = 0.0, yields = (0.0, 1.0)),
            _step(:bg, :b, ()),
        ),
    )
    r = _references(m, reshape([1.0, 0.0, 0.0], 3, 1))
    @test r[2] == [0.0]
end

@testset "representability" begin
    m = _material((_sol(:a), _gas()), (_step(:r, :a, (); basis = 1e-7),))
    @test_throws "below 1e-6" _references(m, reshape([1.0, 0.0], 2, 1))
    m = _material((_sol(:a), _gas()), (_step(:r, :a, (); basis = 1e-6),))
    @test only(_references(m, reshape([1.0, 0.0], 2, 1))[1]) ≈ 1e-6
    # A vanishing fraction of a mixture is not a normalisation.
    m = _material((_sol(:a), _sol(:b), _gas()), (_step(:r, :b, ()),))
    @test_throws "below 1e-12" _references(m, reshape([1.0, 1e-13, 0.0], 3, 1))
    # A long chain of small yields underflows the same test.
    components = (_sol(:a), _sol(:b), _sol(:c), _sol(:d), _gas())
    m = _material(
        components,
        (
            _step(:r1, :a, (:b => 1e-5,); n = 1.0),
            _step(:r2, :b, (:c => 1e-5,); n = 1.0),
            _step(:r3, :c, (:d => 1e-5,); n = 1.0),
            _step(:r4, :d, ()),
        ),
    )
    @test_throws "below 1e-12" _references(m, reshape([1.0, 0, 0, 0, 0], 5, 1))
end

@testset "an inactive reaction" begin
    m = _material(
        (_sol(:a), _sol(:b), _sol(:c), _gas()),
        (_step(:live, :a, ()), _step(:dead, :b, (:c => 1.0,); n = 0.5)),
    )
    r = _references(m, reshape([1.0, 0.0, 0.0, 0.0], 4, 1))
    @test r[2] == [0.0]
end

@testset "cycles" begin
    # A cycle among the ancestors is refused, `basis` the remedy.
    cycle(basis) = _material(
        (_sol(:a), _sol(:b), _sol(:c), _gas()),
        (
            _step(:ab, :a, (:b => 1.0,); n = 1.0),
            _step(:ba, :b, (:a => 0.5,); n = 1.0),
            _step(:bc, :b, (:c => 0.5,); n = 1.0),
            _step(:cg, :c, (); basis),
        ),
    )
    @test_throws "cycle" cycle(NaN)
    @test_throws "a → b → a" cycle(NaN)
    m = cycle(1.0)
    @test only(_references(m, reshape([1.0, 0.0, 0.0, 0.0], 4, 1))[4]) == 1.0
    # An all-first-order cyclic network is never validated.
    m = _material(
        (_sol(:a), _sol(:b), _gas()),
        (_step(:ab, :a, (:b => 1.0,); n = 1.0), _step(:ba, :b, (:a => 1.0,); n = 1.0)),
    )
    @test !needs_lineage(m)
end

@testset "references in setup" begin
    bc = BoundaryConditionSet(
        thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
        mass = (
            top = AmbientGasBC(GasReservoir(P = 0.0, T = 600.0)),
            bottom = ImpermeableBC(),
        ),
    )
    problem(material, ξ) = PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, 4),
        material = material,
        bc_set = bc,
        T_initial = z -> 600.0,
        ξ_initial = ξ,
        tspan = (0.0, 1.0),
    )
    ab = (_sol(:a), _sol(:b), _gas())
    layered(z) = z < 0.005
    ξ = [z -> layered(z) ? 400.0 : 100.0, z -> layered(z) ? 100.0 : 400.0, z -> 0.0]
    kernel(cp) = only(k for k in cp.def.kernels if k isa ReactionKernel)
    # Two layers: each reaction's own lineage, cell by cell.
    cp = setup(problem(_material(ab, (_step(:ra, :a, ()), _step(:rb, :b, ()))), ξ))
    @test kernel(cp).reference ≈ cp.m0[1:2, :]
    # The reference is a constant of the cell: the state has no row for it.
    first_order = setup(problem(_material(ab, (_step(:ra, :a, (); n = 1.0),)), ξ))
    second_order = setup(problem(_material(ab, (_step(:ra, :a, ()),)), ξ))
    @test length(second_order.u0) == length(first_order.u0)
    @test kernel(second_order) isa ReactionKernel{true}
    @test kernel(first_order) isa ReactionKernel{false}
    @test all(iszero, kernel(first_order).reference)
    # Multi-D meshes refuse.
    mesh = box_mesh(0:0.5:1.0, 0:0.5:1.0)
    multi = PyrolysisProblem(
        mesh = mesh,
        material = _material(
            (_sol(:a),),
            (Reaction(:r, :a => :a; A = 1.0, E = 0.0, h = 0.0, n = 2.0),),
        ),
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
        ),
        T_initial = z -> 600.0,
        ξ_initial = [z -> 500.0],
        tspan = (0.0, 1.0),
    )
    @test_throws "a reaction order other than one" setup(multi)
end

@testset "an order other than one keeps the cells' material" begin
    gas = GaseousComponent(:gas, M = 0.03, c = 1000.0, k = 0.02)
    solid = SolidComponent(:a; φ = 0.3, ρ = 500.0, c = 1000.0, k = 0.2, λ = 1e-6)
    problem(n) = PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, 4),
        material = Material(
            name = :order,
            components = (solid, gas),
            reactions = (Reaction(:burn, 1 => 2; A = 1.0, E = 0.0, h = 0.0, n),),
        ),
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (
                top = AmbientGasBC(GasReservoir(P = 0.0, T = 600.0)),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> 600.0,
        ξ_initial = [z -> 500.0, z -> 0.0],
        tspan = (0.0, 1.0),
    )
    core = Pyrolysis.Integration.CoreProblem
    for n in (0.0, 0.5, 2.0)
        @test setup(problem(n); use_ale = true) isa core
        @test setup(problem(n)) isa core
        @test_throws "solve without deactivation" setup(problem(n); deactivation = true)
        @test_throws "solve without deactivation" setup(
            problem(n);
            use_ale = true,
            deactivation = true,
        )
        @test_throws "give no correction" setup(
            problem(n);
            use_ale = true,
            correction = EquidistributionCorrection(layers = :all),
        )
    end
    # A first-order deck is untouched by the rule.
    cp = setup(problem(1.0); use_ale = true, deactivation = true)
    @test cp.def.closure.correction isa LayerCorrection
end
