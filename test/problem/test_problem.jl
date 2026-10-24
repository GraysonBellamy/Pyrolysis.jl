using Test
using StaticArrays
using Pyrolysis
using Pyrolysis.Problem: initial_fields, cell_centers, validate_problem, is_valid_problem
using Pyrolysis.Problem: step_profile, linear_profile, gaussian_profile, tanh_profile
using Pyrolysis.Materials: n_components
using Pyrolysis.Meshes: Mesh, line_mesh, box_mesh
using Pyrolysis.States: mass_index, energy_index, history_index, has_chi
using Pyrolysis.Integration: setup

_gas() = GaseousComponent(:g, M = 0.03, c = 1000.0, k = 0.03, λ = 1e-5)
_material() = Material(;
    name = :inert,
    components = (SolidComponent(:s, ρ = 1000.0, c = 1000.0, k = 0.2), _gas()),
    reactions = (),
)
const _BCS = BoundaryConditionSet(thermal = (top = HeatFluxBC(1e4), bottom = AdiabaticBC()))

_T(z) = tanh_profile(
    z;
    z_center = 0.004,
    width = 0.001,
    value_left = 300.0,
    value_right = 600.0,
)
const _Ξ = [
    z -> gaussian_profile(
        z;
        z_center = 0.005,
        width = 0.002,
        amplitude = 100.0,
        baseline = 900.0,
    ),
    z ->
        linear_profile(z; z_min = 0.001, z_max = 0.009, value_min = 0.0, value_max = 2.0) + step_profile(z; z_interface = 0.006, value_left = 0.0, value_right = 0.5),
]

_problem(mesh; material = _material(), kw...) = PyrolysisProblem(;
    mesh,
    material,
    bc_set = _BCS,
    T_initial = _T,
    ξ_initial = _Ξ,
    tspan = (0.0, 1.0),
    kw...,
)

# The cell values the old mesh's initial-condition application wrote for these
# profiles on `create_stretched_mesh(0.01, 7, 1.15)`, recorded from it.
const _T_RECORDED = [
    300.81152818041426,
    335.68691935063157,
    540.4419862283842,
    596.1780671192503,
    599.70379216582,
    599.9682534715243,
    599.9954500045953,
]
const _Ξ_RECORDED = [
    914.1534953246537 960.6174308561601 998.864661595512 984.1467260597474 946.9403406177331 920.2056284327261 907.5338548885343
    0.011261264794142816 0.49970623810493164 0.9244409975056177 1.7937755708975183 2.114936069499171 2.394206068283218 2.5
]

@testset "problem construction: the arguments kept, the conditions and output times filled" begin
    mesh = create_uniform_mesh(0.01, 20)
    material = _material()
    keyword = PyrolysisProblem(; mesh, material, bc_set = _BCS, tspan = (0.0, 10.0))
    @test keyword.mesh === mesh
    @test keyword.material === material
    @test keyword.tspan == (0.0, 10.0)
    @test n_components(keyword.material) == 2
    positional = PyrolysisProblem(mesh, material, _BCS, (0.0, 10.0))
    @test positional.mesh === mesh && positional.material === material
    @test positional.tspan == (0.0, 10.0)

    # A set given only the top's thermal condition: an adiabatic bottom and
    # impermeable faces.
    partial = PyrolysisProblem(;
        mesh,
        material,
        bc_set = BoundaryConditionSet(thermal = (top = HeatFluxBC(5e4),)),
        tspan = (0.0, 10.0),
    )
    @test partial.bc_set.thermal.top isa HeatFluxBC
    @test partial.bc_set.thermal.bottom isa AdiabaticBC
    @test partial.bc_set.mass.top isa ImpermeableBC
    @test partial.bc_set.mass.bottom isa ImpermeableBC
    @test validate_problem(partial)
    @test is_valid_problem(partial)

    # Without output times, times from the start to the end of the span.
    auto = PyrolysisProblem(; mesh, material, bc_set = _BCS, tspan = (0.0, 100.0))
    @test length(auto.output_times) > 1 && issorted(auto.output_times)
    @test first(auto.output_times) == 0.0
    @test last(auto.output_times) == 100.0

    # Both display forms name the type, the material and its components.
    text = sprint(show, positional)
    @test occursin("PyrolysisProblem", text)
    @test occursin("inert", text)
    @test occursin("20 cells", text)
    full = sprint(show, MIME("text/plain"), positional)
    @test occursin("PyrolysisProblem", full)
    @test occursin("Material: inert", full)
    @test occursin("Components: 2", full)
end

@testset "initial fields: the profiles at the cell centres" begin
    mesh = create_stretched_mesh(0.01, 7, 1.15)
    T0, ξ0 = initial_fields(_problem(mesh), mesh)
    @test T0 == _T_RECORDED
    @test ξ0 == _Ξ_RECORDED
    z = [p[1] for p in mesh.positions]
    @test cell_centers(mesh) == [(z[i] + z[i+1]) / 2 for i = 1:7]
    # Radial meshes evaluate at the same midpoints.
    radial = create_stretched_mesh(0.01, 7, 1.15; coordinates = Cylindrical())
    @test initial_fields(_problem(radial), radial) == (T0, ξ0)

    # On more dimensions at the centroids, as `SVector`s.
    box = box_mesh([0.0, 1.0, 3.0], [0.0, 2.0]; fixed = true)
    q = _problem(
        box;
        T_initial = x -> 300.0 + x[1] + 10x[2],
        ξ_initial = [x -> 1000.0 * x[1], x -> 0.0],
    )
    T2, ξ2 = initial_fields(q, box)
    @test cell_centers(box) ≈ [SVector(0.5, 1.0), SVector(2.0, 1.0)]
    @test T2 ≈ [310.5, 312.0]
    @test ξ2 ≈ [500.0 2000.0; 0.0 0.0]
end

@testset "the set-up state holds the initial fields at the cell centres" begin
    mesh = create_stretched_mesh(0.01, 7, 1.15)
    z = [p[1] for p in mesh.positions]
    V = diff(z)
    problem = _problem(mesh)
    cp = setup(problem)
    layout = cp.def.layout
    @test cp.T0 == _T_RECORDED
    for i = 1:7
        @test cp.u0[energy_index(layout, i)] == _T_RECORDED[i]
        for j = 1:2
            @test cp.u0[mass_index(layout, i, j)] ≈ _Ξ_RECORDED[j, i] * V[i] rtol = 1e-12
        end
    end
    # The solve starts from it.
    sol = Pyrolysis.Integration.solve(problem)
    @test sol.retcode === :Success
    @test sol.T[:, 1] == _T_RECORDED
    @test sol.ξ[:, :, 1] ≈ _Ξ_RECORDED rtol = 1e-12

    # Under a lateral shrinkage law each cell's conversion history starts at
    # G = 0 and M its initial solid mass, the gas it holds excluded.
    shrinking = Material(;
        name = :shrinking,
        components = (SolidComponent(:s, ρ = 1000.0, c = 1000.0, k = 0.2), _gas()),
        reactions = (Reaction(:decomp, 1 => 2, A = 1e5, E = 5e4, h = 1e5),),
        lateral_shrinkage_law = χ -> 1.0 - 0.2χ,
    )
    cpχ = setup(_problem(mesh; material = shrinking); use_ale = true)
    layoutχ = cpχ.def.layout
    @test has_chi(layoutχ)
    for i = 1:7
        @test cpχ.u0[history_index(layoutχ, i, 1)] == 0
        @test cpχ.u0[history_index(layoutχ, i, 2)] ≈ _Ξ_RECORDED[1, i] * V[i] rtol = 1e-12
    end
end

@testset "validation on the native mesh: the same refusals" begin
    mesh = create_uniform_mesh(0.01, 20)
    good = _problem(mesh)
    @test validate_problem(good)
    @test is_valid_problem(good)
    @test validate_problem(good; strict = true)

    # The mesh: a position per node, finite.
    short = Mesh(mesh.topology, mesh.positions[1:(end-1)])
    @test_throws ErrorException validate_problem(_problem(short))
    nan = copy(mesh.positions)
    nan[3] = SVector(NaN)
    @test_throws ErrorException validate_problem(_problem(Mesh(mesh.topology, nan)))

    # The boundary tags: `:top` and `:bottom`, each with a thermal condition.
    x = collect(range(0.0, 0.01; length = 5))
    err = try
        validate_problem(_problem(line_mesh(x; top = :surface)))
        ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("[:top]", err)
    @test occursin(":surface", err)

    # A per-species ambient concentration with the wrong component count.
    bad_mass = BoundaryConditionSet(
        thermal = (top = HeatFluxBC(1e4), bottom = AdiabaticBC()),
        mass = (
            top = ConvectiveMassBC(h_m = 0.05, ξ_inf = (0.0, 0.1, 0.2)),
            bottom = ImpermeableBC(),
        ),
    )
    @test_throws ErrorException validate_problem(_problem(mesh; bc_set = bad_mass))

    # The time span, the first step and the initial-condition count.
    @test_throws ErrorException validate_problem(_problem(mesh; tspan = (10.0, 5.0)))
    @test_throws ErrorException validate_problem(_problem(mesh; dt_initial = 2.0))
    @test_throws ErrorException validate_problem(_problem(mesh; ξ_initial = [z -> 1.0]))

    # A permeable matrix with no pore space at the cell centres.
    dense = Material(;
        name = :wood,
        components = (
            SolidComponent(:wood; ρ = 500.0, c = 1500.0, k = 0.15, ε = 0.9, κ = 1e-14),
            _gas(),
        ),
        reactions = (Reaction(:pyrolysis, 1 => 2, A = 1e10, E = 150e3, h = -300e3),),
    )
    darcy(ξ; material = dense) = _problem(mesh; material, ξ_initial = [ξ, z -> 0.0])
    err = try
        validate_problem(darcy(z -> 500.0))
        ""
    catch e
        sprint(showerror, e)
    end
    # The message names the fix.
    @test occursin("ρ_intrinsic", err)
    @test occursin("background_gas_concentration", err)
    @test occursin("porosity", err)
    @test !is_valid_problem(darcy(z -> 500.0))
    @test validate_problem(darcy(z -> 0.0))
    # An intrinsic density gives the matrix its pore space; a region that
    # starts at it is still refused.
    porous = Material(;
        name = :wood,
        components = (
            SolidComponent(
                :wood;
                ρ = 500.0,
                c = 1500.0,
                k = 0.15,
                ε = 0.9,
                κ = 1e-14,
                ρ_intrinsic = 1500.0,
            ),
            _gas(),
        ),
        reactions = dense.reactions,
    )
    @test validate_problem(darcy(z -> 500.0; material = porous))
    @test_throws ErrorException validate_problem(
        darcy(z -> z < 0.005 ? 1500.0 : 500.0; material = porous),
    )
    # Permeable only in a char the initial state holds none of: no flow to
    # refuse in the dense cells.
    charring = Material(;
        name = :charring,
        components = (
            SolidComponent(:virgin; ρ = 500.0, c = 1500.0, k = 0.15, ε = 0.9),
            SolidComponent(
                :char;
                ρ = 100.0,
                c = 1100.0,
                k = 0.1,
                ε = 0.95,
                κ = 1e-12,
                ρ_intrinsic = 1950.0,
            ),
            _gas(),
        ),
        reactions = (
            Reaction(
                :pyrolysis,
                1 => (2, 3),
                A = 1e10,
                E = 150e3,
                h = -300e3,
                yields = (0.3, 0.7),
            ),
        ),
    )
    @test validate_problem(
        _problem(mesh; material = charring, ξ_initial = [z -> 500.0, z -> 0.0, z -> 0.0]),
    )

    # Strict: the initial conditions must evaluate at a cell centre.
    throwing = _problem(mesh; T_initial = z -> error("no temperature"))
    @test validate_problem(throwing)
    @test_throws ErrorException validate_problem(throwing; strict = true)

    # The solve validates by default.
    @test_throws ErrorException Pyrolysis.Integration.solve(
        _problem(mesh; tspan = (1.0, 0.0)),
    )
end

@testset "validation warns when Fickian gas transport lacks a condensed transfer coefficient" begin
    function ablating(; λ_solid = 0.0, κ = 0.0, reactive = true)
        solid = SolidComponent(
            :solid,
            ρ = 1000.0,
            c = 1000.0,
            k = 1.0,
            λ = λ_solid,
            κ = κ,
            ρ_intrinsic = 2000.0,
        )
        gas = GaseousComponent(:gas, M = 0.05, c = 1000.0, k = 0.02, λ = 1e-4)
        reactions = reactive ? (Reaction(:ablate, 1 => 2; A = 1e10, E = 1e5, h = 1e6),) : ()
        return Material(; name = :ablation, components = (solid, gas), reactions)
    end
    mesh = create_uniform_mesh(0.01, 10)
    deck(material) = _problem(mesh; material, ξ_initial = [z -> 1000.0, z -> 0.0])
    @test_logs (:warn, r"gas-transfer coefficient") validate_problem(deck(ablating()))
    # A condensed coefficient, a permeable matrix or no gas production keeps
    # validation silent.
    @test_logs validate_problem(deck(ablating(λ_solid = 1e-6)))
    @test_logs validate_problem(deck(ablating(κ = 1e-14)))
    @test_logs validate_problem(deck(ablating(reactive = false)))
end

@testset "unsupported stoichiometry cannot pass validation" begin
    for reactants in ((), ((1, 1.0), (2, 1.0)))
        rxn = Reaction(
            :unsupported,
            reactants,
            ((2, 1.0),),
            1.0,
            0.0,
            Pyrolysis.Materials.to_property(0.0),
            0.0,
            Inf,
        )
        @test isnan(Pyrolysis.Materials.mass_balance_error(rxn))
        @test_logs (:warn, r"unverifiable mass balance") @test !Pyrolysis.Materials.verify_stoichiometry(
            rxn,
        )
        material = Material(
            name = :unsupported,
            components = _material().components,
            reactions = (rxn,),
        )
        problem = _problem(create_uniform_mesh(0.01, 2); material)
        @test_throws r"requires exactly one reactant" validate_problem(problem)
    end
    balanced = Reaction(:balanced, 1 => 2; A = 1.0, E = 0.0, h = 0.0)
    @test Pyrolysis.Materials.verify_stoichiometry(balanced)
    unbalanced = Reaction(
        :unbalanced,
        1 => (1, 2);
        A = 1.0,
        E = 0.0,
        h = 0.0,
        yields = (0.2, 0.9),
        validate_mass = false,
    )
    @test_logs (:warn, r"mass balance error") @test !Pyrolysis.Materials.verify_stoichiometry(
        unbalanced,
    )
end
