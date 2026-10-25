# Reproduce guide figures with the documentation environment, without Plots:
#   julia --project=docs scripts/figures/generate.jl
# Optional --first / --second restrict generation to figures 1–8 / 9–15.
# --from-data redraws the tutorial from its committed CSVs without solving it.
using Pyrolysis
using Pyrolysis.Integration: setup
using Pyrolysis.Closures: EquidistributionCorrection, NoCorrection, activation
using Pyrolysis.States
using Pyrolysis.ColoredJacobian: dense_jacobian, missing_entries
using SparseArrays
using SHA

include("svg.jl")
include("first_batch.jl")
include("quickstart.jl")

function small_problem(; law = nothing)
    material = Material(
        name = :FigureMaterial,
        components = (
            SolidComponent(:virgin; ρ = 1000.0, c = 1400.0, k = 0.2, α = 0.2, λ = 1e-6),
            GaseousComponent(:gas; M = 0.03, c = 1100.0, k = 0.03, λ = 1e-5),
        ),
        reactions = (
            Reaction(:decomposition, :virgin => :gas; A = 5e9, E = 1.4e5, h = 8e5),
        ),
        lateral_shrinkage_law = law,
    )
    return PyrolysisProblem(
        mesh = create_stretched_mesh(0.01, 6, 1.1),
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = RadiativeFluxBC(3e4) + RadiativeBC(T_inf = 300.0),
                bottom = AdiabaticBC(),
            ),
            mass = (top = DiffusiveMassBC(0.02), bottom = ImpermeableBC()),
        ),
        T_initial = z -> 500.0 + 2e4z,
        ξ_initial = [z -> 1000.0, z -> 0.1],
        tspan = (0.0, 10.0),
    )
end

if !("--second" in ARGS)
    cp = setup(small_problem(); use_ale = true, correction = EquidistributionCorrection())
    r = collect(range(0, 0.5; length = 201))
    a = [activation(cp.def.closure.correction, x, 1.0) for x in r]
    @assert a[1] == 1 && a[end] == 0
    first_batch((r, a))
    sol, selected = "--from-data" in ARGS ? saved_quickstart_data() : quickstart_data()
    quickstart_figures(sol, selected)
end

if !("--first" in ARGS)
    include("second_batch.jl")
    second_batch()
end
