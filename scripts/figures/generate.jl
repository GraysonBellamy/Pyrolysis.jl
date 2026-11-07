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
            SolidComponent(
                :virgin;
                ρ = 1000.0,
                φ = 0.3,
                c = 1400.0,
                k = 0.2,
                α = 0.2,
                λ = 1e-6,
            ),
            GaseousComponent(:gas; M = 0.03, c = 1100.0, k = 0.03),
        ),
        reactions = (
            Reaction(:decomposition, :virgin => :gas; A = 5e9, E = 1.4e5, h = 8e5),
        ),
        lateral_shrinkage_law = law,
    )
    return PyrolysisProblem(
        mesh = stretched_mesh(0.01, 6, 1.1),
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = RadiativeFluxBC(3e4) + RadiativeBC(T_inf = 300.0),
                bottom = AdiabaticBC(),
            ),
            mass = (
                top = AmbientGasBC(GasReservoir(P = 0.0, T = 300.0)),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> 500.0 + 2e4z,
        ξ_initial = [z -> 1000.0, z -> 0.1],
        tspan = (0.0, 10.0),
    )
end

# The worked slab of Technical Reference §8.8 at its consistent start: 10 mm,
# α = 100 m⁻¹, ε = 0.9, a diffuse 50 kW/m² heater and surroundings at 300 K,
# a mirror back. Returns the inward depth of every face and its two fluxes over εq.
function radiation_slab()
    slab = SolidComponent(:slab, ρ = 1000.0, c = 1500.0, k = 0.2, ε = 0.9, α = 0.1)
    problem = PyrolysisProblem(
        mesh = uniform_mesh(0.01, 10),
        material = Material(name = :slab, components = (slab,), reactions = ()),
        bc_set = BoundaryConditionSet(
            thermal = (
                top = RadiativeFluxBC(50e3) + RadiativeBC(T_inf = 300.0),
                bottom = AdiabaticBC(),
            ),
        ),
        T_initial = z -> 300.0,
        ξ_initial = [z -> 1000.0],
        tspan = (0.0, 1.0),
    )
    cp = setup(problem; radiation_model = RADIATIVE_TRANSFER)
    layout, u = cp.def.layout, cp.u0
    faces = Pyrolysis.Meshes.n_nodes(layout):-1:1
    flux(k) = [u[radiation_index(layout, f, k)] / (0.9 * 50e3) for f in faces]
    return (collect(range(0, 1; length = length(faces))), flux(1), flux(2))
end

if !("--second" in ARGS)
    cp = setup(small_problem(); use_ale = true, correction = EquidistributionCorrection())
    r = collect(range(0, 0.5; length = 201))
    a = [activation(cp.def.closure.correction, x, 1.0) for x in r]
    @assert a[1] == 1 && a[end] == 0
    first_batch((r, a), radiation_slab())
    sol, selected = "--from-data" in ARGS ? saved_quickstart_data() : quickstart_data()
    quickstart_figures(sol, selected)
end

if !("--first" in ARGS)
    include("second_batch.jl")
    second_batch()
end
