# What deactivation costs the PMMA cone of `examples/06` at order zero, against
# what the mesh costs: the deck without deactivation at 72 and 144 cells, and with
# it at 72. A run passes when deactivation moves the integrated gas generation and
# the mass by no more than the refinement does. The pointwise peak is recorded
# only: at order zero the rate is a staircase of single cells, and that figure
# reads the time between two runs' steps.
#
# Run from the repository root, once per form:
#   PYROLYSIS_ENERGY_STATE=temperature julia --project=test --startup-file=no --threads=1 \
#       test/gates/burnout_order_zero.jl
# Writes test/gates/measurements/burnout_order_zero_<form>.md.

using Pyrolysis
using Pyrolysis.Closures: NoCorrection
using Printf
include(joinpath(@__DIR__, "stamp.jl"))

const FORM = Symbol(get(ENV, "PYROLYSIS_ENERGY_STATE", "temperature"))
const T_END = 2400.0

function cone(N)
    lines = readlines(joinpath(@__DIR__, "..", "..", "examples", "06_cone_calorimeter.jl"))
    a = findfirst(startswith("pmma_cp = CallableProperty("), lines)
    b = findfirst(startswith("solution = solve("), lines)
    kept = String[]
    skip = false
    for line in lines[a:(b-1)]
        if skip
            startswith(line, ")") && (skip = false)
            continue
        end
        if startswith(line, "println(")
            skip = !endswith(rstrip(line), ")")
            continue
        end
        startswith(line, "flush(") || push!(kept, line)
    end
    m = Module(:cone)
    Core.eval(m, :(using Pyrolysis))
    Core.eval(m, :(const SMOKE = false))
    Core.eval(m, Meta.parseall(join(kept, "\n")))
    p = Core.eval(m, :(problem))
    r = only(p.material.reactions)
    material = Material(
        name = p.material.name,
        components = p.material.components,
        reactions = (
            Reaction(
                r.name,
                ((r.reactants[1][1], 0.0),),
                r.products,
                r.A,
                r.E,
                r.heat,
                r.T_min,
                r.T_max;
                gate_width = r.gate_width,
                validate_mass = false,
                basis = r.basis,
            ),
        ),
        mixing = p.material.mixing,
        lateral_shrinkage_law = p.material.lateral_shrinkage_law,
        depletion_limiter = p.material.depletion_limiter,
        gas_viscosity = p.material.gas_viscosity,
        gas_transport = p.material.gas_transport,
    )
    return PyrolysisProblem(
        mesh = create_uniform_mesh(5.8e-3, N),
        material = material,
        bc_set = p.bc_set,
        T_initial = p.T_initial,
        ξ_initial = p.ξ_initial,
        tspan = (0.0, T_END),
        dt_initial = p.dt_initial,
        output_times = collect(0.0:1.0:T_END),
    )
end

run_deck(problem; deactivation) = solve(
    problem;
    abstol = 1e-11,
    reltol = 1e-9,
    pressure_tolerance = 1e-4,
    dtmin = 1e-14,
    use_ale = true,
    min_thickness = 0.01e-3,
    correction = NoCorrection(),
    deactivation,
    progress = false,
    energy_state = FORM,
)

# Over 30 s to the end: the largest and the integrated difference of the gas
# generation, over the peak and over the mass lost, and the largest mass difference.
function difference(ref, test)
    window = findall(t -> t >= 30.0, ref.t)
    t = ref.t[window]
    r = ref.extras.gas_generation_rate[window]
    d = abs.(test.extras.gas_generation_rate[window] .- r)
    M = ref.extras.total_mass
    lost = M[1] - M[window[end]]
    integral = sum((d[j] + d[j+1]) / 2 * (t[j+1] - t[j]) for j = 1:(length(t)-1))
    return (
        peak = maximum(d) / maximum(r),
        integral = integral / lost,
        mass = maximum(abs.(test.extras.total_mass[window] .- M[window])) / M[1],
    )
end

# The decks come from evaluating the example, so they are built first: a solve in
# the same statement could not call the laws that evaluation defines.
deck_72, deck_144 = cone(72), cone(144)
coarse = run_deck(deck_72; deactivation = false)
fine = run_deck(deck_144; deactivation = false)
burned = run_deck(deck_72; deactivation = true)
all(completed, (coarse, fine, burned)) || error("a run did not reach the end")
cost = difference(coarse, burned)
mesh = difference(coarse, fine)
pass = cost.integral <= mesh.integral && cost.mass <= mesh.mass

lines = [
    "# Order zero against the mesh, $(FORM) form",
    "",
    COMMIT_STAMP,
    "",
    "| difference | peak | integrated | mass |",
    "|---|---:|---:|---:|",
    @sprintf(
        "| deactivation, 72 cells | %.3f %% | %.3f %% | %.3e |",
        100cost.peak,
        100cost.integral,
        cost.mass
    ),
    @sprintf(
        "| 72 to 144 cells | %.3f %% | %.3f %% | %.3e |",
        100mesh.peak,
        100mesh.integral,
        mesh.mass
    ),
    "",
    "Deactivation: $(length(burned.events)) events. Result: $(pass ? "pass" : "fail").",
]
out = joinpath(@__DIR__, "measurements", "burnout_order_zero_$(FORM).md")
write(out, join(lines, "\n") * "\n")
println(join(lines, "\n"))
pass || exit(1)
