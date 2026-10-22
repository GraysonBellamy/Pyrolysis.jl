# The old core's bit-level fingerprint (plan §9.3): an edit to the old core's
# source (§2.3: additive, behind the old defaults) must leave its kinetics and
# its solves bit-identical. Hashed are 3,072 reaction rates and heats through
# `Physics.evaluate_reactions!` (reaction orders 0, 0.5, 1 and 2; depletion
# thresholds 1e-3, 0.1 and 10 kg/m³, the limiter on and off; 16 solid
# concentrations from −1 to 1e3 kg/m³ at four temperatures) and the saved
# temperatures, concentrations and times of eight old-core solves (V-RX-01,
# 02a, 02c, 02d, 03, 04a, V-CPL-01 and V-CPL-04 at their CI resolutions), each
# against its digest on 17bc477, the tree before the finite-volume core's
# old-core edits. The digests are Julia's `hash`, which depends on the Julia
# version; the reference digests are Julia 1.13.0's.
#
# Run from the repository root before and after an old-core edit:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/oldcore_fingerprint.jl
# Writes test/newcore/gates/measurements/oldcore_fingerprint.md.

ENV["PYROLYSIS_CORE"] = "old"
using Pyrolysis
include(joinpath(@__DIR__, "..", "..", "verification", "casedefs.jl"))

const OUT = joinpath(@__DIR__, "measurements")
mkpath(OUT)

# The digests on 17bc477.
const REFERENCE_RATES = 0xc8fea394cd668f5a
const REFERENCE_SOLVES = (
    "V-RX-01" => UInt64(189106276520408973),
    "V-RX-02a" => UInt64(2020822377329231978),
    "V-RX-02c" => UInt64(15262768053745139010),
    "V-RX-02d" => UInt64(6136562261384834374),
    "V-RX-03" => UInt64(15903990052688888789),
    "V-RX-04a" => UInt64(13431261678381534594),
    "V-CPL-01" => UInt64(11013119082032119509),
    "V-CPL-04" => UInt64(2838787475738871859),
)

# The hash of the rates and heats over the grid, and their number.
function rate_fingerprint()
    vals = Float64[]
    for n in (0.0, 0.5, 1.0, 2.0), threshold in (1e-3, 1e-1, 10.0), enabled in (true, false)
        solid = SolidComponent(:s, ρ = 1000.0, c = 1000.0, k = 0.2)
        gas = GaseousComponent(:g, M = 0.03, c = 1000.0, k = 0.03, λ = 1e-5)
        mat = Material(;
            name = :fp,
            components = (solid, gas),
            reactions = (Reaction(:r, 1 => 2; A = 1e8, E = 1e5, h = 1e5, n = n),),
            depletion_limiter = Pyrolysis.Materials.TanhLimiter(; threshold, enabled),
        )
        rates, heats = zeros(1), zeros(1)
        for ξ in (
                -1.0,
                -1e-3,
                0.0,
                1e-6,
                1e-4,
                1e-3,
                1e-2,
                0.05,
                0.1,
                0.5,
                1.0,
                2.0,
                5.0,
                10.0,
                100.0,
                1000.0,
            ),
            T in (300.0, 450.0, 700.0, 1200.0)

            Pyrolysis.Physics.evaluate_reactions!(rates, heats, (ξ, 1.0), T, mat)
            push!(vals, rates[1], heats[1])
        end
    end
    return hash(vals), length(vals)
end

include(joinpath(@__DIR__, "stamp.jl"))
started = time()
rate_digest, n_rates = rate_fingerprint()
rows = [(
    item = "reaction rates and heats ($(n_rates) values)",
    reference = REFERENCE_RATES,
    digest = rate_digest,
    detail = "",
)]
cases = load_verification_cases()
for (id, reference) in REFERENCE_SOLVES
    case = only(c for c in cases if c.id == id)
    problem = Base.invokelatest(case.build, case.ci.n_cells)
    sol = Base.invokelatest(_case_solve, case, problem; case.solve_kwargs...)
    push!(
        rows,
        (
            item = id,
            reference,
            digest = hash((sol.T, sol.ξ, sol.t)),
            detail = "$(sol.retcode), $(length(sol.t)) saved times",
        ),
    )
    println(id, ": ", sol.retcode, " ", rows[end].digest == reference ? "equal" : "DIFFERS")
    flush(stdout)
end
identical = all(r -> r.digest == r.reference, rows)

_hex(h) = "0x" * string(h; base = 16, pad = 16)
open(joinpath(OUT, "oldcore_fingerprint.md"), "w") do io
    println(io, "# The old core's bit-level fingerprint\n")
    println(
        io,
        "Produced by `test/newcore/gates/oldcore_fingerprint.jl` at $(COMMIT_STAMP). ",
        "Julia $(VERSION). The hash of $(n_rates) reaction rates and heats through ",
        "`Physics.evaluate_reactions!` (orders 0, 0.5, 1, 2; depletion thresholds 1e-3, ",
        "0.1 and 10 kg/m³, the limiter on and off; 16 solid concentrations from −1 to ",
        "1e3 kg/m³ at 300, 450, 700 and 1200 K) and `hash((T, ξ, t))` of eight old-core ",
        "solves at their CI resolutions, against the same digests on 17bc477, the tree ",
        "before the finite-volume core's old-core edits (Julia 1.13.0).\n",
    )
    println(io, "| item | 17bc477 | this tree | equal | solve |")
    println(io, "|---|---|---|---|---|")
    for r in rows
        println(
            io,
            "| $(r.item) | `$(_hex(r.reference))` | `$(_hex(r.digest))` | ",
            "$(r.digest == r.reference) | $(r.detail) |",
        )
    end
    println(
        io,
        identical ? "\nEvery digest equals its reference: the old core is bit-identical." :
        "\nA digest differs from its reference: the old core is not bit-identical.",
    )
end
println(
    identical ? "every digest equals its reference" : "a digest differs",
    "; wrote $(OUT); $(round(Int, time() - started)) s",
)
