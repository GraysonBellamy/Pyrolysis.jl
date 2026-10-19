# The iteration matrix on a moving mesh, behind gate item 2.
#
# OrdinaryDiffEq forms the iteration matrix `W = J − M/(γ dt)` with
# `jacobian2W!`. With a mass matrix that is a multiple of the identity it
# copies `J` and shifts the diagonal, so `W` keeps `J`'s stored structure;
# with any other mass matrix (a `Diagonal` with zeros on the algebraic rows,
# as both cores have on a moving mesh) it broadcasts, and the sparse
# broadcast stores only the entries that compute to nonzero. `W`'s structure
# then follows the state, and KLU, whose symbolic analysis is reused only on
# an unchanged structure, analyses afresh at every step.
#
# 1. The structure: the stored entries of `J` and of `W` at every accepted
#    step of the new core on the `ale` deck (100 cells, solid gas-transfer
#    coefficient 1e-6 m²/s, reltol 1e-7), and how many distinct
#    factorizations the integrator made.
# 2. The cost: each core's default on the matched-accuracy decks of
#    `audit_decks.jl`, warm minimum of three, first with OrdinaryDiffEq's
#    `jacobian2W!` and then with a structure-keeping one substituted in this
#    process (copy `J`, shift the diagonal by `M/(γ dt)`; the same matrix, its
#    stored zeros kept). The substitution changes nothing in the package; it
#    measures what a structure-keeping `jacobian2W!` upstream would give both
#    cores. Run this script in its own process.
#
# Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/w_formation.jl
# Writes test/newcore/gates/measurements/w_formation.md.

ENV["PYROLYSIS_CORE"] = "old"
const _HARNESS =
    normpath(joinpath(@__DIR__, "..", "..", "..", "scripts", "mesh_phase0", "baselines"))
include(joinpath(_HARNESS, "jacobian_common.jl"))

using SparseArrays: SparseMatrixCSC, nnz
using LinearAlgebra: Diagonal
using OrdinaryDiffEq: DiscreteCallback

const OUT = joinpath(@__DIR__, "measurements")
const DIFFERENTIATION = Base.root_module(
    Base.PkgId(
        Base.UUID("4302a76b-040a-498a-8c04-15b101fed76b"),
        "OrdinaryDiffEqDifferentiation",
    ),
)
# (kind, cells, solid gas-transfer coefficient, the new core's matched reltol)
const DECKS = (
    (:eulerian, 100, 1e-6, 1e-7),
    (:ale, 25, 1e-6, 1e-7),
    (:ale, 100, 1e-6, 1e-7),
    (:ale_darcy, 25, nothing, 1e-6),
    (:ale_darcy, 100, nothing, 3e-7),
)

function structure()
    problem, kw = problem_for(:ale, 100; solid_λ = 1e-6)
    J_nnz, W_nnz, factorizations = Int[], Int[], UInt[]
    probe = DiscreteCallback(
        (u, t, integ) -> true,
        integ -> begin
            c = integ.cache
            push!(J_nnz, nnz(c.J))
            push!(W_nnz, nnz(c.W))
            push!(factorizations, objectid(c.linsolve.cacheval))
        end;
        save_positions = (false, false),
    )
    Pyrolysis.Integration.solve(
        problem;
        use_ale = true,
        kw.radiation_model,
        reltol = 1e-7,
        abstol = 1e-9,
        callback = probe,
    )
    return (;
        steps = length(J_nnz),
        J = extrema(J_nnz),
        W = extrema(W_nnz),
        W_distinct = length(unique(W_nnz)),
        factorizations = length(unique(factorizations)),
    )
end

function timings(kind, n, solid_λ, reltol)
    problem, kw = problem_for(kind, n; solid_λ)
    new() = Pyrolysis.Integration.solve(
        problem;
        use_ale = kw.use_ale,
        kw.radiation_model,
        reltol,
        abstol = reltol / 100,
    )
    old() = solve(problem; kw..., reltol = 1e-6, abstol = 1e-8)
    best(f) = (f(); minimum(@elapsed(f()) for _ = 1:3))
    return (; new = best(new), old = best(old), T = new().T)
end

include(joinpath(@__DIR__, "stamp.jl"))
load_start = Sys.loadavg()
started = time()
probe = structure()
println(probe)
before = Dict(d => timings(d...) for d in DECKS)
# The structure-keeping formation, in this process only.
@eval DIFFERENTIATION function jacobian2W!(
    W::$(SparseMatrixCSC),
    mass_matrix::$(Diagonal),
    dtgamma::Number,
    J::$(SparseMatrixCSC),
)::Nothing
    copyto!(W, J)
    invdtgamma = inv(dtgamma)
    @inbounds for i in axes(W, 1)
        m = mass_matrix.diag[i]
        iszero(m) || (W[i, i] -= m * invdtgamma)
    end
    return nothing
end
after = Dict(d => timings(d...) for d in DECKS)
kept = structure()
println(kept)
load_end = Sys.loadavg()

_g(x; d = 3) = string(round(x; sigdigits = d))
open(joinpath(OUT, "w_formation.md"), "w") do io
    println(io, "# The iteration matrix on a moving mesh\n")
    println(
        io,
        "Produced by `test/newcore/gates/w_formation.jl` at $(COMMIT_STAMP). Julia $(VERSION), ",
        "$(Threads.nthreads()) thread, $(Sys.CPU_THREADS) logical CPUs; load average ",
        "$(round.(load_start; digits = 1)) at the start and $(round.(load_end; digits = 1)) ",
        "at the end.\n",
    )
    println(
        io,
        "The new core on `ale` at 100 cells (solid λ 1e-6, reltol 1e-7): $(probe.steps) ",
        "accepted steps; `J` stores $(probe.J[1]) to $(probe.J[2]) entries, `W` ",
        "$(probe.W[1]) to $(probe.W[2]) ($(probe.W_distinct) distinct counts); ",
        "$(probe.factorizations) distinct KLU factorizations, each a fresh symbolic ",
        "analysis. With the structure-keeping formation `W` stores $(kept.W[1]) to ",
        "$(kept.W[2]) entries and the solve makes $(kept.factorizations) KLU ",
        "factorization, refactored in place at every other step.\n",
    )
    println(
        io,
        "Warm minimum of three [s]: the new core at its matched tolerance, the old core ",
        "at its default, with OrdinaryDiffEq's `jacobian2W!` (broadcast) and with a ",
        "structure-keeping one. ΔT: the largest temperature difference the substitution ",
        "makes to the new core's solution [K].\n",
    )
    println(
        io,
        "| deck | new, broadcast | new, kept | old, broadcast | old, kept | ratio, broadcast | ratio, kept | new kept / old broadcast | ΔT [K] |",
    )
    println(io, "|" * repeat("---|", 9))
    for d in DECKS
        b, a = before[d], after[d]
        deck = d[3] === nothing ? "$(d[1]) $(d[2])" : "$(d[1]) $(d[2]), solid λ $(d[3])"
        println(
            io,
            "| $(deck) | $(_g(b.new)) | $(_g(a.new)) | $(_g(b.old)) | $(_g(a.old)) | ",
            "$(_g(b.new / b.old)) | $(_g(a.new / a.old)) | $(_g(a.new / b.old)) | ",
            "$(_g(maximum(abs, a.T .- b.T))) |",
        )
    end
end
println(
    "wrote $(OUT); load $(load_start) -> $(load_end); $(round(Int, time() - started)) s",
)
