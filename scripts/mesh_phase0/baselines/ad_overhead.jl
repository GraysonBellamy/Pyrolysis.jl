# Phase 0 item 5: cost of one ForwardDiff dual pass through the current
# residual (chunk widths 4..32) against the Float64 residual and the
# structured Jacobian assembly, plus a real colored-AD Jacobian on the
# Eulerian decks checked against the structured exact Jacobian. The projected
# colored-AD cost ceil(K/C)·t_dual(C) uses the lifted model's color count K
# (both species stencils) with the CURRENT residual's dual cost; the lifted
# residual itself is not timed. Written to results/ad_overhead.csv. Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 scripts/mesh_phase0/baselines/ad_overhead.jl
const T_START = time()
include(joinpath(@__DIR__, "jacobian_common.jl"))
using ForwardDiff
using ForwardDiff: Dual, Partials

const OUT_CSV = joinpath(RESULTS, "ad_overhead.csv")
const DECKS = (:eulerian, :ale, :ale_darcy)
const SIZES = (100, 400)
const CHUNKS = (4, 8, 12, 16, 24, 32)
const T_EVAL = 1.0

struct PhaseTag end

# One-hot seed in slot `p` (0 ≤ p ≤ C; p = 0 gives zero partials).
@inline function seed_dual(::Type{Dual{Tag,Float64,C}}, v, p) where {Tag,C}
    return Dual{Tag}(v, Partials{C,Float64}(ntuple(m -> ifelse(m == p, 1.0, 0.0), Val(C))))
end

function dual_vectors(::Type{Tu}, u) where {Tu}
    u_d = Tu[seed_dual(Tu, u[i], mod1(i, ForwardDiff.npartials(Tu))) for i in eachindex(u)]
    return u_d, similar(u_d)
end

"""
Colored sparse Jacobian: seed each column's one-hot partial in the slot of its
color (chunked `C` colors per pass), run the dual residual, scatter the
partials into the stored entries of `S`.
"""
function colored_jacobian!(
    Jc,
    S,
    coloring,
    K,
    ::Type{Tu},
    u,
    u_d,
    du_d,
    def,
    t,
    ws_ad,
) where {Tu}
    C = ForwardDiff.npartials(Tu)
    passes = cld(K, C)
    rv = rowvals(S)
    for q = 1:passes
        for j in eachindex(u)
            k = coloring[j]
            chunk = (k - 1) ÷ C + 1
            u_d[j] = seed_dual(Tu, u[j], chunk == q ? (k - 1) % C + 1 : 0)
        end
        residual!(du_d, u_d, def, t, ws_ad)
        for j in eachindex(u)
            k = coloring[j]
            (k - 1) ÷ C + 1 == q || continue
            p = (k - 1) % C + 1
            for idx in nzrange(S, j)
                Jc.nzval[idx] = ForwardDiff.partials(du_d[rv[idx]], p)
            end
        end
    end
    return Jc
end

const HEADER = (
    "deck",
    "n",
    "NC",
    "state_len",
    "quantity",
    "C",
    "median_s",
    "min_s",
    "nevals",
    "repeats",
    "per_color_overhead",
    "lifted_colors_c0",
    "lifted_colors_c_active",
    "projected_c0_s",
    "projected_c_active_s",
    "colors",
    "passes",
    "max_abs_diff",
    "max_abs_ref",
    "note",
)

function emit!(
    io,
    deck,
    n,
    NC,
    N,
    quantity,
    C,
    tm;
    per_color = "",
    lifted_c0 = "",
    lifted_cact = "",
    projected_c0 = "",
    projected_cact = "",
    colors = "",
    passes = "",
    maxdiff = "",
    maxref = "",
    note = "",
)
    csv_row(
        io,
        deck,
        n,
        NC,
        N,
        quantity,
        C,
        tm.median,
        tm.min,
        tm.nevals,
        tm.repeats,
        per_color,
        lifted_c0,
        lifted_cact,
        projected_c0,
        projected_cact,
        colors,
        passes,
        maxdiff,
        maxref,
        note,
    )
    @printf(
        "%-10s n=%4d %-22s C=%3s median=%10.3e s min=%10.3e s %s\n",
        deck,
        n,
        quantity,
        string(C),
        tm.median,
        tm.min,
        note
    )
    flush(stdout)
end

function fail!(io, deck, n, NC, N, quantity, C, err)
    msg = replace(sprint(showerror, err), '\n' => ' ')
    csv_row(
        io,
        deck,
        n,
        NC,
        N,
        quantity,
        C,
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "FAILED: " * first(msg, 300),
    )
    println("FAILED $deck n=$n $quantity C=$C: ", first(msg, 300))
    flush(stdout)
end

open(OUT_CSV, "w") do io
    csv_row(io, HEADER...)
    for deck in DECKS, n in SIZES
        local def, ws, u, NC, n_cells
        try
            def, ws, u, NC, n_cells = build_deck(deck, n)
        catch err
            fail!(io, deck, n, "", "", "build", "", err)
            continue
        end
        hot_state!(u, n_cells)
        N = length(u)
        mode = def.mode
        chi = chi_mode(mode) isa WithChi
        bl = radiation_mode(mode) isa BeerLambert
        # Lifted color counts for both species stencils (c0: 3 cells; c_active: 5).
        K_c0 =
            color_pattern(lifted_pattern(n_cells, NC; chi, bl, species_stencil = 3)[1]).colors
        K_cact =
            color_pattern(lifted_pattern(n_cells, NC; chi, bl, species_stencil = 5)[1]).colors
        du = similar(u)

        # Float64 residual.
        t_res = time_op(() -> residual!(du, u, def, T_EVAL, ws))
        emit!(
            io,
            deck,
            n_cells,
            NC,
            N,
            "residual",
            "",
            t_res;
            lifted_c0 = K_c0,
            lifted_cact = K_cact,
        )

        # Structured assembly, exact and reduced.
        t_struct = Dict{Symbol,Any}()
        for (name, strategy) in
            (("structured_exact", DirectSolve()), ("structured_reduced", ApproxSolve()))
            try
                spec = JacobianSpec(Structured(; strategy), def, ws)
                J = Pyrolysis.Jacobian._alloc_target(spec, N)
                tm = time_op(() -> compute_jacobian!(J, spec, u, def, T_EVAL, ws))
                t_struct[Symbol(name)] = (spec, J, tm)
                emit!(io, deck, n_cells, NC, N, name, "", tm; note = "nnz=$(nnz(J))")
            catch err
                fail!(io, deck, n_cells, NC, N, name, "", err)
            end
        end

        # Dual passes.
        for C in CHUNKS
            Tu = Dual{PhaseTag,Float64,C}
            try
                ws_ad = build_ad_workspace(def, ws, Tu)
                ws_ad.I_radiation = ws.I_radiation
                u_d, du_d = dual_vectors(Tu, u)
                tm = time_op(() -> residual!(du_d, u_d, def, T_EVAL, ws_ad))
                per_color = (tm.median - t_res.median) / (C * t_res.median)
                emit!(
                    io,
                    deck,
                    n_cells,
                    NC,
                    N,
                    "dual",
                    C,
                    tm;
                    per_color = per_color,
                    lifted_c0 = K_c0,
                    lifted_cact = K_cact,
                    projected_c0 = cld(K_c0, C) * tm.median,
                    projected_cact = cld(K_cact, C) * tm.median,
                    note = @sprintf("dual/residual=%.2f", tm.median / t_res.median)
                )
            catch err
                fail!(io, deck, n_cells, NC, N, "dual", C, err)
            end
        end

        # Real colored Jacobian on the Eulerian deck (current :full pattern),
        # checked against the structured exact Jacobian at the same state.
        if deck === :eulerian && haskey(t_struct, :structured_exact)
            try
                spec, J_exact, _ = t_struct[:structured_exact]
                S = current_pattern(ws, :full)
                c = color_pattern(S)
                K = c.colors
                C = K <= 32 ? K : 16
                Tu = Dual{PhaseTag,Float64,C}
                ws_ad = build_ad_workspace(def, ws, Tu)
                ws_ad.I_radiation = ws.I_radiation
                u_d, du_d = dual_vectors(Tu, u)
                Jc = SparseMatrixCSC{Float64,Int}(
                    S.m,
                    S.n,
                    copy(S.colptr),
                    copy(S.rowval),
                    zeros(nnz(S)),
                )
                (Jc.colptr == J_exact.colptr && Jc.rowval == J_exact.rowval) ||
                    error("colored pattern differs from the structured target pattern")
                compute_jacobian!(J_exact, spec, u, def, T_EVAL, ws)
                colored_jacobian!(
                    Jc,
                    S,
                    c.coloring,
                    K,
                    Tu,
                    u,
                    u_d,
                    du_d,
                    def,
                    T_EVAL,
                    ws_ad,
                )
                maxdiff = maximum(abs, Jc.nzval .- J_exact.nzval)
                maxref = maximum(abs, J_exact.nzval)
                tm = time_op(
                    () -> colored_jacobian!(
                        Jc,
                        S,
                        c.coloring,
                        K,
                        Tu,
                        u,
                        u_d,
                        du_d,
                        def,
                        T_EVAL,
                        ws_ad,
                    ),
                )
                emit!(
                    io,
                    deck,
                    n_cells,
                    NC,
                    N,
                    "colored_ad_current_full",
                    C,
                    tm;
                    colors = K,
                    passes = cld(K, C),
                    maxdiff = maxdiff,
                    maxref = maxref,
                    note = @sprintf(
                        "colored/structured_exact=%.2f; ordering=%s",
                        tm.median / t_struct[:structured_exact][3].median,
                        c.ordering
                    )
                )
            catch err
                fail!(io, deck, n_cells, NC, N, "colored_ad_current_full", "", err)
            end
        end
    end
end
@printf("wrote %s in %.1f s\n", OUT_CSV, time() - T_START)
