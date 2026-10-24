# Phase 0 item 4 check: do the lifted-pattern model's dependency sets contain
# what the executable residual actually reads? Small ALE decks (6 and 8
# cells) for :ale, :ale_darcy (Beer–Lambert) and :ale_darcy_withchi, and the
# Eulerian deck, are differentiated with dense ForwardDiff at a perturbed,
# reacting state on a non-uniform mesh. For every cell row (T, ξ_j, χ) the
# cells, nodes and χ columns with |∂/∂u| > 1e-12·(row max) are extracted and
# compared with the model's sets (containment, not equality). Written to
# results/lifted_containment.csv. Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 scripts/mesh_phase0/baselines/lifted_containment.jl
#
# Two mesh-velocity configurations are measured:
#   override=none  the production residual. Its cell rows never read the node
#                  velocity (the prefix sum feeds only the z rows, which the
#                  lifted model replaces by local balance rows), so the
#                  measured sets are the local physics and geometry alone.
#   override=zero  `Workspace.w_override = (z, t) -> 0`. This zeroes the z
#                  rows, but the verification-only TVD temperature advection
#                  then fires with w_rel = 0 − w_physical, which is the prefix
#                  sum, so the T rows gain the cells and nodes BELOW them. The
#                  reach ABOVE the cell (second neighbour and its far node)
#                  comes from the limiter's own stencil and is reported
#                  against the c_active model. The production verdict is the
#                  override=none run.
# Under Beer–Lambert the current residual's T rows read every cell above
# through the optical depth; the lifted model carries that through the face
# intensities I_f. The local sets are therefore measured with I_radiation = 0
# and the radiation path separately as the support of J(I₀) − J(0), which
# must lie within the closure of the intensity chain (cells i..n, nodes
# i..n+1).
const T_START = time()
include(joinpath(@__DIR__, "jacobian_common.jl"))
using Pyrolysis.Jacobian: DenseAD_ForwardDiff
using Pyrolysis.Problem: temperature_range, concentration_range, z_range, chi_range
using Pyrolysis.Materials: n_components

const OUT_CSV = joinpath(RESULTS, "lifted_containment.csv")
const T_EVAL = 1.0
const REL_TOL = 1e-12
const N_CELLS = (6, 8)
const I0 = 50e3

# Perturbed, reacting state: hot_state! on T and ξ, a non-uniform monotone
# mesh (interior nodes displaced by 0.15 h·sin), and a χ profile when the deck
# carries one.
function perturbed_state!(u, layout, n)
    hot_state!(u, n)
    zr = z_range(layout)
    if !isempty(zr)
        z = view(u, zr)
        h = (z[end] - z[1]) / n
        for k = 2:n
            z[k] += 0.15h * sin(2π * (k - 1) / n)
        end
    end
    for (m, idx) in enumerate(chi_range(layout))
        u[idx] = 0.1 + 0.3m / n
    end
    return u
end

# Dense ForwardDiff Jacobian of the current residual with the workspace's
# override and incident intensity as set; the Dual workspace built by the spec
# is checked to carry the same override.
function dense_jacobian(def, ws, u; w_override, I_radiation)
    ws.w_override = w_override
    I_radiation === nothing || (ws.I_radiation = I_radiation)
    spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
    spec.cache.ws_ad.w_override === ws.w_override ||
        error("the Dual workspace does not carry w_override")
    J = Pyrolysis.Jacobian._alloc_target(spec, length(u))
    compute_jacobian!(J, spec, u, def, T_EVAL, ws)
    ws.w_override = nothing
    return J
end

cell_cols(layout, m) =
    (m, (first(concentration_range(layout, j)) + m - 1 for j = 1:n_components(layout))...)

# Cells, nodes and χ columns a row of `J` reads above the threshold, plus the
# count of nonzero entries at or below it.
function measured_dependencies(J, layout, row)
    n = layout.n_cells
    rowmax = maximum(abs, view(J, row, :))
    thr = REL_TOL * rowmax
    cells = Int[]
    nodes = Int[]
    chis = Int[]
    sub = 0
    for m = 1:n
        hit = false
        for col in cell_cols(layout, m)
            v = abs(J[row, col])
            v > thr && (hit = true)
            0 < v <= thr && (sub += 1)
        end
        hit && push!(cells, m)
    end
    for (k, col) in enumerate(z_range(layout))
        v = abs(J[row, col])
        v > thr && push!(nodes, k)
        0 < v <= thr && (sub += 1)
    end
    for (m, col) in enumerate(chi_range(layout))
        v = abs(J[row, col])
        v > thr && push!(chis, m)
        0 < v <= thr && (sub += 1)
    end
    return (; cells, nodes, chis, sub, rowmax)
end

# Cells and nodes a row of a Bool pattern over the current layout reads.
function pattern_dependencies(S, layout, row)
    n = layout.n_cells
    cells = [m for m = 1:n if any(c -> S[row, c], cell_cols(layout, m))]
    nodes = [k for (k, col) in enumerate(z_range(layout)) if S[row, col]]
    return (; cells, nodes)
end

# Compact run encoding of a sorted index set: "2-4", "1-3 6-8", "" when empty.
function runs(v::AbstractVector{Int})
    isempty(v) && return ""
    out = String[]
    a = b = v[1]
    for x in v[2:end]
        if x == b + 1
            b = x
        else
            push!(out, a == b ? string(a) : "$a-$b")
            a = b = x
        end
    end
    push!(out, a == b ? string(a) : "$a-$b")
    return join(out, ' ')
end

# Offsets of a set relative to cell i: (left, right) = (i − min, max − i).
offsets(v, i) = isempty(v) ? (0, 0) : (i - minimum(v), maximum(v) - i)
# Node offsets relative to the cell's own nodes i and i+1.
node_offsets(v, i) = isempty(v) ? (0, 0) : (i - minimum(v), maximum(v) - (i + 1))

const HEADER = (
    "deck",
    "n",
    "NC",
    "override",
    "I_radiation",
    "row_kind",
    "cell",
    "measured_cells",
    "measured_nodes",
    "measured_chi_cols",
    "model_cells_c0",
    "model_nodes_c0",
    "model_cells_c_active",
    "model_nodes_c_active",
    "current_full_cells",
    "current_full_nodes",
    "contained_c0",
    "contained_c_active",
    "upper_contained_c_active",
    "radiation_path_cells",
    "radiation_path_nodes",
    "radiation_path_in_chain",
    "subthreshold_nonzeros",
    "rowmax",
    "note",
)

# Rows of the current layout with their model counterparts: (kind, cell, row,
# model_row(L)).
function cell_rows(layout, L)
    n = layout.n_cells
    NC = n_components(layout)
    out = Tuple{String,Int,Int,Function}[]
    for i = 1:n
        push!(out, ("T", i, i, L -> lT(L, i)))
        for j = 1:NC
            push!(
                out,
                (
                    "xi$j",
                    i,
                    first(concentration_range(layout, j)) + i - 1,
                    L -> lξ(L, j, i),
                ),
            )
        end
        if !isempty(chi_range(layout))
            push!(out, ("chi", i, first(chi_range(layout)) + i - 1, L -> lG(L, i)))
        end
    end
    return out
end

function run_configuration!(io, deck, n, def, ws, u, J_local, J_rad, override, I_rad, note)
    layout = ws.layout
    NC = n_components(layout)
    chi = chi_mode(def.mode) isa WithChi
    bl = radiation_mode(def.mode) isa BeerLambert
    is_ale = mesh_mode(def.mode) isa ALE
    S0, L0 = lifted_pattern(n, NC; chi, bl, species_stencil = 3)
    S5, L5 = lifted_pattern(n, NC; chi, bl, species_stencil = 5)
    S_full = current_pattern(ws, :full)
    all_c0 = true
    all_c5 = true
    all_up5 = true
    all_rad = true
    nfail_c0 = 0
    maxoff = Dict(:cl => 0, :cr => 0, :nl => 0, :nr => 0, :fl => 0, :fr => 0)
    for (kind, i, row, mrow) in cell_rows(layout, L0)
        m = measured_dependencies(J_local, layout, row)
        d0 = model_row_dependencies(S0, L0, mrow(L0))
        d5 = model_row_dependencies(S5, L5, mrow(L5))
        pf = pattern_dependencies(S_full, layout, row)
        chi_ok0 = isempty(m.chis) || d0.A
        chi_ok5 = isempty(m.chis) || d5.A
        c0 = issubset(m.cells, d0.cells) && issubset(m.nodes, d0.nodes) && chi_ok0
        c5 = issubset(m.cells, d5.cells) && issubset(m.nodes, d5.nodes) && chi_ok5
        up5 =
            (isempty(m.cells) || maximum(m.cells) <= maximum(d5.cells)) &&
            (isempty(m.nodes) || maximum(m.nodes) <= maximum(d5.nodes)) &&
            chi_ok5
        all_c0 &= c0
        all_c5 &= c5
        all_up5 &= up5
        nfail_c0 += !c0
        cl, cr = offsets(m.cells, i)
        nl, nr = node_offsets(m.nodes, i)
        fl, fr = offsets(pf.cells, i)
        maxoff[:cl] = max(maxoff[:cl], cl)
        maxoff[:cr] = max(maxoff[:cr], cr)
        maxoff[:nl] = max(maxoff[:nl], nl)
        maxoff[:nr] = max(maxoff[:nr], nr)
        maxoff[:fl] = max(maxoff[:fl], fl)
        maxoff[:fr] = max(maxoff[:fr], fr)
        rad_cells = ""
        rad_nodes = ""
        rad_ok = ""
        if J_rad !== nothing
            D = J_rad .- J_local
            r = measured_dependencies(D, layout, row)
            rad_cells = runs(r.cells)
            rad_nodes = runs(r.nodes)
            ok =
                all(c -> i <= c <= n, r.cells) &&
                all(k -> i <= k <= n + 1, r.nodes) &&
                (isempty(r.chis) || d0.A)
            all_rad &= ok
            rad_ok = string(ok)
        end
        csv_row(
            io,
            deck,
            n,
            NC,
            override,
            I_rad,
            kind,
            i,
            runs(m.cells),
            runs(m.nodes),
            runs(m.chis),
            runs(d0.cells),
            runs(d0.nodes),
            runs(d5.cells),
            runs(d5.nodes),
            runs(pf.cells),
            runs(pf.nodes),
            c0,
            c5,
            up5,
            rad_cells,
            rad_nodes,
            rad_ok,
            m.sub,
            @sprintf("%.3e", m.rowmax),
            "",
        )
    end
    summary_note =
        note * (
            is_ale ?
            @sprintf(
                "; measured reach: cells i-%d..i+%d, nodes i-%d..i+%d; current :full pattern cells i-%d..i+%d",
                maxoff[:cl],
                maxoff[:cr],
                maxoff[:nl],
                maxoff[:nr] + 1,
                maxoff[:fl],
                maxoff[:fr]
            ) :
            @sprintf(
                "; measured reach: cells i-%d..i+%d; current :full pattern cells i-%d..i+%d",
                maxoff[:cl],
                maxoff[:cr],
                maxoff[:fl],
                maxoff[:fr]
            )
        )
    csv_row(
        io,
        deck,
        n,
        NC,
        override,
        I_rad,
        "summary",
        "",
        "i-$(maxoff[:cl])..i+$(maxoff[:cr])",
        is_ale ? "i-$(maxoff[:nl])..i+$(maxoff[:nr]+1)" : "",
        "",
        "i-1..i+1",
        is_ale ? "i-1..i+2" : "",
        "i-2..i+2",
        is_ale ? "i-2..i+3" : "",
        "i-$(maxoff[:fl])..i+$(maxoff[:fr])",
        "",
        all_c0,
        all_c5,
        all_up5,
        "",
        "",
        J_rad === nothing ? "" : string(all_rad),
        "",
        "",
        summary_note,
    )
    @printf(
        "%-20s n=%d override=%-5s I=%-7s contained_c0=%s contained_c_active=%s upper_c_active=%s rad_in_chain=%s (%d/%d rows fail c0) %s\n",
        deck,
        n,
        override,
        string(I_rad),
        all_c0,
        all_c5,
        all_up5,
        J_rad === nothing ? "-" : string(all_rad),
        nfail_c0,
        length(cell_rows(layout, L0)),
        summary_note
    )
    flush(stdout)
    return (; all_c0, all_c5, all_up5, all_rad)
end

open(OUT_CSV, "w") do io
    csv_row(io, HEADER...)
    for deck in (:eulerian, :ale, :ale_darcy, :ale_darcy_withchi), n in N_CELLS
        local def, ws, u, NC, n_cells
        try
            def, ws, u, NC, n_cells = build_deck(deck, n)
        catch err
            msg = replace(sprint(showerror, err), '\n' => ' ')
            csv_row(io, deck, n, "", "", "", "BUILD_FAILED", ntuple(_ -> "", 18)..., msg)
            println("BUILD FAILED for $deck n=$n: ", msg)
            flush(stdout)
            continue
        end
        perturbed_state!(u, ws.layout, n_cells)
        bl = radiation_mode(def.mode) isa BeerLambert
        is_ale = mesh_mode(def.mode) isa ALE
        I_default = ws.I_radiation
        overrides =
            is_ale ? (("none", nothing), ("zero", (z, t) -> 0.0)) : (("none", nothing),)
        for (olabel, ovr) in overrides
            try
                if bl
                    J0 = dense_jacobian(def, ws, u; w_override = ovr, I_radiation = 0.0)
                    J1 = dense_jacobian(def, ws, u; w_override = ovr, I_radiation = I0)
                    ws.I_radiation = I_default
                    run_configuration!(
                        io,
                        deck,
                        n_cells,
                        def,
                        ws,
                        u,
                        J0,
                        J1,
                        olabel,
                        0.0,
                        "local sets at I_radiation=0; radiation path = support of J($(I0)) - J(0)",
                    )
                else
                    J = dense_jacobian(def, ws, u; w_override = ovr, I_radiation = nothing)
                    run_configuration!(
                        io,
                        deck,
                        n_cells,
                        def,
                        ws,
                        u,
                        J,
                        nothing,
                        olabel,
                        I_default,
                        is_ale ? "surface absorption" : "Eulerian (no node columns)",
                    )
                end
            catch err
                msg = replace(sprint(showerror, err), '\n' => ' ')
                csv_row(
                    io,
                    deck,
                    n_cells,
                    NC,
                    olabel,
                    "",
                    "FAILED",
                    ntuple(_ -> "", 18)...,
                    first(msg, 400),
                )
                println("FAILED $deck n=$n override=$olabel: ", first(msg, 400))
                flush(stdout)
            end
        end
    end
end
@printf("wrote %s in %.1f s\n", OUT_CSV, time() - T_START)
