# Phase 0 item 4: column color counts of the current production Jacobian
# patterns against the lifted-pattern model, written to
# results/jacobian_colors.csv. Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 scripts/mesh_phase0/baselines/jacobian_colors.jl
const T_START = time()
include(joinpath(@__DIR__, "jacobian_common.jl"))
include(joinpath(ROOT, "test", "verification", "casedefs.jl"))
include(joinpath(ROOT, "test", "verification", "cases", "V-BO-01.jl"))

const OUT_CSV = joinpath(RESULTS, "jacobian_colors.csv")
const DECKS = (:eulerian, :ale, :ale_darcy, :ale_darcy_withchi, :vbo01)
const SIZES = (25, 100, 400)
const MODEL_NC = (2, 3, 5)

function emit!(io, deck, n, NC, pattern, stencil, S, nloc, formula, note)
    t0 = time()
    c = color_pattern(S)
    bw = pattern_bandwidth(S, nloc)
    csv_row(
        io,
        deck,
        n,
        NC,
        pattern,
        stencil,
        size(S, 1),
        c.nnz,
        bw,
        c.max_row_degree,
        c.colors,
        c.ordering,
        formula,
        note,
    )
    @printf(
        "%-22s n=%4d NC=%d %-20s %-8s rows=%5d nnz=%8d bw=%5d maxdeg=%5d colors=%5d formula=%3s (%s) [%.1fs]\n",
        deck,
        n,
        NC,
        pattern,
        stencil,
        size(S, 1),
        c.nnz,
        bw,
        c.max_row_degree,
        c.colors,
        string(formula),
        c.ordering,
        time() - t0
    )
    flush(stdout)
    return c
end

# Lifted model for both species stencils; returns Dict(stencil => colors).
function emit_lifted!(io, deck, n, NC; chi, bl)
    out = Dict{Int,Int}()
    for s in STENCILS
        S, L = lifted_pattern(n, NC; chi, bl, species_stencil = s)
        h = (s - 1) ÷ 2
        note =
            "cell-major interleaved; cell rows read every unknown of cells i-$(h)..i+$(h), " *
            "nodes z_(i-$(h))..z_(i+$(h+1)), w_i, w_(i+1)" *
            (bl ? ", I_i, I_(i+1)" : "") *
            (
                chi ?
                ", A; chi as the inventories G,M; dense A,theta_A rows excluded from the coloring (columns kept); bandwidth over local columns" :
                ""
            ) *
            "; radial has identical locality (areas and volumes are local functions of the node positions), no separate radial rows"
        c = emit!(
            io,
            deck,
            n,
            NC,
            "lifted_planar",
            stencil_label(s),
            S,
            ncols_local(L),
            lifted_color_formula(NC, s; chi, bl),
            note,
        )
        out[s] = c.colors
    end
    return out
end

open(OUT_CSV, "w") do io
    csv_row(
        io,
        "deck",
        "n",
        "NC",
        "pattern",
        "stencil",
        "rows",
        "nnz",
        "bandwidth",
        "max_row_degree",
        "greedy_colors",
        "ordering",
        "model_formula",
        "note",
    )
    for deck in DECKS
        sizes = deck === :vbo01 ? (0,) : SIZES
        for n in sizes
            local def, ws, u, NC, n_cells
            try
                def, ws, u, NC, n_cells = build_deck(deck, n)
            catch err
                msg = sprint(showerror, err)
                csv_row(
                    io,
                    deck,
                    n,
                    "",
                    "BUILD_FAILED",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    replace(msg, '\n' => ' '),
                )
                println("BUILD FAILED for $deck n=$n: ", msg)
                flush(stdout)
                continue
            end
            mode = def.mode
            is_ale = mesh_mode(mode) isa ALE
            chi = chi_mode(mode) isa WithChi
            bl = radiation_mode(mode) isa BeerLambert
            modenote = "mode=" * replace(string(typeof(mode)), "Pyrolysis.Problem." => "")
            S_full = current_pattern(ws, :full)
            emit!(
                io,
                deck,
                n_cells,
                NC,
                "current_full",
                "",
                S_full,
                size(S_full, 2),
                "",
                modenote,
            )
            if is_ale
                S_red = current_pattern(ws, :ale_reduced)
                emit!(
                    io,
                    deck,
                    n_cells,
                    NC,
                    "current_ale_reduced",
                    "",
                    S_red,
                    size(S_red, 2),
                    "",
                    modenote,
                )
            end
            emit_lifted!(io, deck, n_cells, NC; chi, bl)
        end
    end
    # Lifted-pattern model alone: NC and n ladders for the three physics
    # combinations (surface absorption; Beer–Lambert; Beer–Lambert + χ).
    for (deck, chi, bl) in (
        ("model_sa", false, false),
        ("model_bl", false, true),
        ("model_bl_chi", true, true),
    )
        for NC in MODEL_NC, n in SIZES
            emit_lifted!(io, deck, n, NC; chi, bl)
        end
    end
end
@printf("wrote %s in %.1f s\n", OUT_CSV, time() - T_START)
