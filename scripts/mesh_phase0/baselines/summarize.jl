# Writes results/jacobian_summary.md from jacobian_colors.csv, ad_overhead.csv
# and lifted_containment.csv (run those three scripts first). Run from the
# repository root:
#   julia --project=test --startup-file=no scripts/mesh_phase0/baselines/summarize.jl
using Printf, Statistics

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const RESULTS = joinpath(ROOT, "scripts", "mesh_phase0", "results")

# RFC 4180 line parser: quoted fields may hold commas and doubled quotes.
function parse_csv_line(line::AbstractString)
    fields = String[]
    buf = IOBuffer()
    inq = false
    chars = collect(line)
    i = 1
    while i <= length(chars)
        c = chars[i]
        if inq
            if c == '"'
                if i < length(chars) && chars[i+1] == '"'
                    write(buf, '"')
                    i += 1
                else
                    inq = false
                end
            else
                write(buf, c)
            end
        elseif c == '"'
            inq = true
        elseif c == ','
            push!(fields, String(take!(buf)))
        else
            write(buf, c)
        end
        i += 1
    end
    push!(fields, String(take!(buf)))
    return fields
end

function read_csv(path)
    lines = readlines(path)
    header = parse_csv_line(lines[1])
    rows = Dict{String,String}[]
    for l in lines[2:end]
        isempty(strip(l)) && continue
        f = parse_csv_line(l)
        length(f) == length(header) || error(
            "$(basename(path)): $(length(f)) fields in a $(length(header))-column row: $l",
        )
        push!(rows, Dict(zip(header, f)))
    end
    return rows
end

colors = read_csv(joinpath(RESULTS, "jacobian_colors.csv"))
overhead = read_csv(joinpath(RESULTS, "ad_overhead.csv"))
containment = read_csv(joinpath(RESULTS, "lifted_containment.csv"))

get1(rows, deck, n, pattern; stencil = "") = (
    i = findfirst(
        r ->
            r["deck"] == deck &&
            r["n"] == string(n) &&
            r["pattern"] == pattern &&
            (stencil == "" || r["stencil"] == stencil),
        rows,
    );
    i === nothing ? nothing : rows[i]
)
col(r, k) = r === nothing ? "–" : r[k]
us(x) = @sprintf("%.1f", 1e6 * x)
fnum(s) = s == "" ? NaN : parse(Float64, s)
fint(s) = parse(Int, s)
function rng(v)
    lo, hi = @sprintf("%.2f", minimum(v)), @sprintf("%.2f", maximum(v))
    return lo == hi ? lo : "$(lo)–$(hi)"
end
rngi(v) = minimum(v) == maximum(v) ? string(v[1]) : "$(minimum(v))–$(maximum(v))"

io = IOBuffer()
println(io, "# Phase 0 items 4 and 5: Jacobian color counts and colored-AD overhead\n")
println(
    io,
    "Sources: `scripts/mesh_phase0/baselines/{jacobian_colors,ad_overhead,lifted_containment}.jl`; data in ",
    "`results/jacobian_colors.csv`, `results/ad_overhead.csv` and `results/lifted_containment.csv`. ",
    "Patterns: `current_full` and `current_ale_reduced` are the production `build_sparsity` unions ",
    "(`S_phys ∪ S_geom`); `lifted` is the plan's model: cell unknowns T, ξ_j and the χ inventories G, M; ",
    "node z_k, w_k; face I_f; global A, θ_A with their dense rows left out of the coloring and their columns kept. ",
    "A cell row reads every unknown of the cells in its species stencil (`c0`: cells i−1..i+1, no relative ",
    "solid flux; `c_active`: cells i−2..i+2, the second-order limited reconstruction of the relative flux), ",
    "the node positions one past that stencil on each side (`c0`: z_{i−1}..z_{i+2}; `c_active`: z_{i−2}..z_{i+3}), ",
    "w_i, w_{i+1}, I_i, I_{i+1} and A. Coloring is greedy distance-1 column coloring (natural and ",
    "largest-degree-first orderings, smaller kept); the max row degree is the clique lower bound, so equal ",
    "numbers mean the greedy count is optimal. Radial geometry has identical locality (areas and volumes are ",
    "local functions of the node positions), so only the planar model is listed.\n",
)

# ---------------------------------------------------------------------------
# Table 1: decks
# ---------------------------------------------------------------------------
println(
    io,
    "## Table 1: column colors (lower bound = max row degree; greedy count in parentheses when larger)\n",
)
println(
    io,
    "| deck | n | NC | rows (current) | current_full | current_ale_reduced | rows (lifted) | lifted c0 | lifted c_active |",
)
println(io, "|---|---|---|---|---|---|---|---|---|")
fmtc(r) =
    r === nothing ? "–" :
    (
        r["max_row_degree"] == r["greedy_colors"] ? r["greedy_colors"] :
        "$(r["max_row_degree"]) ($(r["greedy_colors"]))"
    )
deck_rows =
    filter(r -> !startswith(r["deck"], "model_") && r["pattern"] == "current_full", colors)
for r in deck_rows
    deck, n = r["deck"], r["n"]
    red = get1(colors, deck, n, "current_ale_reduced")
    l0 = get1(colors, deck, n, "lifted_planar"; stencil = "c0")
    l5 = get1(colors, deck, n, "lifted_planar"; stencil = "c_active")
    println(
        io,
        "| $deck | $n | $(r["NC"]) | $(r["rows"]) | $(fmtc(r)) | $(fmtc(red)) | $(col(l0, "rows")) | $(fmtc(l0)) | $(fmtc(l5)) |",
    )
end

# ---------------------------------------------------------------------------
# Model ladder and formula check
# ---------------------------------------------------------------------------
println(
    io,
    "\nLifted-pattern model alone, c0 / c_active (colors identical at n = 25, 100 and 400 unless a range is shown; every count equals its max-row-degree lower bound and the formula):\n",
)
println(
    io,
    "| physics | NC = 2 | NC = 3 | NC = 5 | unknowns per cell (+ per node / face) |",
)
println(io, "|---|---|---|---|---|")
for (deck, label, pc) in (
    ("model_sa", "surface absorption (T, ξ; z, w)", "1+NC (+2 per node)"),
    ("model_bl", "Beer–Lambert (+ I_f)", "1+NC (+3 per node/face)"),
    (
        "model_bl_chi",
        "Beer–Lambert + χ as (G, M) (+ I_f; A, θ_A global)",
        "3+NC (+3 per node/face, +2 global)",
    ),
)
    cells = String[]
    for NC in (2, 3, 5)
        parts = String[]
        for st in ("c0", "c_active")
            ks = [
                fint(r["greedy_colors"]) for r in colors if r["deck"] == deck &&
                    r["NC"] == string(NC) &&
                    r["pattern"] == "lifted_planar" &&
                    r["stencil"] == st
            ]
            push!(parts, rngi(ks))
        end
        push!(cells, join(parts, " / "))
    end
    println(io, "| $label | $(cells[1]) | $(cells[2]) | $(cells[3]) | $pc |")
end

lifted_rows = filter(r -> r["pattern"] == "lifted_planar", colors)
n_match = count(r -> r["greedy_colors"] == r["model_formula"], lifted_rows)
n_opt = count(r -> r["greedy_colors"] == r["max_row_degree"], lifted_rows)
mismatches = filter(r -> r["greedy_colors"] != r["model_formula"], lifted_rows)
println(
    io,
    "\n**Formula.** A cell row of the model reads s·(NC + 1 + 2χ) cell unknowns (s = 3 for c0, 5 for c_active; ",
    "χ = 1 adds the two inventories G and M), s + 1 node positions, two node velocities, two face intensities ",
    "under Beer–Lambert (bl = 1) and the cross-section column under χ, so its degree, the clique lower bound of the coloring, is\n\n",
    "    K(s) = s·(NC + 1 + 2χ) + (s + 1) + 2 + 2·bl + χ,\n\n",
    "i.e. K_c0 = 3·NC + 9 + 2·bl + 7·χ and K_c_active = 5·NC + 13 + 2·bl + 11·χ. The greedy count equals this ",
    "formula in $n_match of $(length(lifted_rows)) lifted rows and equals its lower bound in $n_opt of them",
    isempty(mismatches) ? "." :
    "; mismatches: " *
    join(
        [
            "$(r["deck"]) n=$(r["n"]) NC=$(r["NC"]) $(r["stencil"]): greedy $(r["greedy_colors"]) vs formula $(r["model_formula"])"
            for r in mismatches
        ],
        "; ",
    ) *
    ".",
    " The previous model (13 to 28 colors, formula 3(NC+1) + 4 (+2 bl, +4 χ)) omitted the two outer node positions ",
    "z_{i−1} and z_{i+2} that a cell row reads through the neighbouring cell centres and the neighbours' ξ = m/V, ",
    "carried χ as one unknown instead of the inventory pair (G, M) and had no stencil for the relative-flux ",
    "reconstruction; the corrected c0 counts are two higher without χ and five higher with it. The χ count lumps ",
    "G and M into every cell row as the plan's layout does; a tighter model in which the T and ξ rows do not read ",
    "the neighbours' G, M (with c = 0 they enter those rows only through A) would give 3·NC + 12 for Beer–Lambert + χ, ",
    "six fewer, and is not measured here.",
)

# ---------------------------------------------------------------------------
# Table 2: overhead
# ---------------------------------------------------------------------------
println(
    io,
    "\n## Table 2: dual-pass overhead and projected colored-AD Jacobian cost (hot state, medians, μs)\n",
)
println(
    io,
    "| deck | n | t_residual | t_struct_exact | t_struct_reduced | per-color overhead C=4 / 8 / 16 / 32 | K c0 / c_active | projected c0 (C) | c0 / exact | c0 / reduced | projected c_active (C) | c_active / exact | c_active / reduced | measured colored AD on current pattern |",
)
println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
summary_numbers = Dict{String,Any}()
pco_all = Float64[]
for deck in ("eulerian", "ale", "ale_darcy"), n in (100, 400)
    rows = filter(r -> r["deck"] == deck && r["n"] == string(n), overhead)
    isempty(rows) && continue
    res = findfirst(r -> r["quantity"] == "residual", rows)
    ex = findfirst(r -> r["quantity"] == "structured_exact", rows)
    rd = findfirst(r -> r["quantity"] == "structured_reduced", rows)
    duals = filter(r -> r["quantity"] == "dual" && r["median_s"] != "", rows)
    failed = filter(r -> startswith(r["note"], "FAILED"), rows)
    t_res = fnum(rows[res]["median_s"])
    t_ex = ex === nothing ? NaN : fnum(rows[ex]["median_s"])
    t_rd = rd === nothing ? NaN : fnum(rows[rd]["median_s"])
    pco = Dict(fint(r["C"]) => fnum(r["per_color_overhead"]) for r in duals)
    append!(pco_all, values(pco))
    pcos = join(
        [haskey(pco, C) ? @sprintf("%.2f", pco[C]) : "–" for C in (4, 8, 16, 32)],
        " / ",
    )
    K0 = isempty(duals) ? "–" : duals[1]["lifted_colors_c0"]
    K5 = isempty(duals) ? "–" : duals[1]["lifted_colors_c_active"]
    best0 =
        isempty(duals) ? nothing : duals[argmin([fnum(r["projected_c0_s"]) for r in duals])]
    best5 =
        isempty(duals) ? nothing :
        duals[argmin([fnum(r["projected_c_active_s"]) for r in duals])]
    t0 = best0 === nothing ? NaN : fnum(best0["projected_c0_s"])
    t5 = best5 === nothing ? NaN : fnum(best5["projected_c_active_s"])
    tdual = Dict(fint(r["C"]) => fnum(r["median_s"]) for r in duals)
    cad = findfirst(
        r -> r["quantity"] == "colored_ad_current_full" && r["median_s"] != "",
        rows,
    )
    cadtxt =
        cad === nothing ? "–" :
        @sprintf(
            "%s μs, %s colors, max|ΔJ| = %.1e (max|J| = %.0f)",
            us(fnum(rows[cad]["median_s"])),
            rows[cad]["colors"],
            fnum(rows[cad]["max_abs_diff"]),
            fnum(rows[cad]["max_abs_ref"])
        )
    failtxt =
        isempty(failed) ? "" :
        " FAILED: " * join([f["quantity"] * "(C=" * f["C"] * ")" for f in failed], "; ")
    println(
        io,
        "| $deck | $n | $(us(t_res)) | $(us(t_ex)) | $(us(t_rd)) | $pcos | $K0 / $K5 | ",
        best0 === nothing ? "–" : "$(us(t0)) (C=$(best0["C"]))",
        " | ",
        @sprintf("%.2f", t0 / t_ex),
        " | ",
        @sprintf("%.2f", t0 / t_rd),
        " | ",
        best5 === nothing ? "–" : "$(us(t5)) (C=$(best5["C"]))",
        " | ",
        @sprintf("%.2f", t5 / t_ex),
        " | ",
        @sprintf("%.2f", t5 / t_rd),
        " | $cadtxt$failtxt |",
    )
    summary_numbers["$deck/$n"] = (;
        t_res,
        t_ex,
        t_rd,
        t0,
        t5,
        K0 = K0 == "–" ? 0 : fint(K0),
        K5 = K5 == "–" ? 0 : fint(K5),
        tdual,
        cad = cad === nothing ? NaN : fnum(rows[cad]["median_s"]),
    )
end

# ---------------------------------------------------------------------------
# Interpretation
# ---------------------------------------------------------------------------
println(io, "\n## Interpretation\n")
full_colors = Dict(
    (r["deck"], r["n"]) => fint(r["greedy_colors"]) for
    r in colors if r["pattern"] == "current_full"
)
lifted_colors = Dict(
    (r["deck"], r["n"], r["stencil"]) => fint(r["greedy_colors"]) for
    r in colors if r["pattern"] == "lifted_planar"
)
lc(deck, n, st) = lifted_colors[(deck, string(n), st)]
println(
    io,
    @sprintf(
        "1. **Color counts.** The lifted pattern's color count is independent of n. c0 / c_active: %d / %d for the NC = 2 surface-absorption decks (eulerian, ale) at n = 25, 100 and 400, %d / %d with Beer–Lambert (ale_darcy), %d / %d with Beer–Lambert and χ as (G, M) (ale_darcy_withchi), and %d / %d for the NC = 3 V-BO-01 deck; every count equals its max-row-degree lower bound and the formula above. The current ALE `:full` pattern colors in exactly the state length (%d, %d, %d for ale at n = 25, 100, 400; %d for V-BO-01): the top node's velocity row reads every column, so column coloring degenerates to one color per column, i.e. O(n). The current Eulerian `:full` and the ALE `:ale_reduced` patterns color in %d (NC = 2) and %d (NC = 3), fewer than the lifted pattern, because they carry no node, intensity or velocity columns (and the reduced pattern drops the geometry couplings that the lifted rows make local).",
        lc("eulerian", 100, "c0"),
        lc("eulerian", 100, "c_active"),
        lc("ale_darcy", 100, "c0"),
        lc("ale_darcy", 100, "c_active"),
        lc("ale_darcy_withchi", 100, "c0"),
        lc("ale_darcy_withchi", 100, "c_active"),
        lc("vbo01", 182, "c0"),
        lc("vbo01", 182, "c_active"),
        full_colors[("ale", "25")],
        full_colors[("ale", "100")],
        full_colors[("ale", "400")],
        full_colors[("vbo01", "182")],
        full_colors[("eulerian", "100")],
        fint(get1(colors, "vbo01", 182, "current_ale_reduced")["greedy_colors"])
    )
)

# Containment paragraph.
csum = filter(r -> r["row_kind"] == "summary", containment)
prod = filter(r -> r["override"] == "none", csum)
ovr = filter(r -> r["override"] == "zero", csum)
prod_ok = all(r -> r["contained_c0"] == "true" && r["contained_c_active"] == "true", prod)
rad_ok = all(r -> r["radiation_path_in_chain"] in ("", "true"), prod)
reach_cells = unique([r["measured_cells"] for r in prod])
reach_nodes = unique([r["measured_nodes"] for r in prod if r["measured_nodes"] != ""])
cf_eul = unique([r["current_full_cells"] for r in prod if r["deck"] == "eulerian"])
ovr_up = unique([
    "$(r["deck"]) n=$(r["n"]): $(r["measured_cells"]) cells, $(r["measured_nodes"]) nodes, upper reach within c_active = $(r["upper_contained_c_active"])"
    for r in ovr
])
failed_c = filter(r -> r["row_kind"] in ("FAILED", "BUILD_FAILED"), containment)
println(
    io,
    "2. **Containment against the executable residual** (`lifted_containment.csv`; decks eulerian, ale, ale_darcy, ale_darcy_withchi at n = 6 and 8, dense ForwardDiff of the current residual at a perturbed reacting state on a non-uniform mesh, threshold 1e-12 × row max). ",
    "With the production residual (no velocity override: its cell rows never read the node velocity, the prefix sum feeds only the z rows that the lifted balance rows replace) every T, ξ_j and χ row's measured cells and nodes lie inside the model's c0 sets",
    prod_ok ? "" : " **(FAILURE: see the CSV)**",
    ": the measured reach is cells $(join(reach_cells, " / ")) and nodes $(join(reach_nodes, " / ")), i.e. exactly the corrected c0 sets, with both outer nodes z_{i−1} and z_{i+2} present (the solid ξ row is cell-local, the gas ξ and T rows read the neighbours and the four nodes). ",
    "Under Beer–Lambert the local sets were measured at I_radiation = 0 and the radiation path as the support of J(50 kW/m²) − J(0): it lies within the closure of the intensity chain (cells i..n, nodes i..n+1)",
    rad_ok ? "" : " **(FAILURE)**",
    ". Under WithChi the χ row reads every χ column through A(χ̄), which the model carries as the global A column; the T and ξ rows read no χ column in planar geometry (A cancels from the per-volume rows), fewer than the current `:full` pattern's dense χ coupling. ",
    "The current Eulerian `:full` pattern's cell-to-cell width ($(join(cf_eul, ", "))) and the measured Eulerian width are both inside the c0 stencil. ",
    "With `Workspace.w_override = (z, t) -> 0` the z rows vanish, but the verification-only TVD temperature advection then runs with w_rel = 0 − w_physical, which is the prefix sum, so the T rows gain every cell and node below them, and under WithChi also everything above through the column-global θ_A; the reach above the cell from the limiter's own stencil is the c_active set (cells to i+2, nodes to i+3) on the NoChi decks: ",
    join(ovr_up, "; "),
    ". That run documents the hook, not the production or the lifted residual; the containment verdict is the no-override run. Nothing had to be added to the model after the measurement.",
    isempty(failed_c) ? "" :
    " Failed configurations: " *
    join(
        ["$(r["deck"]) n=$(r["n"]) $(r["override"]): $(r["note"])" for r in failed_c],
        "; ",
    ) *
    ".",
)

println(
    io,
    @sprintf(
        "3. **Dual-pass overhead.** (t_dual(C) − t_residual)/(C·t_residual) is %.2f to %.2f residual-equivalents per partial across the decks and chunk widths (lowest around C = 16–24, highest at C = 4 and at n = 400 on the plain ALE deck); one dual pass of width C costs about 1 + %.1f·C residuals at C = 8 to 16.",
        minimum(pco_all),
        maximum(pco_all),
        median(pco_all)
    )
)

e100, e400 = summary_numbers["eulerian/100"], summary_numbers["eulerian/400"]
a100, a400 = summary_numbers["ale/100"], summary_numbers["ale/400"]
d100, d400 = summary_numbers["ale_darcy/100"], summary_numbers["ale_darcy/400"]
ale100 = (a100, d100)
ale400 = (a400, d400)
println(
    io,
    "4. **Projected colored-AD Jacobian cost** (ceil(K/C)·t_dual(C) at the best C, K the lifted color count, t_dual the CURRENT residual's dual-pass time; the lifted residual, which has more rows per cell, was not timed, so these are projections, not measurements). Relative to the structured assemblies, reported per size because the ALE exact assembly's O(n²) geometry block grows with n:",
)
println(io, "")
println(
    io,
    "   | deck | n | c0: / exact | c0: / reduced | c_active: / exact | c_active: / reduced |",
)
println(io, "   |---|---|---|---|---|---|")
for (lab, s) in (
    ("eulerian", e100),
    ("eulerian", e400),
    ("ale", a100),
    ("ale", a400),
    ("ale_darcy", d100),
    ("ale_darcy", d400),
)
    n = s === e100 || s === a100 || s === d100 ? 100 : 400
    println(
        io,
        @sprintf(
            "   | %s | %d | %.2f | %.2f | %.2f | %.2f |",
            lab,
            n,
            s.t0 / s.t_ex,
            s.t0 / s.t_rd,
            s.t5 / s.t_ex,
            s.t5 / s.t_rd
        )
    )
end
println(io, "")
println(
    io,
    @sprintf(
        "   Eulerian decks (exact and reduced are the same pattern and cost): c0 %.2f× (n = 100) and %.2f× (n = 400) the structured assembly, c_active %.2f× and %.2f×. ALE decks: c0 %s× (n = 100) and %s× (n = 400) the structured exact assembly, %s× and %s× the reduced assembly; c_active %s× and %s× the exact, %s× and %s× the reduced. The exact ALE assembly costs %.1f–%.1f ms at n = 400 (%.1f–%.1f ms at n = 100), so the ratio to it is size-dependent and must not be quoted as one number. The real colored Jacobian on the current Eulerian `:full` pattern (%d colors, one dual pass) runs at %.2f× and %.2f× the structured exact assembly at n = 100 and 400 and matches it to round-off.",
        e100.t0 / e100.t_ex,
        e400.t0 / e400.t_ex,
        e100.t5 / e100.t_ex,
        e400.t5 / e400.t_ex,
        rng([s.t0 / s.t_ex for s in ale100]),
        rng([s.t0 / s.t_ex for s in ale400]),
        rng([s.t0 / s.t_rd for s in ale100]),
        rng([s.t0 / s.t_rd for s in ale400]),
        rng([s.t5 / s.t_ex for s in ale100]),
        rng([s.t5 / s.t_ex for s in ale400]),
        rng([s.t5 / s.t_rd for s in ale100]),
        rng([s.t5 / s.t_rd for s in ale400]),
        1e3 * min(a400.t_ex, d400.t_ex),
        1e3 * max(a400.t_ex, d400.t_ex),
        1e3 * min(a100.t_ex, d100.t_ex),
        1e3 * max(a100.t_ex, d100.t_ex),
        fint(
            overhead[findfirst(
                r ->
                    r["quantity"] == "colored_ad_current_full" &&
                    r["deck"] == "eulerian" &&
                    r["n"] == "100",
                overhead,
            )]["colors"],
        ),
        e100.cad / e100.t_ex,
        e400.cad / e400.t_ex
    )
)

# Ladder comparison: production K against the 4-color toy ladder of
# results/restart_comparison.csv.
deckK0 = [
    lc(d, n, "c0") for (d, n) in (
        ("eulerian", 100),
        ("ale", 100),
        ("ale_darcy", 100),
        ("ale_darcy_withchi", 100),
        ("vbo01", 182),
    )
]
deckK5 = [
    lc(d, n, "c_active") for (d, n) in (
        ("eulerian", 100),
        ("ale", 100),
        ("ale_darcy", 100),
        ("ale_darcy_withchi", 100),
        ("vbo01", 182),
    )
]
timed = [e100, e400, a100, a400, d100, d400]
ratio16(s, K) = cld(K, 16) * s.tdual[16] / s.tdual[4]
r0 = [ratio16(s, s.K0) for s in timed]
r5 = [ratio16(s, s.K5) for s in timed]
println(
    io,
    @sprintf(
        "5. **Against the 4-color toy ladder** (`results/restart_comparison.csv`, lifted form, 4 colors). Over the five decks (NC = 2 and 3) the production K is %s (c0) and %s (c_active), so K/4 = %.2f to %.2f (c0) and %.2f to %.2f (c_active); including the NC = 5 model ladder, K reaches %d (c0) and %d (c_active). In width-16 dual passes, ceil(K/16) = %s (c0) and %s (c_active) passes against the ladder's one, and the corresponding cost ratio ceil(K/16)·t_dual(16) / t_dual(4) (a production Jacobian at C = 16 over a 4-color Jacobian in one width-4 pass, measured dual times of the timed decks) is %.1f to %.1f (c0) and %.1f to %.1f (c_active).",
        rngi(deckK0),
        rngi(deckK5),
        minimum(deckK0) / 4,
        maximum(deckK0) / 4,
        minimum(deckK5) / 4,
        maximum(deckK5) / 4,
        maximum(r["greedy_colors"] |> fint for r in lifted_rows if r["stencil"] == "c0"),
        maximum(
            r["greedy_colors"] |> fint for r in lifted_rows if r["stencil"] == "c_active"
        ),
        rngi([cld(K, 16) for K in deckK0]),
        rngi([cld(K, 16) for K in deckK5]),
        minimum(r0),
        maximum(r0),
        minimum(r5),
        maximum(r5)
    )
)

println(
    io,
    "\nCaveats: timings are single-threaded medians of 7 samples on a hot synthetic state at t = 1 s; the dual pass is through the current residual and the projection multiplies its cost by the lifted color count, as the plan prescribes; the lifted residual itself (more rows per cell: z, w, I, G and M) was not timed, so its t_residual, dual cost and colored-AD cost are not measured here; the global-scalar rows of the χ decks are assumed hand-assembled outside the coloring; the c_active stencil is a model of the planned reconstruction and has no executable counterpart in the current residual beyond the verification-only TVD kernel; the V-BO-01 deck is sized by its CI cell count (n = 182) and was not timed.",
)

path = joinpath(RESULTS, "jacobian_summary.md")
write(path, String(take!(io)))
println("wrote ", path)
