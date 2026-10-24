# Phase 0 items 4 and 5: Jacobian color counts and colored-AD overhead

Sources: `scripts/mesh_phase0/baselines/{jacobian_colors,ad_overhead,lifted_containment}.jl`; data in `results/jacobian_colors.csv`, `results/ad_overhead.csv` and `results/lifted_containment.csv`. Patterns: `current_full` and `current_ale_reduced` are the production `build_sparsity` unions (`S_phys ∪ S_geom`); `lifted` is the plan's model: cell unknowns T, ξ_j and the χ inventories G, M; node z_k, w_k; face I_f; global A, θ_A with their dense rows left out of the coloring and their columns kept. A cell row reads every unknown of the cells in its species stencil (`c0`: cells i−1..i+1, no relative solid flux; `c_active`: cells i−2..i+2, the second-order limited reconstruction of the relative flux), the node positions one past that stencil on each side (`c0`: z_{i−1}..z_{i+2}; `c_active`: z_{i−2}..z_{i+3}), w_i, w_{i+1}, I_i, I_{i+1} and A. Coloring is greedy distance-1 column coloring (natural and largest-degree-first orderings, smaller kept); the max row degree is the clique lower bound, so equal numbers mean the greedy count is optimal. Radial geometry has identical locality (areas and volumes are local functions of the node positions), so only the planar model is listed.

## Table 1: column colors (lower bound = max row degree; greedy count in parentheses when larger)

| deck | n | NC | rows (current) | current_full | current_ale_reduced | rows (lifted) | lifted c0 | lifted c_active |
|---|---|---|---|---|---|---|---|---|
| eulerian | 25 | 2 | 75 | 9 | – | 127 | 15 | 23 |
| eulerian | 100 | 2 | 300 | 9 | – | 502 | 15 | 23 |
| eulerian | 400 | 2 | 1200 | 9 | – | 2002 | 15 | 23 |
| ale | 25 | 2 | 101 | 101 | 9 | 127 | 15 | 23 |
| ale | 100 | 2 | 401 | 401 | 9 | 502 | 15 | 23 |
| ale | 400 | 2 | 1601 | 1601 | 9 | 2002 | 15 | 23 |
| ale_darcy | 25 | 2 | 101 | 101 | 9 | 153 | 17 | 25 |
| ale_darcy | 100 | 2 | 401 | 401 | 9 | 603 | 17 | 25 |
| ale_darcy | 400 | 2 | 1601 | 1601 | 9 | 2403 | 17 | 25 |
| ale_darcy_withchi | 25 | 2 | 126 | 126 | 9 | 205 | 24 | 36 |
| ale_darcy_withchi | 100 | 2 | 501 | 501 | 9 | 805 | 24 | 36 |
| ale_darcy_withchi | 400 | 2 | 2001 | 2001 | 9 | 3205 | 24 | 36 |
| vbo01 | 182 | 3 | 911 | 911 | 12 | 1277 | 20 | 30 |

Lifted-pattern model alone, c0 / c_active (colors identical at n = 25, 100 and 400 unless a range is shown; every count equals its max-row-degree lower bound and the formula):

| physics | NC = 2 | NC = 3 | NC = 5 | unknowns per cell (+ per node / face) |
|---|---|---|---|---|
| surface absorption (T, ξ; z, w) | 15 / 23 | 18 / 28 | 24 / 38 | 1+NC (+2 per node) |
| Beer–Lambert (+ I_f) | 17 / 25 | 20 / 30 | 26 / 40 | 1+NC (+3 per node/face) |
| Beer–Lambert + χ as (G, M) (+ I_f; A, θ_A global) | 24 / 36 | 27 / 41 | 33 / 51 | 3+NC (+3 per node/face, +2 global) |

**Formula.** A cell row of the model reads s·(NC + 1 + 2χ) cell unknowns (s = 3 for c0, 5 for c_active; χ = 1 adds the two inventories G and M), s + 1 node positions, two node velocities, two face intensities under Beer–Lambert (bl = 1) and the cross-section column under χ, so its degree, the clique lower bound of the coloring, is

    K(s) = s·(NC + 1 + 2χ) + (s + 1) + 2 + 2·bl + χ,

i.e. K_c0 = 3·NC + 9 + 2·bl + 7·χ and K_c_active = 5·NC + 13 + 2·bl + 11·χ. The greedy count equals this formula in 80 of 80 lifted rows and equals its lower bound in 80 of them. The previous model (13 to 28 colors, formula 3(NC+1) + 4 (+2 bl, +4 χ)) omitted the two outer node positions z_{i−1} and z_{i+2} that a cell row reads through the neighbouring cell centres and the neighbours' ξ = m/V, carried χ as one unknown instead of the inventory pair (G, M) and had no stencil for the relative-flux reconstruction; the corrected c0 counts are two higher without χ and five higher with it. The χ count lumps G and M into every cell row as the plan's layout does; a tighter model in which the T and ξ rows do not read the neighbours' G, M (with c = 0 they enter those rows only through A) would give 3·NC + 12 for Beer–Lambert + χ, six fewer, and is not measured here.

## Table 2: dual-pass overhead and projected colored-AD Jacobian cost (hot state, medians, μs)

| deck | n | t_residual | t_struct_exact | t_struct_reduced | per-color overhead C=4 / 8 / 16 / 32 | K c0 / c_active | projected c0 (C) | c0 / exact | c0 / reduced | projected c_active (C) | c_active / exact | c_active / reduced | measured colored AD on current pattern |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| eulerian | 100 | 12.9 | 89.8 | 89.2 | 0.41 / 0.36 / 0.33 / 0.34 | 15 / 23 | 80.1 (C=16) | 0.89 | 0.90 | 111.9 (C=24) | 1.25 | 1.25 | 70.5 μs, 9 colors, max|ΔJ| = 2.8e-14 (max|J| = 209) |
| eulerian | 400 | 38.0 | 278.0 | 276.7 | 0.70 / 0.58 / 0.50 / 0.53 | 15 / 23 | 344.0 (C=16) | 1.24 | 1.24 | 492.0 (C=24) | 1.77 | 1.78 | 288.6 μs, 9 colors, max|ΔJ| = 6.8e-13 (max|J| = 2241) |
| ale | 100 | 11.7 | 1438.3 | 83.3 | 0.65 / 0.53 / 0.44 / 0.45 | 15 / 23 | 92.9 (C=16) | 0.06 | 1.12 | 130.2 (C=24) | 0.09 | 1.56 | – |
| ale | 400 | 43.1 | 25584.6 | 326.5 | 0.73 / 0.61 / 0.51 / 0.56 | 15 / 23 | 396.6 (C=16) | 0.02 | 1.21 | 564.9 (C=24) | 0.02 | 1.73 | – |
| ale_darcy | 100 | 19.4 | 2083.9 | 157.7 | 0.59 / 0.46 / 0.37 / 0.36 | 17 / 25 | 184.2 (C=24) | 0.09 | 1.17 | 244.3 (C=32) | 0.12 | 1.55 | – |
| ale_darcy | 400 | 72.6 | 35343.6 | 633.9 | 0.65 / 0.52 / 0.42 / 0.46 | 17 / 25 | 803.5 (C=24) | 0.02 | 1.27 | 1125.0 (C=16) | 0.03 | 1.77 | – |

## Interpretation

1. **Color counts.** The lifted pattern's color count is independent of n. c0 / c_active: 15 / 23 for the NC = 2 surface-absorption decks (eulerian, ale) at n = 25, 100 and 400, 17 / 25 with Beer–Lambert (ale_darcy), 24 / 36 with Beer–Lambert and χ as (G, M) (ale_darcy_withchi), and 20 / 30 for the NC = 3 V-BO-01 deck; every count equals its max-row-degree lower bound and the formula above. The current ALE `:full` pattern colors in exactly the state length (101, 401, 1601 for ale at n = 25, 100, 400; 911 for V-BO-01): the top node's velocity row reads every column, so column coloring degenerates to one color per column, i.e. O(n). The current Eulerian `:full` and the ALE `:ale_reduced` patterns color in 9 (NC = 2) and 12 (NC = 3), fewer than the lifted pattern, because they carry no node, intensity or velocity columns (and the reduced pattern drops the geometry couplings that the lifted rows make local).
2. **Containment against the executable residual** (`lifted_containment.csv`; decks eulerian, ale, ale_darcy, ale_darcy_withchi at n = 6 and 8, dense ForwardDiff of the current residual at a perturbed reacting state on a non-uniform mesh, threshold 1e-12 × row max). With the production residual (no velocity override: its cell rows never read the node velocity, the prefix sum feeds only the z rows that the lifted balance rows replace) every T, ξ_j and χ row's measured cells and nodes lie inside the model's c0 sets: the measured reach is cells i-1..i+1 and nodes i-1..i+2, i.e. exactly the corrected c0 sets, with both outer nodes z_{i−1} and z_{i+2} present (the solid ξ row is cell-local, the gas ξ and T rows read the neighbours and the four nodes). Under Beer–Lambert the local sets were measured at I_radiation = 0 and the radiation path as the support of J(50 kW/m²) − J(0): it lies within the closure of the intensity chain (cells i..n, nodes i..n+1). Under WithChi the χ row reads every χ column through A(χ̄), which the model carries as the global A column; the T and ξ rows read no χ column in planar geometry (A cancels from the per-volume rows), fewer than the current `:full` pattern's dense χ coupling. The current Eulerian `:full` pattern's cell-to-cell width (i-1..i+1) and the measured Eulerian width are both inside the c0 stencil. With `Workspace.w_override = (z, t) -> 0` the z rows vanish, but the verification-only TVD temperature advection then runs with w_rel = 0 − w_physical, which is the prefix sum, so the T rows gain every cell and node below them, and under WithChi also everything above through the column-global θ_A; the reach above the cell from the limiter's own stencil is the c_active set (cells to i+2, nodes to i+3) on the NoChi decks: ale n=6: i-4..i+2 cells, i-4..i+3 nodes, upper reach within c_active = true; ale n=8: i-6..i+2 cells, i-6..i+3 nodes, upper reach within c_active = true; ale_darcy n=6: i-4..i+2 cells, i-4..i+3 nodes, upper reach within c_active = true; ale_darcy n=8: i-6..i+2 cells, i-6..i+3 nodes, upper reach within c_active = true; ale_darcy_withchi n=6: i-4..i+4 cells, i-4..i+5 nodes, upper reach within c_active = false; ale_darcy_withchi n=8: i-6..i+6 cells, i-6..i+7 nodes, upper reach within c_active = false. That run documents the hook, not the production or the lifted residual; the containment verdict is the no-override run. Nothing had to be added to the model after the measurement.
3. **Dual-pass overhead.** (t_dual(C) − t_residual)/(C·t_residual) is 0.32 to 0.73 residual-equivalents per partial across the decks and chunk widths (lowest around C = 16–24, highest at C = 4 and at n = 400 on the plain ALE deck); one dual pass of width C costs about 1 + 0.5·C residuals at C = 8 to 16.
4. **Projected colored-AD Jacobian cost** (ceil(K/C)·t_dual(C) at the best C, K the lifted color count, t_dual the CURRENT residual's dual-pass time; the lifted residual, which has more rows per cell, was not timed, so these are projections, not measurements). Relative to the structured assemblies, reported per size because the ALE exact assembly's O(n²) geometry block grows with n:

   | deck | n | c0: / exact | c0: / reduced | c_active: / exact | c_active: / reduced |
   |---|---|---|---|---|---|
   | eulerian | 100 | 0.89 | 0.90 | 1.25 | 1.25 |
   | eulerian | 400 | 1.24 | 1.24 | 1.77 | 1.78 |
   | ale | 100 | 0.06 | 1.12 | 0.09 | 1.56 |
   | ale | 400 | 0.02 | 1.21 | 0.02 | 1.73 |
   | ale_darcy | 100 | 0.09 | 1.17 | 0.12 | 1.55 |
   | ale_darcy | 400 | 0.02 | 1.27 | 0.03 | 1.77 |

   Eulerian decks (exact and reduced are the same pattern and cost): c0 0.89× (n = 100) and 1.24× (n = 400) the structured assembly, c_active 1.25× and 1.77×. ALE decks: c0 0.06–0.09× (n = 100) and 0.02× (n = 400) the structured exact assembly, 1.12–1.17× and 1.21–1.27× the reduced assembly; c_active 0.09–0.12× and 0.02–0.03× the exact, 1.55–1.56× and 1.73–1.77× the reduced. The exact ALE assembly costs 25.6–35.3 ms at n = 400 (1.4–2.1 ms at n = 100), so the ratio to it is size-dependent and must not be quoted as one number. The real colored Jacobian on the current Eulerian `:full` pattern (9 colors, one dual pass) runs at 0.79× and 1.04× the structured exact assembly at n = 100 and 400 and matches it to round-off.
5. **Against the 4-color toy ladder** (`results/restart_comparison.csv`, lifted form, 4 colors). Over the five decks (NC = 2 and 3) the production K is 15–24 (c0) and 23–36 (c_active), so K/4 = 3.75 to 6.00 (c0) and 5.75 to 9.00 (c_active); including the NC = 5 model ladder, K reaches 33 (c0) and 51 (c_active). In width-16 dual passes, ceil(K/16) = 1–2 (c0) and 2–3 (c_active) passes against the ladder's one, and the corresponding cost ratio ceil(K/16)·t_dual(16) / t_dual(4) (a production Jacobian at C = 16 over a 4-color Jacobian in one width-4 pass, measured dual times of the timed decks) is 2.2 to 4.3 (c0) and 4.1 to 4.8 (c_active).

Caveats: timings are single-threaded medians of 7 samples on a hot synthetic state at t = 1 s; the dual pass is through the current residual and the projection multiplies its cost by the lifted color count, as the plan prescribes; the lifted residual itself (more rows per cell: z, w, I, G and M) was not timed, so its t_residual, dual cost and colored-AD cost are not measured here; the global-scalar rows of the χ decks are assumed hand-assembled outside the coloring; the c_active stencil is a model of the planned reconstruction and has no executable counterpart in the current residual beyond the verification-only TVD kernel; the V-BO-01 deck is sized by its CI cell count (n = 182) and was not timed.
