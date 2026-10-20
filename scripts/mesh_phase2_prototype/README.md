# Phase 2 baselines

Planning evidence for `docs/Mesh_Framework_Phase2_Plan.md` (Appendix A). Nothing
here is package code, and nothing here is run by `Pkg.test()`.

`baselines.jl` holds the measurements of the plan's appendix that override no
package method:

```sh
julia --project=test --startup-file=no --threads=1 scripts/mesh_phase2_prototype/baselines.jl <name>
```

with `<name>` one of `vbo_old`, `vbo_new` (the five burnout cases through the
suite's routing), `vbo_ladder_old`, `vbo01_qois_old` (the old core at 30, 60
and 120 PMMA cells), `vbo02_restart` (the restart policy after an event),
`df_cores` (the Douglas fir deck of `examples/07_douglas_fir_cone.jl` on both
cores), `integrators` (Rosenbrock and SDIRK candidates on the audit decks) and
`tolweights` (the tolerance weights of the algebraic and inventory rows).
`cases.jl` is the case runner `vbo_old` and `vbo_new` use.

The correction these baselines bracket is the package's
`EquidistributionCorrection`; the gate scripts of `test/newcore/gates/` measure
it (`burnout_ladder.jl`, `burnout_timing.jl`).

The plan's tables were taken on a shared machine; timings move by up to 40 %
between runs, so compare the two cores back to back in one process and repeat.
