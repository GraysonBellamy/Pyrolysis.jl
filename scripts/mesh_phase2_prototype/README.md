# Phase 2 prototype and baselines

Planning evidence for `docs/Mesh_Framework_Phase2_Plan.md` (Appendix A). Nothing
here is package code, and nothing here is run by `Pkg.test()`.

`common.jl` prototypes the correction `c = w − v_s` of the mesh-framework design
on the finite-volume core by evaluating replacement methods into the package's
modules: the per-node correction with its activation and monitor (evaluated once
per residual evaluation, by the relative-flux kernel, as the plan's `ds.c`; a
face whose relative velocity is exactly zero is skipped), the sweep term
in the balance rows and the consistent sweep (`Closures`), the relative-flux
kernel at the relative velocity `−c` (`Kernels`), the Jacobian pattern widened
for the velocity rows (`ColoredJacobian`). It also builds a `CoreProblem` with
the relative-flux kernel installed on a reacting material, which
`Integration.setup` refuses. Cells deactivate by the package's block-and-band
rule, the band given to `build_vbo` as a multiple `κ` of the threshold. The
tables of the plan's Appendix A were taken on `bb2860b`, before the package had
that rule, with the prototype's own version of it, which merged a run cell by
cell. The overrides persist in the Julia session: run each experiment in its
own process, and never in a session that runs the package's tests.

```sh
julia --project=test --startup-file=no --threads=1 scripts/mesh_phase2_prototype/run.jl <experiment>
```

with `<experiment>` one of `check` (containment of the widened pattern and the
colored Jacobian against the dense one), `sweep` (the response time `τ` on
V-BO-01 with the plain block rule), `band` (the near-threshold band and `τ` at
30, 60 and 120 PMMA cells), `repeat` (three solves of one problem), `ladder`
(grids, thresholds, the monitor's mass term, the activation, the energy form),
`order` (the band runs in another order), `continuous` (event location by a
continuous callback on V-BO-02, which also overrides the deactivation callback)
and `scaling` (the response time scaled by the square of the cell count across
the grid ladder).

`baselines.jl` holds the measurements that override nothing and may share a
session:

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

The plan's tables were taken on a shared machine; timings move by up to 40 %
between runs, so compare the two cores back to back in one process and repeat.
