# The iteration matrix on a moving mesh

Produced by `test/newcore/gates/w_formation.jl`. Julia 1.13.0, 1 thread, 22 logical CPUs; load average [21.5, 21.5, 22.1] at the start and [22.0, 21.8, 22.1] at the end.

The new core on `ale` at 100 cells (solid λ 1e-6, reltol 1e-7): 90 accepted steps; `J` stores 4780 to 4780 entries, `W` 2414 to 3623 (56 distinct counts); 90 distinct KLU factorizations, each a fresh symbolic analysis. With the structure-keeping formation `W` stores 4780 to 4780 entries and the solve makes 1 KLU factorization, refactored in place at every other step.

Warm minimum of three [s]: the new core at its matched tolerance, the old core at its default, with OrdinaryDiffEq's `jacobian2W!` (broadcast) and with a structure-keeping one. ΔT: the largest temperature difference the substitution makes to the new core's solution [K].

| deck | new, broadcast | new, kept | old, broadcast | old, kept | ratio, broadcast | ratio, kept | new kept / old broadcast | ΔT [K] |
|---|---|---|---|---|---|---|---|---|
| eulerian 100, solid λ 1.0e-6 | 0.0581 | 0.0586 | 0.04 | 0.0409 | 1.45 | 1.43 | 1.47 | 0.0 |
| ale 25, solid λ 1.0e-6 | 0.0276 | 0.0198 | 0.014 | 0.0136 | 1.97 | 1.45 | 1.41 | 5.5e-11 |
| ale 100, solid λ 1.0e-6 | 0.114 | 0.0842 | 0.0448 | 0.0437 | 2.55 | 1.93 | 1.88 | 5.88e-11 |
| ale_darcy 25 | 0.0558 | 0.0347 | 0.0337 | 0.0333 | 1.66 | 1.04 | 1.03 | 7.12e-6 |
| ale_darcy 100 | 0.205 | 0.181 | 0.134 | 0.137 | 1.53 | 1.32 | 1.36 | 5.22e-11 |
