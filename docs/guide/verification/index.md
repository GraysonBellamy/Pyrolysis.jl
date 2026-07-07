# Verification suite report

Exact-solution verification of Pyrolysis.jl: every case solves a
problem with a closed-form (or manufactured) reference solution and
reports error norms, pinned tolerances with provenance, and observed
convergence orders. Methodology, case selection, and the dual-mode
(CI assert / report) architecture are described in
`docs/Verification_Suite_Plan_2026-07-06.md`; the assert-mode twin of
this report runs in `Pkg.test()` from the same case definitions
(`test/verification/`).

Figures follow the FDS Verification Guide convention: lines are the
exact/reference solution, markers are the solver. Tolerances are
pinned at ~3× the observed error at pin time (the provenance column).

Where a case has a published analog in another solver, it carries a
"Comparison with other codes" section: FDS 6.11.0 and Gpyro 0.8200
were run locally on the same problems (decks and outputs committed
with provenance in `test/verification/reference/`), and their errors
against the same exact solutions are overlaid — muted gray, log
scale — together with a wall-time comparison at every resolution
rung.

Generated 2026-07-07 by
`scripts/verification/run_reports.jl`. Regenerate with:
`julia --project=scripts/verification scripts/verification/run_reports.jl`

| case | title | tags | QoIs within tol | worst error/tol | observed order |
| --- | --- | --- | --- | --- | --- |
| [V-CPL-01](V-CPL-01.md) | Neumann Stefan problem via apparent heat capacity (Gpyro params) | moving_front, latent_heat, state_dependent_properties, coupled | 9/9 | 0.3349 | — |
| [V-CPL-01b](V-CPL-01.md) | Neumann Stefan problem via gated reaction (±50 K rescaled) | moving_front, latent_heat, kinetics, temperature_gate | 9/9 | 0.3505 | — |
| [V-CPL-02](V-CPL-02.md) | Landau ablation: steady surface regression under ALE | ale, moving_boundary, kinetics, surface_regression, coupled | 6/6 | 0.6247 | 2.493 |
| [V-CPL-03](V-CPL-03.md) | shrink/swell exact endpoints (6 permutations + WithChi law) | ale, volume_change, swelling, chi, mass_conservation | 18/18 | 0.3269 | — |
| [V-CPL-04](V-CPL-04.md) | Frank-Kamenetskii slab: subcritical steady state + criticality bracket | reaction_conduction, bifurcation, exothermic, coupled_jacobian | 6/6 | 0.3826 | 2.474 |
| [V-CPL-05](V-CPL-05.md) | TGA-to-gasification consistency ladder (1D → 0D thin limit) | kinetics, consistency_ladder, thin_limit, bc_coupling | 4/4 | 0.334 | — |
| [V-GT-01](V-GT-01.md) | Gaussian gas-pulse diffusion (plane source) | gas_transport, fickian, mass_conservation | 7/7 | 0.2298 | 2.012 |
| [V-GT-02](V-GT-02.md) | steady compressible Darcy outflow (P² profile) | gas_transport, darcy, pressure_closure, steady | 4/4 | 0.1471 | 0.8798 |
| [V-HC-01](V-HC-01.md) | semi-infinite solid, convective surface BC (erfc) | heat_conduction, convective_bc, surface_newton | 9/9 | 0.01337 | 1.993 |
| [V-HC-02a](V-HC-02.md) | finite slab, convective BC, Bi = 0.1 | heat_conduction, convective_bc, biot_sweep | 8/8 | 0.003794 | 2 |
| [V-HC-02b](V-HC-02.md) | finite slab, convective BC, Bi = 1.0 | heat_conduction, convective_bc, biot_sweep | 8/8 | 0.01816 | 2 |
| [V-HC-02c](V-HC-02.md) | finite slab, convective BC, Bi = 10.0 | heat_conduction, convective_bc, biot_sweep | 8/8 | 0.02029 | 2.001 |
| [V-HC-02d](V-HC-02.md) | finite slab, convective BC, Bi = 100.0 | heat_conduction, convective_bc, biot_sweep | 8/8 | 0.02125 | 2 |
| [V-HC-03](V-HC-03.md) | finite slab, constant flux, adiabatic back | heat_conduction, heat_flux_bc, energy_balance | 5/5 | 0.02041 | 2 |
| [V-HC-04a](V-HC-04.md) | two-layer slab, k-ratio 10.0:1 | heat_conduction, laminate, harmonic_mean, steady | 4/4 | 0.003336 | — |
| [V-HC-04b](V-HC-04.md) | two-layer slab, k-ratio 100.0:1 | heat_conduction, laminate, harmonic_mean, steady | 4/4 | 0.04082 | — |
| [V-HC-04c](V-HC-04.md) | two-layer slab, k-ratio 10.0:1 (off-face interface) | heat_conduction, laminate, harmonic_mean, steady | 5/5 | 0.7869 | — |
| [V-HC-05](V-HC-05.md) | steady Kirchhoff profile with k(T) = k₀(1+βT) | heat_conduction, temperature_dependent_k, kirchhoff, steady | 2/3 | 21.04 | 2 |
| [V-MMS-01](V-MMS-01.md) | manufactured solution, coupled kinetics + T-props + gas transport | mms, coupled, kinetics, gas_transport, order_verification | 5/5 | 0.005304 | 1.992 |
| [V-MMS-02](V-MMS-02.md) | ALE free-stream preservation under prescribed mesh motion | mms, ale, gcl, free_stream | 4/4 | 0.2693 | — |
| [V-MMS-03](V-MMS-03.md) | ALE manufactured solution, prescribed smooth mesh motion | mms, ale, gcl, order_verification | 3/3 | 0.1033 | 2.156 |
| [V-MMS-04](V-MMS-04.md) | Induced-motion ALE manufactured solution (reaction-driven mesh) | mms, ale, kinetics, order_verification, approx_jacobian | 6/6 | 0.3324 | 2.014 |
| [V-RAD-01](V-RAD-01.md) | in-depth absorption + conduction (sub-surface peak) | radiation, beer_lambert, convective_bc, coupled | 11/11 | 0.1733 | 1.995 |
| [V-RX-01](V-RX-01.md) | two-step series chain A→B→C, isothermal | kinetics, multi_step, isothermal | 5/5 | 0.2777 | — |
| [V-RX-02a](V-RX-02.md) | zero-order decay (exact burnout time) | kinetics, reaction_order, isothermal | 2/2 | 0.002842 | — |
| [V-RX-02b](V-RX-02.md) | second-order decay (ThermaKin C+C analog) | kinetics, reaction_order, isothermal | 2/2 | 0.345 | — |
| [V-RX-02c](V-RX-02.md) | contracting-sphere R3 decay (FDS NaHCO₃, 450 K) | kinetics, reaction_order, isothermal | 2/2 | 0.2772 | — |
| [V-RX-02d](V-RX-02.md) | zero-order decay, depletion limiter ON (guard) | kinetics, depletion_limiter, isothermal | 2/2 | 0.001152 | — |
| [V-RX-03](V-RX-03.md) | competing parallel reactions with yield splits | kinetics, parallel_reactions, stoichiometry, energy | 6/6 | 0.2767 | — |
| [V-RX-04a](V-RX-04.md) | thermally-thin radiative cooling (implicit exact) | radiation, radiative_bc, surface_newton, thermally_thin | 2/2 | 0.3244 | — |
| [V-RX-04b](V-RX-04.md) | steady radiative–convective surface equilibrium | radiation, radiative_bc, convective_bc, surface_newton, steady | 3/3 | 0.1074 | — |

Cases without report-mode output (no `report` metadata yet): V-CPL-02b, V-CPL-06a, V-CPL-06b.
