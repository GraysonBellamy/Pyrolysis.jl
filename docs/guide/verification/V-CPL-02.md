# V-CPL-02 — Landau ablation: steady surface regression under ALE

**Tags:** `ale`, `moving_boundary`, `kinetics`, `surface_regression`, `coupled`

**References:**

- H.G. Landau, Q. Appl. Math. 8:81 (1950)

## Problem statement

H.G. Landau, Q. Appl. Math. 8:81 (1950). An analytical ALE check.
Semi-infinite solid, constant net
surface flux q″, material removed at ablation temperature T_p with heat of
ablation ΔH_p. Exact steady traveling wave (exact/landau.jl):

  v_ss = q″/(ρ[c(T_p − T₀) + ΔH_p]),  T(ζ) = T₀ + (T_p − T₀)e^(−v_ss ζ/α)

Realization: single reaction solid → gas, the gas outside the strain
(volume loss at the solid density ⇒ surface regression) under ALE with cell
deactivation. The kinetics-controlled front approximates the fixed-T_p
ablation limit; (A, E) are chosen on the leading-order Laplace locus

  A = v² E (T_p − T₀) e^(E/RT_p) / (α R T_p²)

so the designed regression rate is v_target at ablation temperature ≈ T_p.
At finite E the front self-adjusts: the primary QoI therefore compares the
measured regression rate v against the Landau balance evaluated at the
MEASURED surface temperature — an identity that must hold for the steady
wave regardless of kinetics — and the design locus is only a loose check.

With c_gas = c_solid the reaction enthalpy is constant. Integrating the
steady first law from stationary virgin solid at T₀ to the escaping gas
gives q″ = ρ v [c(T_s − T₀) + ΔH], independently of the front's kinetics.
This identity requires a stationary cold far field and negligible changes
in stored energy and gas holdup. The reaction gate leaves the virgin solid
inert. The broader-front companion uses a longer domain and later window;
both cases check the back temperature and report mesh refinement.

The half-cell surface-temperature closure and finite front resolution
contribute spatial error. The Laplace design locus remains a loose check;
it is not the analytical reference for the measured regression rate.
Prescribed net heat flux has no blowing correction. The depletion limiter
is active as front cells approach deactivation.

## V-CPL-02 — Landau ablation: steady surface regression under ALE

### Quantities of interest (n = 240)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| steady energy balance: v vs Landau v_ss(measured T_s) | 0.001954 | 0.001954 | 0.0002337 | 0.0021 | yes | steady first law with stationary virgin solid and c_gas = c_solid |
| regression consistency: thickness slope vs gas generation | 0.001954 | 0.001954 | 0.0002491 | 0.0021 | yes | gas outside the strain ⇒ dL/dt = −ṁ‴_gas/ρ exactly |
| mass ledger closure over plateau window | 0.778 | 0.7818 | 0.00484 | 0.016 | yes | total-mass drop vs ∫ ṁ″ dt; residual = holdup drift |
| in-depth profile decay rate (surface-attached frame) | 0.001969 | 0.001954 | 0.007773 | 0.027 | yes | ln(T−T₀) slope over ζ/δ ∈ [0.75, 1.7], 12 cells |
| design-locus regression rate (loose) | 0.001954 | 0.002 | 0.02306 | 0.1 | yes | leading-order Laplace locus only pins T_p to O(RT_p²/E); informational, tightens as E ↑ (see V-CPL-02b) |
| back-face rise (semi-infinite validity) | 0.0003584 | 0 | 0.0003584 | 0.01 | yes | guards the semi-infinite assumption |

![V-CPL-02 mlr_history](assets/V-CPL-02/V-CPL-02_mlr_history.png)

![V-CPL-02 mlr_history_error](assets/V-CPL-02/V-CPL-02_mlr_history_error.png)

![V-CPL-02 thickness_history](assets/V-CPL-02/V-CPL-02_thickness_history.png)

![V-CPL-02 thickness_history_error](assets/V-CPL-02/V-CPL-02_thickness_history_error.png)

![V-CPL-02 profile_surface_frame](assets/V-CPL-02/V-CPL-02_profile_surface_frame.png)

![V-CPL-02 profile_surface_frame_error](assets/V-CPL-02/V-CPL-02_profile_surface_frame_error.png)

![V-CPL-02 convergence](assets/V-CPL-02/V-CPL-02_convergence.png)

![V-CPL-02 param_E](assets/V-CPL-02/V-CPL-02_param_E.png)

### Convergence

| n_cells | h | wall (s) | balance | consistency |
| --- | --- | --- | --- | --- |
| 60 | 0.01667 | 0.1804 | 0.0002449 | 0.0001098 |
| 100 | 0.01 | 0.4438 | 0.0001066 | 0.0005932 |
| 160 | 0.00625 | 1.069 | 0.0002021 | 0.0006322 |
| 240 | 0.004167 | 2.371 | 0.0002337 | 0.0002491 |

Wall time: median of up to five warm runs of problem construction and solve on one pinned core, after a discarded compiling run.

Observed order -0.08591 (balance).

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| use_ale | true |
| min_thickness | 0.0001 |
| deactivation | true |
| abstol | 1.0e-7 |
| reltol | 1.0e-5 |

## V-CPL-02b — Landau ablation: kinetics-limit insensitivity variant (E/2)

### Quantities of interest (n = 320)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| steady energy balance: v vs Landau v_ss(measured T_s) | 0.001914 | 0.001913 | 0.000466 | 0.002 | yes | identity must close at E = 50 kJ/mol as at 100 (locus check) |
| regression consistency: thickness slope vs gas generation | 0.001914 | 0.001914 | 3.559e-05 | 0.0017 | yes | gas outside the strain ⇒ dL/dt = −ṁ‴_gas/ρ exactly |
| stationary far-field temperature | 3.481e-07 | 0 | 3.481e-07 | 0.0001 | yes | cold reaction gate and domain length qualify the steady reference |

![V-CPL-02b mlr_history](assets/V-CPL-02b/V-CPL-02b_mlr_history.png)

![V-CPL-02b mlr_history_error](assets/V-CPL-02b/V-CPL-02b_mlr_history_error.png)

![V-CPL-02b convergence](assets/V-CPL-02b/V-CPL-02b_convergence.png)

### Convergence

| n_cells | h | wall (s) | balance | consistency |
| --- | --- | --- | --- | --- |
| 120 | 0.008333 | 0.762 | 0.0004264 | 6.497e-05 |
| 200 | 0.005 | 1.999 | 0.000453 | 7.265e-05 |
| 320 | 0.003125 | 5.208 | 0.000466 | 3.559e-05 |

Wall time: median of up to five warm runs of problem construction and solve on one pinned core, after a discarded compiling run.

Observed order -0.09103 (balance).

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| use_ale | true |
| min_thickness | 0.0001 |
| deactivation | true |
| abstol | 1.0e-7 |
| reltol | 1.0e-5 |
