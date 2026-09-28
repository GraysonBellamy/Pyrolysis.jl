# V-MMS-03 — ALE manufactured solution, prescribed smooth mesh motion

**Tags:** `mms`, `ale`, `gcl`, `order_verification`

**References:**

- docs/Verification_Suite_Plan_2026-07-06.md Tier 3 (ALE-MMS deliverable)
- Thomas & Lombard, AIAA J. 17 (1979) — geometric conservation law

## Problem statement

The moving-mesh order-verification companion to V-MMS-01 (and the "hardest
and most valuable single test" of the suite plan's Tier 3): non-trivial
manufactured fields solved on a mesh compressed to 70 % of its thickness by
the prescribed uniform-stretch velocity

  w*(z, t) = (ṡ/s)·z,   s(t) = 1 − 0.15·(1 − cos(2πt/20)),  z_i(t) = z_i(0)·s(t)

through the `Workspace.w_override` hook. On a grid point moving with w, the
grid time derivative is ∂tT* + w·∂zT*; the discrete ALE advection term
supplies the w·∂zT* part, so the same Eulerian manufactured sources —
evaluated at the LIVE cell centers read from the z-block — make T*, ξ* the
exact solution on the moving mesh (see exact/mms.jl header, "ALE"). The
error norms therefore measure exactly the moving-mesh transport truncation:
ALE advection (van-Leer-limited upwind reconstruction) + lab-frame
operators on a nonuniform-in-time mesh.

No reactions: mesh motion is prescribed, not volume-change-driven, so the
case isolates GCL/advection from kinetics (V-MMS-01 covers those). The
solid field is spatially UNIFORM but time-varying: condensed phases are
mesh-bound in this ALE formulation (v_solid = w ⇒ no advection operator,
see compute_ale_species_advection!), so a z-varying solid target under
prescribed non-material motion would have no operator to balance its
w·∂zξ* term (constraint documented in exact/mms.jl). The time variation
keeps the solid source-injection path active; T and the gas field carry
the full space–time structure and exercise the moving-mesh advection.

Expected order: the limiter's slope reconstruction is formally 2nd order on
smooth data but clips at extrema (the sin peak crosses the domain), so the
observed order sits between 1.5 and 2 — recorded at pin time; a drop below
min_order flags a moving-mesh transport regression.

Jacobian: dense FD (`jacobian = :fd`) — the structured geometry Jacobian
assumes reaction-driven mesh velocity (see `Workspace.w_override` note).

### Quantities of interest (n = 128)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| L∞ T error on moving mesh, t=10.0 s | 0.003857 | — | 0.003857 | 0.25 | yes | manufactured T* on 30%-compressed mesh |
| relative L∞ ξ_solid (mesh-bound) error, t=10.0 s | 6.204e-12 | — | 6.204e-12 | 1e-10 | yes | uniform solid: verifies the Lagrangian (no-advection) solid path + source injection |
| relative L∞ ξ_gas error, t=10.0 s | 0.0006197 | — | 0.0006197 | 0.006 | yes | gas is advected relative to the mesh: moving-mesh transport error channel |

![V-MMS-03 T_profiles](assets/V-MMS-03/V-MMS-03_T_profiles.png)

![V-MMS-03 T_profiles_error](assets/V-MMS-03/V-MMS-03_T_profiles_error.png)

![V-MMS-03 xi_profiles](assets/V-MMS-03/V-MMS-03_xi_profiles.png)

![V-MMS-03 xi_profiles_error](assets/V-MMS-03/V-MMS-03_xi_profiles_error.png)

![V-MMS-03 xi_gas_profile](assets/V-MMS-03/V-MMS-03_xi_gas_profile.png)

![V-MMS-03 xi_gas_profile_error](assets/V-MMS-03/V-MMS-03_xi_gas_profile_error.png)

![V-MMS-03 error_profiles](assets/V-MMS-03/V-MMS-03_error_profiles.png)

![V-MMS-03 thickness_history](assets/V-MMS-03/V-MMS-03_thickness_history.png)

![V-MMS-03 thickness_history_error](assets/V-MMS-03/V-MMS-03_thickness_history_error.png)

![V-MMS-03 convergence](assets/V-MMS-03/V-MMS-03_convergence.png)

### Convergence

| n_cells | h | wall (s) | L2_T | Linf_T |
| --- | --- | --- | --- | --- |
| 16 | 0.0625 | 0.1125 | 0.2321 | 0.3061 |
| 32 | 0.03125 | 0.6345 | 0.05006 | 0.06913 |
| 64 | 0.01562 | 1.222 | 0.01107 | 0.01579 |
| 128 | 0.007812 | 8.375 | 0.002632 | 0.003857 |

Observed order 2.156 (L2_T), expected 2.0.

### Solver configuration

| setting | value |
| --- | --- |
| integrator | KenCarp4 (default) |
| abstol | 1.0e-11 |
| reltol | 1.0e-9 |

