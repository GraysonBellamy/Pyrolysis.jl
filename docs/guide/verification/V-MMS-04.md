# V-MMS-04 — Induced-motion ALE manufactured solution (reaction-driven mesh)

**Tags:** `mms`, `ale`, `kinetics`, `order_verification`, `approx_jacobian`

**References:**

- docs/Verification_Suite_Plan_2026-07-06.md Tier 3 + cross-cutting requirement 3

## Problem statement

The Tier-3 completion case: V-MMS-02/03 verify GCL + moving-mesh transport
under PRESCRIBED motion (`w_override`) and deliberately bypass the θ
prefix-sum mesh velocity, the dilation correction −ξθ with θ ≠ 0, and the
structured geometry Jacobian. Here the kinetics are manufactured so the
rates are analytic on the manufactured fields and the motion is genuinely
REACTION-DRIVEN — no `w_override` — closing the kinetics→θ→w→geometry
chain and making both `DirectSolve` AND `ApproxSolve` valid backends. The
fixed-resolution ApproxSolve solve whose exact-error QoI sits next to
DirectSolve's is cross-cutting requirement 3 of the suite plan (quantify
ApproxSolve against a known-exact answer).

## Manufactured design (see exact/mms.jl, "Induced mesh motion")

Constant intrinsic densities, γ_solid = 1 / γ_gas = 0 (canonical), single
first-order reaction s → g with the depletion limiter numerically inactive
over the manufactured range (tanh(ξ_s/1.0) ≡ 1.0 in Float64 for
ξ_s ≳ 20 kg/m³; here ξ_s ≥ 500). T* is spatially UNIFORM (time-varying):
V-MMS-01/03 already order-verify the spatial T operators — this case's
target is the kinetics→geometry chain. The solid is z-VARYING,
ξ_s* = a(t) + b(t)·sin(πz/L): representable here (and only here) because
the mesh moves at the solid material velocity, so the material-frame source
S_s = ∂tξ* + w*·∂zξ* − ω̇* + ξ*·θ* has the moving-point chain rule built
in (`mms_species_source_material` via the `MMSMotion` companion).

  θ*(z,t) = −rate*/ρ_s = −k(T*(t))·ξ_s*(z,t)/ρ_s   (volume_change_rate)
  w*(z,t) = ∫₀ᶻ θ* dζ = −(k(t)/ρ_s)·[a(t)·z + b(t)·(L/π)(1 − cos(πz/L))]

The case builds w* with `mms_mesh_velocity` (32-node Gauss–Legendre on the
package kinetics — machine precision, requirement 4); the closed form above
is the independent hand-check in exact/selftests.jl. θ* < 0 throughout ⇒
monotone ~8 % compression over 20 s, every node inside the maximal domain
[0, L].

## Moving-boundary BC positions

Induced motion means the exact top-face trajectory solves the 1D node ODE
dz*/dt = w*(z*, t), z*(0) = L — not analytic here, so `build` precomputes it
with Vern9 at 1e-14 tolerances (requirement-4 compliant) and the BC closures
read the dense interpolant. The bottom node is pinned (w*(0,t) = 0).

## Jacobian caveat (do not claim exactness)

With sources evaluated at live cell centers, du += S(z(u), t) has a real
∂S/∂z chain the analytic `Structured(DirectSolve())` Jacobian does not
carry (driver design decision 1's ∂S/∂u ≡ 0 holds only on a fixed mesh).
Harmless for solution correctness and observed order (inexact Newton
matrix), and the DirectSolve-vs-ApproxSolve comparison stays clean because
both share the omission.

### Quantities of interest (n = 128)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| L∞ T error (DirectSolve), t=20.0 s | 1.079e-08 | — | 1.079e-08 | 1e-06 | yes | z-uniform T* — kinetics/source channel, not spatial T operators |
| relative L∞ ξ_solid (z-varying, material frame) error, t=20.0 s | 4.088e-08 | — | 4.088e-08 | 2e-06 | yes | solid rides the reaction-driven mesh: prefix-sum w + dilation −ξθ + material-frame source |
| relative L∞ ξ_gas error, t=20.0 s | 0.0001505 | — | 0.0001505 | 0.0065 | yes | gas advected relative to the reacting, compressing mesh |
| relative surface-position error, t=20.0 s | 1.008e-06 | 0.008963 | 1.008e-06 | 5e-05 | yes | \|z_top − z*_top\| / (L − z*_top); recession ≈ 1.04 mm |
| L∞ T error (ApproxSolve), t=20.0 s | 2.351e-07 | — | 2.351e-07 | 1e-06 | yes | requirement 3: ApproxSolve vs known-exact answer |
| relative L∞ ξ_solid error (ApproxSolve), t=20.0 s | 6.647e-07 | — | 6.647e-07 | 2e-06 | yes | requirement 3: ApproxSolve vs known-exact answer |

![V-MMS-04 T_profiles](assets/V-MMS-04/V-MMS-04_T_profiles.png)

![V-MMS-04 T_profiles_error](assets/V-MMS-04/V-MMS-04_T_profiles_error.png)

![V-MMS-04 xi_profiles](assets/V-MMS-04/V-MMS-04_xi_profiles.png)

![V-MMS-04 xi_profiles_error](assets/V-MMS-04/V-MMS-04_xi_profiles_error.png)

![V-MMS-04 xi_gas_profile](assets/V-MMS-04/V-MMS-04_xi_gas_profile.png)

![V-MMS-04 xi_gas_profile_error](assets/V-MMS-04/V-MMS-04_xi_gas_profile_error.png)

![V-MMS-04 error_profiles](assets/V-MMS-04/V-MMS-04_error_profiles.png)

![V-MMS-04 thickness_history](assets/V-MMS-04/V-MMS-04_thickness_history.png)

![V-MMS-04 thickness_history_error](assets/V-MMS-04/V-MMS-04_thickness_history_error.png)

![V-MMS-04 convergence](assets/V-MMS-04/V-MMS-04_convergence.png)

### Convergence

| n_cells | h | wall (s) | L2_xi_s | Linf_xi_g |
| --- | --- | --- | --- | --- |
| 16 | 0.0625 | 0.1284 | 0.0007645 | 0.01281 |
| 32 | 0.03125 | 0.2717 | 0.0001908 | 0.002076 |
| 64 | 0.01562 | 0.8112 | 4.739e-05 | 6.73e-05 |
| 128 | 0.007812 | 3.066 | 1.159e-05 | 0.0003011 |

Observed order 2.014 (L2_xi_s), expected 2.0.

### Solver configuration

| setting | value |
| --- | --- |
| integrator | KenCarp4 (default) |
| abstol | 1.0e-11 |
| reltol | 1.0e-9 |

