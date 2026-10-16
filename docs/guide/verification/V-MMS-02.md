# V-MMS-02 — ALE free-stream preservation under prescribed mesh motion

**Tags:** `mms`, `ale`, `gcl`, `free_stream`

**References:**

- Thomas & Lombard, AIAA J. 17 (1979) — geometric conservation law
- docs/guide/technical_reference/18_verification-validation.md §18.2 (MMS driver, ALE cases)

## Problem statement

The sharpest possible discrete-GCL assert: spatially uniform T*, ξ* with no
reactions and zero-flux manufactured BCs, while the mesh is compressed to
70 % of its initial thickness through the `Workspace.w_override` hook with
the uniform-stretch velocity

  w*(z, t) = (ṡ/s)·z,   s(t) = 1 − 0.15·(1 − cos(2πt/20)),   z_bottom pinned

(node trajectories z_i(t) = z_i(0)·s(t) exactly). A uniform state must stay
uniform to machine precision on the moving mesh: the ALE advection of a
constant field must vanish discretely (limiter slope ≡ 0), the dilation
correction is zero (θ reaction-driven, no reactions), and the conductive/
species fluxes are zero. Any drift is a moving-mesh transport defect that
smooth-field order studies (V-MMS-03) can hide inside their truncation
error.

The node-motion ODE ż_i = w*(z_i,t) is integrated by the same stiff
integrator, so the top-cell-center trajectory doubles as an exactness check
on the prescribed-motion path (temporal tolerance only).

Jacobian: colored/dense FD (`jacobian = :fd`) — the structured geometry
Jacobian assumes reaction-driven mesh velocity and does not model the
override (documented on `Workspace.w_override`).

### Quantities of interest (n = 20)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| free-stream T drift (all cells, all snapshots) | 5.684e-14 | 0 | 5.684e-14 | 1e-10 | yes | uniform field must survive 30% mesh compression |
| free-stream ξ_solid relative drift | 1.895e-16 | 0 | 1.895e-16 | 1e-12 | yes | ALE advection + dilation must vanish for a uniform state |
| free-stream ξ_gas relative drift | 2.22e-15 | 0 | 2.22e-15 | 1e-12 | yes | ALE advection + dilation must vanish for a uniform state |
| top cell-center trajectory vs z₀·s(t_end) | 0.006825 | 0.006825 | 4.039e-07 | 1.5e-06 | yes | node ODE ż = (ṡ/s)z integrated at solver tolerance |

![V-MMS-02 thickness_history](assets/V-MMS-02/V-MMS-02_thickness_history.png)

![V-MMS-02 thickness_history_error](assets/V-MMS-02/V-MMS-02_thickness_history_error.png)

### Solver configuration

| setting | value |
| --- | --- |
| integrator | KenCarp4 (default) |
| abstol | 1.0e-10 |
| reltol | 1.0e-8 |

