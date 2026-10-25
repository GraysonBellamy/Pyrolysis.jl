# V-GT-01 — Gaussian gas-pulse diffusion (plane source)

**Tags:** `gas_transport`, `lumped_transfer`, `mass_conservation`

**References:**

- J. Crank, The Mathematics of Diffusion (plane-source solution)
- Stoliarov & Lyon, FSS 9:1141 (2008), gas-transport verification

## Problem statement

Pure lumped wall transport of gas released inside the solid (Crank/ThermaKin
plane source). Instead of approximating a delta function, the initial
condition is a finite-width Gaussian — itself an exact solution — so the
comparison is IC-discretization-free:

  ξ(z,t) = M_tot/√(2πσ²(t)) · exp(−(z−z₀)²/(2σ²(t))),  σ²(t) = σ₀² + 2dt

At constant physical porosity φ and temperature, J=−(Λ/φ)∂zξ.
The Gaussian variance therefore grows with diffusivity d=Λ/φ. Adiabatic
faces, no reaction/source, and internally carried enthalpy preserve the
initially uniform temperature; this is checked in the QoIs.
The sealed image-source correction is below the spatial error on this domain.

QoIs: profile L∞ (relative to the evolving peak), total-mass conservation,
and the parameter-free variance growth law Var(t) = σ₀² + 2dt.

### Quantities of interest (n = 640)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| isothermal invariant | 5.684e-14 | — | 5.684e-14 | 1e-08 | yes | sealed, source-free uniform-temperature state |
| L∞ ξ-profile error / peak, t=2.0 s | 7.808e-05 | — | 7.808e-05 | 0.004 | yes | Gaussian stays Gaussian; n=160 uniform, σ₀ = 4Δz |
| L∞ ξ-profile error / peak, t=4.0 s | 4.812e-05 | — | 4.812e-05 | 0.004 | yes | Gaussian stays Gaussian; n=160 uniform, σ₀ = 4Δz |
| L∞ ξ-profile error / peak, t=8.0 s | 2.687e-05 | — | 2.687e-05 | 0.004 | yes | Gaussian stays Gaussian; n=160 uniform, σ₀ = 4Δz |
| gas-mass conservation (sealed) | 2.595e-15 | — | 2.595e-15 | 1e-13 | yes | antisymmetric face fluxes telescope exactly |
| variance growth, t=2.0 s | 5e-06 | 5e-06 | 2.563e-06 | 2e-05 | yes | Var(t) = σ₀² + 2dt; relative error reported |
| variance growth, t=4.0 s | 9e-06 | 9e-06 | 3.654e-06 | 2e-05 | yes | Var(t) = σ₀² + 2dt; relative error reported |
| variance growth, t=8.0 s | 1.7e-05 | 1.7e-05 | 5.196e-06 | 2e-05 | yes | Var(t) = σ₀² + 2dt; relative error reported |

![V-GT-01 xi_profiles](assets/V-GT-01/V-GT-01_xi_profiles.png)

![V-GT-01 xi_profiles_error](assets/V-GT-01/V-GT-01_xi_profiles_error.png)

![V-GT-01 variance_history](assets/V-GT-01/V-GT-01_variance_history.png)

![V-GT-01 variance_history_error](assets/V-GT-01/V-GT-01_variance_history_error.png)

![V-GT-01 convergence](assets/V-GT-01/V-GT-01_convergence.png)

### Convergence

| n_cells | h | wall (s) | L2 | Linf |
| --- | --- | --- | --- | --- |
| 40 | 0.025 | 0.01526 | 0.003984 | 0.01214 |
| 80 | 0.0125 | 0.05138 | 0.0009693 | 0.003074 |
| 160 | 0.00625 | 0.1978 | 0.0002408 | 0.0007707 |
| 320 | 0.003125 | 0.8064 | 6.02e-05 | 0.0001928 |
| 640 | 0.001563 | 2.677 | 1.54e-05 | 4.812e-05 |

Observed order 2.004 (L2), expected 2.0.

### Solver configuration

| setting | value |
| --- | --- |
| integrator | Rodas5P (default) |
| abstol | 1.0e-12 |
| reltol | 1.0e-10 |
