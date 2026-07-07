# V-CPL-03 — shrink/swell exact endpoints (6 permutations + WithChi law)

**Tags:** `ale`, `volume_change`, `swelling`, `chi`, `mass_conservation`

**References:**

- FDS Verification Guide `shrink_swell` (endpoint concept)
- Technical reference Ch. 10/11 (volumetric strain model)

## Problem statement

FDS verifies exact final thickness/density for ρ_product ≷ ρ_reactant
permutations. Under OUR volumetric-strain model the endpoint follows from
(technical reference Ch. 10/11; src/Physics/volume_change.jl):

  V̇/V = θ = Σⱼ (γⱼ/ρⱼ) ω̇ⱼ

with ω̇ⱼ the reaction mass rates [kg/m³/s]. For a single reaction
A → Σⱼ yⱼ Bⱼ (Σyⱼ = 1) at rate r: ω̇_A = −r, ω̇_Bⱼ = +yⱼ r, so per cell

  dV = (γ_A/ρ_A − Σⱼ yⱼ γⱼ/ρⱼ) dm_A        (r V dt = −dm_A)

— LINEAR in the converted mass, hence path-independent: the endpoint is
exact regardless of rate history, temperature, or spatial distribution.
Integrating from m_A = ξ_A0·V₀ to 0, with the standard initialization
ξ_A0 = ρ_A (bulk density) and γ_A = 1:

  V_f/V₀ = 1 − γ_A + Σⱼ yⱼ γⱼ ρ_A/ρⱼ  =  Σⱼ yⱼ γⱼ ρ_A/ρⱼ

and a pure solid product ends exactly at bulk density ρ_B/γ_B. All ρ here
are the components' BULK densities (the strain kernel uses bulk ρ, not
ρ_intrinsic — see the SolidComponent docstring note).

NOTE the plan's "γ > 1 intumescent branch" is NOT representable: γ is
validated to [0, 1] (universal volume-fraction semantics shared with
mixture density and conductivity mixing). Swelling is exercised through
ρ_B < ρ_A (P1) and the partial-contribution branch through γ_B ∈ (0, 1)
(P6). Documented as a model limitation.

Permutations (ρ_A = 1000, γ_A = 1, isothermal via E = 0, h = 0, everything
impermeable). Mass QoIs: EVERY permutation asserts total-mass conservation
to integrator tolerance — the gas mesh-motion transport is the conservative
flux form (tech-ref §11.2.3), so held gas in a collapsing impermeable
matrix compresses (ξ_g ↑ as V ↓) instead of being deleted with the
vanishing volume (a non-telescoping pointwise form would lose tens of
percent of total mass in P4 — that configuration is the regression
anchor). Gas-producing permutations additionally assert the SOLID ledger
(char endpoint = y_char·m_A0 exactly).

  P1 swell:            A → B(ρ=500)               V_f/V₀ = 2
  P2 shrink:           A → B(ρ=2000)              V_f/V₀ = 1/2
  P3 neutral:          A → B(ρ=1000)              V_f/V₀ = 1
  P4 char + gas:       A → 0.3 C(ρ=1000) + 0.7 G  V_f/V₀ = 0.3
  P5 mass loss, ΔV=0:  A → 0.5 C(ρ=500)  + 0.5 G  V_f/V₀ = 1
  P6 partial γ:        A → B(ρ=1000, γ=0.5)       V_f/V₀ = 1/2
  P7 WithChi:          P5 chemistry + lateral law A/A₀ = 1 − 0.5·χ̄
       χ̄_f = 0.5 (released dry-mass fraction) ⇒ A_f/A₀ = 0.75 exactly,
       L_f/L₀ = (V_f/V₀)/(A_f/A₀) = 4/3.

P1 (the swelling branch, exercising mesh GROWTH) runs through the standard
harness path; the remaining permutations are solved inside `qois` (each is
a ~1 s, 8-cell isothermal solve) so the whole family stays one case.

First-order kinetics with A_r·t_end = 16 ⇒ residual conversion e⁻¹⁶ ≈ 1e-7
bounds the endpoint error; the depletion limiter is disabled so the decay
stays exactly first-order (with it, the tail below ~1 kg/m³ slows and the
residual at fixed t_end grows ~1000×).

### Quantities of interest (n = 8)

| QoI | value | exact | error | tolerance | within tol | provenance |
| --- | --- | --- | --- | --- | --- | --- |
| P1: final thickness ratio | 2 | 2 | 7.102e-07 | 3e-06 | yes | path-independent endpoint of θ = Σγⱼω̇ⱼ/ρⱼ (case header) |
| P1: total mass conservation | 2 | 2 | 6.539e-07 | 2e-06 | yes | solid-only permutation; impermeable everywhere |
| P2: final thickness ratio | 0.5 | 0.5 | 2.981e-08 | 3e-06 | yes | path-independent endpoint of θ = Σγⱼω̇ⱼ/ρⱼ (case header) |
| P2: total mass conservation | 2 | 2 | 1.423e-07 | 2e-06 | yes | solid-only permutation; impermeable everywhere |
| P3: final thickness ratio | 1 | 1 | 0 | 3e-06 | yes | path-independent endpoint of θ = Σγⱼω̇ⱼ/ρⱼ (case header) |
| P3: total mass conservation | 2 | 2 | 0 | 2e-06 | yes | solid-only permutation; impermeable everywhere |
| P4: final thickness ratio | 0.3 | 0.3 | 1.571e-07 | 3e-06 | yes | path-independent endpoint of θ = Σγⱼω̇ⱼ/ρⱼ (case header) |
| P4: total mass conservation | 2 | 2 | 1.055e-07 | 2e-06 | yes | held gas compresses under collapse (flux form); the P4 corner is the Finding-2 regression anchor |
| P4: final char mass (solid ledger) | 0.6 | 0.6 | 2.18e-07 | 2e-06 | yes | char endpoint is path-independent; complements the total-mass anchor |
| P5: final thickness ratio | 1 | 1 | 0 | 3e-06 | yes | path-independent endpoint of θ = Σγⱼω̇ⱼ/ρⱼ (case header) |
| P5: total mass conservation | 2 | 2 | 6.661e-16 | 2e-06 | yes | held gas compresses under collapse (flux form); the P4 corner is the Finding-2 regression anchor |
| P5: final char mass (solid ledger) | 1 | 1 | 1.125e-07 | 2e-06 | yes | char endpoint is path-independent; complements the total-mass anchor |
| P6: final thickness ratio | 0.5 | 0.5 | 2.981e-08 | 3e-06 | yes | path-independent endpoint of θ = Σγⱼω̇ⱼ/ρⱼ (case header) |
| P6: total mass conservation | 2 | 2 | 1.423e-07 | 2e-06 | yes | solid-only permutation; impermeable everywhere |
| P7: final thickness ratio | 1.333 | 1.333 | 6.951e-08 | 3e-06 | yes | path-independent endpoint of θ = Σγⱼω̇ⱼ/ρⱼ (case header) |
| P7: total mass conservation | 2 | 2 | 1.979e-08 | 2e-06 | yes | held gas compresses under collapse (flux form); the P4 corner is the Finding-2 regression anchor |
| P7: final char mass (solid ledger) | 1 | 1 | 1.323e-07 | 2e-06 | yes | char endpoint is path-independent; complements the total-mass anchor |
| P7: final lateral area ratio A_f/A₀ = law(χ̄ = 1/2) | 0.75 | 0.75 | 4.972e-08 | 3e-06 | yes | χ̄_f = released dry-mass fraction = y_gas = 0.5 |

![V-CPL-03 thickness_histories](assets/V-CPL-03/V-CPL-03_thickness_histories.png)

![V-CPL-03 thickness_histories_error](assets/V-CPL-03/V-CPL-03_thickness_histories_error.png)

### Comparison with other codes

The same case was solved with FDS 6.11.0 at 4 cells (2 mm slab; stretched); decks, outputs, and
run provenance are committed under `test/verification/reference/`.
Each code's error against the same exact solution is drawn below
on a log scale, muted gray behind this solver's series — the
signed linear-scale panel above shows where the error lives, this
one compares magnitudes across codes.

![V-CPL-03 thickness_histories_vs_codes_error](assets/V-CPL-03/V-CPL-03_thickness_histories_vs_codes_error.png)

Wall times at every ladder rung against the reference runs.
Resolutions and simulated spans differ where noted (details in
`reference/timings.csv`), so cross-code timings are indicative
rather than a controlled benchmark; rungs at a matched
resolution are directly comparable.

![V-CPL-03 solve_time](assets/V-CPL-03/V-CPL-03_solve_time.png)

Solution overlays including the other codes' points: [thickness_histories_vs_codes](assets/V-CPL-03/V-CPL-03_thickness_histories_vs_codes.png).

### Solver configuration

| setting | value |
| --- | --- |
| integrator | KenCarp4 (default) |
| use_ale | true |
| min_thickness | 1.0e-5 |
| abstol | 1.0e-10 |
| reltol | 1.0e-8 |

