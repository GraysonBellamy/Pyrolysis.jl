#!/usr/bin/env python3
"""
Derivation of the White-Pine variable-cross-section (shrinkage) calibration —
the VOLUME-grounded re-derivation that supersedes the area-only fit documented
in the older `WP_SR2` block of wp_properties.jl.

WHY THIS EXISTS
---------------
The previous calibration fit a single empirical *radial* area law
    A_f/A_0 = 1 - (1 - s_r^2) * chi_bar
to three diameter anchors and declared the axial (thickness) change "out of
scope." But the CAPA tests give mass, diameter AND thickness per flux, so we
actually measure the full volume V = (pi/4) D^2 L and the residue bulk density
rho = m/V. That turns shrinkage from one fitted area law into a closed,
checkable system:

    * VOLUME is not free — the solver derives it from component densities via
      dV/V = sum_j (gamma_j / rho_j) dxi_j   (src/Physics/volume_change.jl).
      So V_f is a PREDICTION we validate against m/rho, not a fit target.
    * The RADIAL share (A_f/A_0) is the only empirically fitted partition.
    * The AXIAL change is then EMERGENT: L_f/L_0 = (V_f/V_0) / (A_f/A_0).
      It is not a separate model and not "out of scope"; it falls out of the
      volume balance. Char that swells axially at high flux is simply the case
      where radial shrinkage outruns the (small) volume loss.

This script:
  1. Reduces the measured (D, L, m) anchors to A_f/A_0, L_f/L_0, V_f/V_0,
     chi_bar_exp = 1 - m_f/m_0, and the pre/post bulk densities.
  2. Validates the density model: the pre-test triangle (m_0, D_0, L_0) must
     reproduce rho_virgin, and the (near-fully-charred) 70 kW residue density
     must reproduce rho_char — an INDEPENDENT check of WP_RHO_CHAR that the
     dilatometry-based value never had.
  3. Re-fits the radial partition law against measured A_f/A_0 vs the
     model-free chi_bar_exp (slope-through-origin LSQ). The linear law's RMS
     already sits at the diameter-measurement noise floor (+/-0.013), so the
     2-parameter saturating form only buys an overfit (RMS far below noise on
     3 points / 1 residual DOF) — the parsimonious linear law stands for the
     ENDPOINTS. Its faster-than-linear saturation is a real but
     endpoint-unconstrained TRAJECTORY signal (see [5]).
  4. Decomposes the strains: radial LINEAR strain is ~flux-independent (~7%);
     ALL the flux dependence lives in the axial direction (shrink at 30 kW ->
     swell at 70 kW = char intumescence). This is why a single radial law +
     density-driven volume is the right factoring of the data.

Mirrors the convention of derive_cone_view_factor.py / derive_surface_emissivity.py:
self-contained, reproducible, regression-asserted.

Run:  uv run --with numpy python derive_shrinkage.py
"""

import numpy as np

# ===========================================================================
# Component densities under test (kg/m^3) — must match wp_properties.jl
# ===========================================================================
WP_RHO_VIRGIN = 357.5
WP_RHO_CHAR = 89.2
RHO_CHAR_UNCERT = 0.06          # +/-6% measurement uncertainty noted in wp_properties.jl

# ===========================================================================
# Measured CAPA anchors (UMD CAPA N2, dry, multi-replicate means)
# ---------------------------------------------------------------------------
# Pre-test diameter D_0 is the per-flux replicate mean; pre-test thickness L_0
# is the nominal 8.5 mm sample (the same value wp_forward.jl uses as
# `sample_thickness`). The pre-test triangle is cross-checked below: with
# L_0 = 8.5 mm, m_0/(rho_virgin * A_0) reproduces ~8.5 mm, confirming L_0.
# ===========================================================================
L0_MM = 8.5  # nominal pre-test thickness (m)

ANCHORS = {
    # flux : (D0_mm, Df_mm, Lf_mm, m0_g, mf_over_m0)
    30: dict(D0=69.51, Df=65.0, Lf=7.60, m0=11.59, mratio=0.477),
    50: dict(D0=69.43, Df=64.4, Lf=8.80, m0=11.93, mratio=0.273),
    70: dict(D0=69.59, Df=64.4, Lf=9.30, m0=11.21, mratio=0.233),
}


def bulk_density(m_g, D_mm, L_mm):
    """Bulk density (kg/m^3) from mass (g) and cylinder geometry (mm)."""
    V_mm3 = np.pi * (D_mm / 2.0) ** 2 * L_mm
    # g / mm^3 -> kg/m^3 :  (1e-3 kg) / (1e-9 m^3) = 1e6
    return (m_g / V_mm3) * 1e6


def reduce_anchors():
    rows = []
    for q, a in ANCHORS.items():
        A0 = np.pi * (a["D0"] / 2.0) ** 2
        Af = np.pi * (a["Df"] / 2.0) ** 2
        Af_A0 = Af / A0
        Lf_L0 = a["Lf"] / L0_MM
        Vf_V0 = Af_A0 * Lf_L0
        chi = 1.0 - a["mratio"]
        rho0 = bulk_density(a["m0"], a["D0"], L0_MM)
        rho_f = bulk_density(a["mratio"] * a["m0"], a["Df"], a["Lf"])
        rows.append(
            dict(q=q, Af_A0=Af_A0, Lf_L0=Lf_L0, Vf_V0=Vf_V0, chi=chi,
                 rho0=rho0, rho_f=rho_f, **a)
        )
    return rows


def binary_char_fraction(rho_f):
    """Rough virgin+char binary estimate of the char mass fraction in the
    residue, from 1/rho_f = Y_c/rho_char + (1-Y_c)/rho_virgin. Ignores
    intermediates, so it is an ILLUSTRATIVE diagnostic only — its job is to
    show the residue density sits between the virgin and char anchors."""
    inv = 1.0 / rho_f
    inv_v = 1.0 / WP_RHO_VIRGIN
    inv_c = 1.0 / WP_RHO_CHAR
    return (inv - inv_v) / (inv_c - inv_v)


def fit_linear_through_origin(chi, y):
    """Slope-through-origin LSQ: y = b * chi, with b = sum(chi*y)/sum(chi^2).
    Returns (s_r2 = 1 - b, b, rms_residual)."""
    chi = np.asarray(chi)
    y = np.asarray(y)
    b = np.sum(chi * y) / np.sum(chi * chi)
    resid = y - b * chi
    rms = float(np.sqrt(np.mean(resid ** 2)))
    return 1.0 - b, b, rms


def fit_saturating(chi, AfA0):
    """2-parameter saturating form A/A0 = s + (1-s) exp(-chi/c). Grid-search c
    (with s solved in closed form per c by linear LSQ). With 3 points and 2
    parameters this near-interpolates (1 residual DOF), so the returned RMS
    drops well below the diameter noise floor — an overfit, not evidence the
    form is needed."""
    chi = np.asarray(chi)
    AfA0 = np.asarray(AfA0)
    best = None
    for c in np.geomspace(0.02, 100.0, 4000):
        e = np.exp(-chi / c)
        # A/A0 = s + (1-s) e = e + s (1 - e)  -> linear in s
        x = 1.0 - e
        # minimize || (AfA0 - e) - s*x ||
        s = np.sum(x * (AfA0 - e)) / np.sum(x * x)
        resid = (e + s * x) - AfA0
        rms = float(np.sqrt(np.mean(resid ** 2)))
        if best is None or rms < best[2]:
            best = (s, c, rms)
    return best  # (s, c, rms)


def main():
    rows = reduce_anchors()

    print("=" * 78)
    print("White-Pine shrinkage — VOLUME-grounded calibration")
    print("=" * 78)

    # ---- 1. Reduced measurements -----------------------------------------
    print("\n[1] Reduced CAPA anchors (measured)")
    print(f"  {'flux':>4} {'D0_mm':>6} {'Df_mm':>6} {'Lf_mm':>6} "
          f"{'A_f/A0':>7} {'L_f/L0':>7} {'V_f/V0':>7} {'chi_exp':>8} "
          f"{'rho0':>6} {'rho_f':>6}")
    for r in rows:
        print(f"  {r['q']:>4} {r['D0']:>6.2f} {r['Df']:>6.2f} {r['Lf']:>6.2f} "
              f"{r['Af_A0']:>7.4f} {r['Lf_L0']:>7.4f} {r['Vf_V0']:>7.4f} "
              f"{r['chi']:>8.3f} {r['rho0']:>6.1f} {r['rho_f']:>6.1f}")

    # ---- 2. Density-model validation -------------------------------------
    print("\n[2] Density-model validation (kg/m^3)")
    rho0_mean = np.mean([r["rho0"] for r in rows])
    print(f"  pre-test bulk rho_0  (mean) = {rho0_mean:7.1f}   "
          f"vs WP_RHO_VIRGIN = {WP_RHO_VIRGIN:.1f}   "
          f"(Delta {100*(rho0_mean-WP_RHO_VIRGIN)/WP_RHO_VIRGIN:+.1f}%)")
    print("    -> confirms the pre-test (m_0, D_0, L_0=8.5 mm) triangle and L_0.")
    r70 = next(r for r in rows if r["q"] == 70)
    drho = 100 * (r70["rho_f"] - WP_RHO_CHAR) / WP_RHO_CHAR
    print(f"  70 kW residue rho_f         = {r70['rho_f']:7.1f}   "
          f"vs WP_RHO_CHAR  = {WP_RHO_CHAR:.1f}   (Delta {drho:+.1f}%)")
    print(f"    -> INDEPENDENT confirmation of char density "
          f"(within the +/-{100*RHO_CHAR_UNCERT:.0f}% char-rho uncertainty);")
    print("       the 70 kW sample is ~fully charred, so m_f/V_f *is* rho_char.")
    print("  partial-conversion residue, rough virgin+char binary Y_char:")
    for r in rows:
        Yc = binary_char_fraction(r["rho_f"])
        print(f"    {r['q']:>2} kW : rho_f={r['rho_f']:6.1f}  "
              f"-> Y_char ~ {Yc:5.2f}  (chi_exp={r['chi']:.3f})")
    print("    (illustrative only — ignores intermediates; shows rho_f sits")
    print("     between the virgin and char anchors, as it must.)")

    # ---- 3. Volume-prediction cross-check --------------------------------
    print("\n[3] Volume prediction  V_f/V_0 = (m_f/m_0)(rho_0/rho_f)")
    print(f"  {'flux':>4} {'V_f/V0 meas':>12} {'(mf/m0)(rho0/rhof)':>20}")
    for r in rows:
        pred = r["mratio"] * (r["rho0"] / r["rho_f"])
        print(f"  {r['q']:>4} {r['Vf_V0']:>12.4f} {pred:>20.4f}")
    print("  (identity check — confirms the reduction is internally consistent;")
    print("   in the solver rho_f is set by the density model, so V_f follows.)")

    # ---- 4. Radial partition law refit -----------------------------------
    chi = [r["chi"] for r in rows]
    AfA0 = [r["Af_A0"] for r in rows]
    y = [1.0 - v for v in AfA0]
    s_r2, b, rms = fit_linear_through_origin(chi, y)
    print("\n[4] Radial partition law  A_f/A_0 = 1 - (1 - s_r^2) * chi_bar")
    print(f"  slope-through-origin LSQ on (chi_exp, 1 - A_f/A0):")
    print(f"    s_r^2 = {s_r2:.4f}   (1 - s_r^2 = {b:.4f})   RMS = {rms:.4f}")
    print("  per-anchor s_r^2 = 1 - (1 - A_f/A0)/chi_exp:")
    for r in rows:
        s_i = 1.0 - (1.0 - r["Af_A0"]) / r["chi"]
        print(f"    {r['q']:>2} kW : s_r^2 = {s_i:.4f}")
    s_sat, c_sat, rms_sat = fit_saturating(chi, AfA0)
    print(f"  saturating form A/A0 = s + (1-s)exp(-chi/c): "
          f"best c = {c_sat:.2f}, s = {s_sat:.4f}, RMS = {rms_sat:.4f}")
    print(f"    -> RMS {rms_sat:.4f} << diameter noise (~0.013): a 2-param")
    print("       overfit on 3 points (1 residual DOF), NOT a justified gain.")
    print("       Linear's RMS already sits AT the noise floor -> p = 1 stands")
    print("       for endpoints. (The flat data does hint at fast radial")
    print("       saturation, but the rate is unconstrained by endpoints; see")
    print("       [5] for why radial is ~flux-independent.)")

    # ---- 5. Strain decomposition (the anisotropy story) ------------------
    print("\n[5] Strain decomposition (positive = shrink)")
    print(f"  {'flux':>4} {'eps_radial(lin)':>15} {'eps_axial':>10} {'eps_vol':>9}")
    for r in rows:
        eps_r = 1.0 - np.sqrt(r["Af_A0"])   # linear radial strain
        eps_a = 1.0 - r["Lf_L0"]            # axial strain (neg = swell)
        eps_v = 1.0 - r["Vf_V0"]
        print(f"  {r['q']:>4} {eps_r:>15.4f} {eps_a:>10.4f} {eps_v:>9.4f}")
    print("  -> radial LINEAR strain ~constant (~7%, barely flux-dependent);")
    print("     axial carries ALL the flux dependence: +shrink at 30 kW ->")
    print("     -swell at 70 kW (char intumescence). Volume loss falls with")
    print("     flux. Hence: fit radial once, let density drive volume, and")
    print("     read axial off the balance.")

    # ---- regression asserts ----------------------------------------------
    assert abs(s_r2 - 0.8005) < 0.01, f"s_r^2 drifted: {s_r2}"
    assert abs(rho0_mean - WP_RHO_VIRGIN) / WP_RHO_VIRGIN < 0.05, "rho_0 off"
    assert abs(r70["rho_f"] - WP_RHO_CHAR) / WP_RHO_CHAR < RHO_CHAR_UNCERT, \
        "70 kW residue density inconsistent with char"
    print("\nAll regression checks passed.")
    print("=" * 78)


if __name__ == "__main__":
    main()
