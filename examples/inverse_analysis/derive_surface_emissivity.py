#!/usr/bin/env python3
"""
Spectral derivation of white-pine radiative surface properties from the FSRI
integrating-sphere reflectivity spectra (data/FSRI_Wood*_REFLECTIVITY_*.csv).

Each file is room-temperature directional-hemispherical reflectivity ρ(ν̃) vs
wavenumber ν̃ [cm⁻¹], 600–7497 cm⁻¹ (≈ 1.33–16.7 µm). Opaque sample ⇒ spectral
absorptance = spectral emittance: α(λ) = ε(λ) = 1 − ρ(λ) (Kirchhoff, per λ).

Two band-integrated quantities come from the SAME spectrum via different Planck
weightings (this is the whole point — they are NOT the same number):

  ABSORPTIVITY for cone radiation : weight by Planck at the cone effective
      temperature T_cone (~940–1160 K, near-IR).  → used in  ε·q_ext
  EMISSIVITY for surface emission : weight by Planck at the surface temperature
      T_s (~300–1000 K, mid/far-IR).              → used in  −ε·σ·T_s⁴

    ε_band(T) = ∫ ε(ν̃) B(ν̃,T) dν̃ / ∫ B(ν̃,T) dν̃,
    B(ν̃,T) ∝ ν̃³ / (exp(c2·ν̃/T) − 1),  c2 = 1.438777 cm·K.

ε(λ) is held flat outside the measured band (the far-IR tail matters only below
~500 K, where emission is small; band-coverage fraction is reported).

Spectra groups:
  virgin                      : FSRI_Wood_REFLECTIVITY_R*            (n=3)
  char formed in CAPA at q_f  : FSRI_Wood_Char_CAPA{20,40,60}_*      (char optical
      properties depend on formation flux q_f). The production model pools 40+60
      (n=11) as its single char table.

Run:  uv run --with numpy python derive_surface_emissivity.py
"""

import glob
import numpy as np

C2 = 1.438777  # second radiation constant [cm·K] for ν̃ in cm⁻¹

# Published tables currently in wp_properties.jl (for the validation assert)
PUB_VIRGIN = {30000: 0.838, 50000: 0.807, 70000: 0.783}
PUB_CHAR = {30000: 0.822, 50000: 0.831, 70000: 0.838}  # 40+60 pooled


def load_group(pattern):
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(pattern)
    arrs = [np.loadtxt(f, delimiter=",", skiprows=1) for f in files]
    nu = arrs[0][:, 0]
    for a in arrs:
        assert np.allclose(a[:, 0], nu), f"grid mismatch in {pattern}"
    return nu, np.mean([a[:, 1] for a in arrs], axis=0), len(files)


def planck_band_average(nu_d, R_d, T):
    """ε band-averaged at temperature T, with flat ε extrapolation outside the
    measured band. Returns (ε_band, in-band blackbody-energy fraction)."""
    eps_d = 1.0 - R_d
    lo = np.arange(50.0, nu_d[0])
    hi = np.arange(nu_d[-1] + 1.0, 12000.0)
    nu = np.concatenate([lo, nu_d, hi])
    eps = np.concatenate([np.full_like(lo, eps_d[0]), eps_d, np.full_like(hi, eps_d[-1])])
    x = C2 * nu / T
    B = np.where(x < 700, nu**3 / np.expm1(np.clip(x, 1e-9, 700)), 0.0)
    den = np.trapezoid(B, nu)
    m = (nu >= nu_d[0]) & (nu <= nu_d[-1])
    return np.trapezoid(eps * B, nu) / den, np.trapezoid(B[m], nu[m]) / den


def recover_cone_temperature(nu_v, R_v, eps_target):
    Ts = np.arange(700.0, 1600.0, 1.0)
    vals = np.array([planck_band_average(nu_v, R_v, T)[0] for T in Ts])
    return Ts[np.argmin(np.abs(vals - eps_target))]


if __name__ == "__main__":
    base = "data/"
    virgin = load_group(base + "FSRI_Wood_REFLECTIVITY_R*.csv")
    char = {q: load_group(base + f"FSRI_Wood_Char_CAPA{q}_REFLECTIVITY_R*.csv")
            for q in (20, 40, 60)}
    char_pool = load_group(base + "FSRI_Wood_Char_CAPA[46]0_REFLECTIVITY_R*.csv")

    # --- Validation: recover cone temps from virgin, check they reproduce char ---
    print("VALIDATION — recover cone T from virgin, confirm it reproduces char")
    Tcone = {}
    for q, e in PUB_VIRGIN.items():
        T = recover_cone_temperature(virgin[0], virgin[1], e)
        Tcone[q] = T
        v = planck_band_average(virgin[0], virgin[1], T)[0]
        c = planck_band_average(char_pool[0], char_pool[1], T)[0]
        print(f"  q={q//1000:>2}k  T_cone={T:4.0f}K  virgin={v:.3f}(pub {e:.3f})  "
              f"char={c:.3f}(pub {PUB_CHAR[q]:.3f})")
        assert abs(v - e) < 3e-3 and abs(c - PUB_CHAR[q]) < 3e-3, "spectra do not reproduce tables"
    print("  ✓ both published tables reproduced from raw spectra\n")

    # --- Emission emissivity ε(T_s) ---
    print("EMISSION EMISSIVITY ε(T_s)  [absorptivity tables for comparison: "
          "virgin 0.78–0.84, char 0.82–0.84]")
    print(f"  {'T_s':>5} {'virgin':>7} {'char20':>7} {'char40':>7} {'char60':>7} "
          f"{'char40+60':>9} {'band%':>6}")
    for Ts in (300, 400, 500, 600, 700, 800, 900, 1000):
        ev, fr = planck_band_average(*virgin[:2], Ts)
        ec = [planck_band_average(*char[q][:2], Ts)[0] for q in (20, 40, 60)]
        ep = planck_band_average(*char_pool[:2], Ts)[0]
        print(f"  {Ts:>5} {ev:7.3f} {ec[0]:7.3f} {ec[1]:7.3f} {ec[2]:7.3f} "
              f"{ep:9.3f} {100*fr:5.0f}%")

    # --- Char absorptivity vs formation flux (model uses 40+60 at all q) ---
    print("\nCHAR ABSORPTIVITY @ cone T — sensitivity to FORMATION flux q_f")
    print(f"  {'test q':>6} {'char20':>7} {'char40':>7} {'char60':>7} {'used':>7}")
    for q, T in Tcone.items():
        ec = [planck_band_average(*char[qf][:2], T)[0] for qf in (20, 40, 60)]
        used = planck_band_average(*char_pool[:2], T)[0]
        print(f"  {q//1000:>4}k  {ec[0]:7.3f} {ec[1]:7.3f} {ec[2]:7.3f} {used:7.3f}")
