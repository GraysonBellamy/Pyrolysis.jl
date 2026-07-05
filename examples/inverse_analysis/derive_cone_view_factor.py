#!/usr/bin/env python3
"""
Derivation of WP_F_CONE — the sample->cone-heater view factor used in the
white-pine top-surface radiative BC (wp_forward.jl `build_bcs`, wp_properties.jl
`WP_F_CONE`).

Method: contour-integral (Stokes-theorem) view-factor model for the cone
calorimeter HELICAL-COIL (HC) heater from

    P. E. Pinto, J. Valdivia, A. Singh, X. Xi, J. Cuevas, J. L. Urban,
    "Investigating radiation heat transfer from the cone calorimeter heater:
     A new view factor model and uncertainty quantification,"
    Int. J. Heat Mass Transfer 245 (2025) 126976.
    https://doi.org/10.1016/j.ijheatmasstransfer.2025.126976
    Code: https://github.com/javaldivia/Cone_Calorimeter_View_Factor  (MIT)

The HC kernel `vf_dA_to_HC` below is reproduced verbatim from that repository
(VF_dA_to_HC.py, (c) 2025 Valdivia/Pinto/Singh/Xi/Cuevas/Urban, MIT License) so
this derivation is self-contained and does not depend on an external clone.

Geometry (MaCFP / CAPA white pine):
    - circular sample, D = 70 mm  -> R = 35 mm
    - sample sits 40 mm below the BOTTOM of the cone base plate
    - repo convention: d = d1 + d2, d1 = sample-to-baseplate-bottom = 40 mm,
      d2 = base-plate thickness = 5 mm  ->  d = 45 mm
    - standard ISO 5660 / ASTM E1354 cone: R1 = R2 = 40 mm, H = 65 mm

Two view factors come out:
    F_center   = local VF at r=0  (== the Schmidt-Boelter gauge VF, F_gauge):
                 governs the calibration identity q_ext = F_gauge*eps_c*sigma*T_c^4
                 and is the surroundings-split fraction for a CENTERLINE 1D column.
    F_avg      = area-average of the local VF over the round sample disk:
                 the fraction for an AREA-AVERAGE (whole-sample MLR) interpretation.

WP_F_CONE is set to F_center (the inverse fit targets the center back-face TC, so
the 1D column is the centerline). The ratio F_avg/F_gauge ~ 0.96 quantifies the
~4% radial flux non-uniformity (sample overhangs the cone aperture); it is NOT
folded into q_ext here — that is a separate modeling decision (see the comment
block at WP_F_CONE in wp_properties.jl).

Run:  uv run --with numpy --with scipy python derive_cone_view_factor.py
"""

import warnings
warnings.filterwarnings("ignore")  # scipy dblquad emits benign tolerance warnings
import numpy as np
from scipy.integrate import dblquad

# --- ISO 5660 / ASTM E1354 cone heater dimensions (mm) ----------------------
R1 = 40.0
R2 = 40.0
H = 65.0


def vf_dA_to_HC(h, k, d, p):
    """View factor from a differential disk element at (h, k), a vertical
    distance d below the cone reference plane, to the helical-coil heater.
    `p` is the element diameter (mm); the result is p-independent as p -> 0.
    Verbatim kernel from Pinto et al. (2025), VF_dA_to_HC.py (MIT)."""

    def f1(u, v):
        return (1 / (np.pi**2 * p**2)) * np.log(
            (0.5 * p * np.cos(u) + h - (-R1 - R2 + R2 * v / (19 * np.pi)) * np.cos(v))**2
            + (0.5 * p * np.sin(u) + k - (-R1 - R2 + R2 * v / (19 * np.pi)) * np.sin(v))**2
            + (H * v / (19 * np.pi) + d)**2
        ) * (p * R2 / (38 * np.pi) * np.sin(v - u)
             + 0.5 * p * (-R1 - R2 + (R2 * v) / (19 * np.pi)) * np.cos(v - u))

    f2 = f1  # same integrand, different v-range (see ranges below)

    def f3(u, v):
        return (1 / (np.pi**2 * p**2)) * np.log(
            (0.5 * p * np.cos(u) + h - (np.sqrt(H**2 + R2**2) / 19) * np.cos(v) + R1 + 18 * R2 / 19)**2
            + (0.5 * p * np.sin(u) + k)**2
            + ((np.sqrt(H**2 + R2**2) / 19) * np.sin(v) + H / 19 + d)**2
        ) * (p / 38) * np.sqrt(H**2 + R2**2) * np.sin(u) * np.sin(v)

    def f4(u, v):
        return (1 / (np.pi**2 * p**2)) * np.log(
            (0.5 * p * np.cos(u) + h - (np.sqrt(H**2 + R2**2) / 19) * np.cos(v) - R1 - R2 / 19)**2
            + (0.5 * p * np.sin(u) + k)**2
            + ((np.sqrt(H**2 + R2**2) / 19) * np.sin(v) + 18 * H / 19 + d)**2
        ) * (p / 38) * np.sqrt(H**2 + R2**2) * np.sin(v) * np.sin(u)

    return (dblquad(f1, 17 * np.pi, 0, 0, 2 * np.pi)[0]
            + dblquad(f2, 2 * np.pi, 19 * np.pi, 0, 2 * np.pi)[0]
            + dblquad(f3, np.arctan(H / R2) - np.pi, np.arctan(H / R2), 0, 2 * np.pi)[0]
            + dblquad(f4, np.pi - np.arctan(H / R2), 2 * np.pi - np.arctan(H / R2), 0, 2 * np.pi)[0])


def vf_local(r, d, p=1e-3):
    """Local (differential) view factor at radius r, 4-azimuth mean to wash out
    the helix's non-axisymmetric phase."""
    if r == 0.0:
        return vf_dA_to_HC(0.0, 0.0, d, p)
    pts = [(r, 0.0), (0.0, r), (-r, 0.0), (0.0, -r)]
    return float(np.mean([vf_dA_to_HC(hh, kk, d, p) for (hh, kk) in pts]))


def sample_average(d, R, nr=21):
    """Area-average of the local VF over a disk of radius R:
    F_avg = (2/R^2) integral_0^R VF(r) r dr."""
    rs = np.linspace(0.0, R, nr)
    vf = np.array([vf_local(r, d) for r in rs])
    return 2.0 / R**2 * np.trapezoid(vf * rs, rs), rs, vf


if __name__ == "__main__":
    R_SAMPLE = 35.0     # mm (70 mm diameter)
    GAUGE = 3.175       # mm, Schmidt-Boelter sensor diameter
    d = 45.0            # mm = 40 (below baseplate bottom) + 5 (plate thickness)

    # Sanity check against the paper's published examples (d=30 == std 25 mm gap)
    assert abs(vf_dA_to_HC(0, 0, 30, 1e-5) - 0.6690) < 5e-4, "HC kernel regressed"

    f_gauge = vf_dA_to_HC(0.0, 0.0, d, GAUGE)
    f_center = vf_local(0.0, d)
    f_avg, rs, vf = sample_average(d, R_SAMPLE)

    print(f"HC model | 70 mm round sample | d = {d:.0f} mm "
          f"(40 mm below baseplate bottom + 5 mm plate)")
    print(f"  F_gauge  (SB, center)              = {f_gauge:.4f}")
    print(f"  F_center (centerline 1D column)    = {f_center:.4f}   <- WP_F_CONE")
    print(f"  F_avg    (area-avg over disk)      = {f_avg:.4f}")
    print(f"  F_avg / F_gauge  (flux non-unif.)  = {f_avg / f_gauge:.3f}  "
          f"(~{100 * (1 - f_avg / f_gauge):.0f}% below gauge nominal)")
