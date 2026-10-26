#!/usr/bin/env python3
"""Extract SI CSVs from four Gpyro verification runs and check exact solutions.

Run as `python3 extract_results.py WORKSPACE`, with raw runs under
WORKSPACE/gpyro_run. Results are written to WORKSPACE/gpyro_out.
"""
import argparse
import csv, math, os, shutil

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("workspace", help="directory containing gpyro_run")
workspace = parser.parse_args().workspace
RUN = os.path.join(workspace, "gpyro_run")
OUT = os.path.join(workspace, "gpyro_out")
os.makedirs(OUT, exist_ok=True)
C2K = 273.15

def read_profile(path):
    rows = [r for r in csv.reader(open(path))]
    z = [float(v) for v in rows[0][1:] if v.strip()]
    data = []
    for r in rows[1:]:
        t = float(r[0])
        T = [float(v) + C2K for v in r[1:] if v.strip()]
        assert len(T) == len(z), (len(T), len(z))
        data.append((t, T))
    return z, data

def read_summary(path):
    rows = [r for r in csv.reader(open(path))]
    hdr = [h.strip() for h in rows[0] if h.strip()]
    data = [[float(v) for v in r if v.strip() != ""] for r in rows[1:]]
    return hdr, data

def nearest(data, t_target):
    return min(data, key=lambda p: abs(p[0] - t_target))

report = []

# ---------------------------------------------------------------- V-HC-01
# exact: erfc solution, absorbed flux + convection == convection to Tinf_eff
def erfc(x): return math.erfc(x)
k, rho, c, q, h, T0 = 0.2, 1000.0, 1400.0, 25e3, 20.0, 300.0
alpha = k / (rho * c)
Tinf_eff = T0 + q / h
def T_hc01(x, t):
    if t <= 0: return T0
    s = math.sqrt(alpha * t)
    eta = x / (2 * s)
    arg = h * s / k
    # erfc(eta) - exp(h x/k + h^2 alpha t/k^2) erfc(eta + arg), scaled
    term2 = math.exp(min(h * x / k + arg * arg, 700)) * erfc(eta + arg)
    return T0 + (Tinf_eff - T0) * (erfc(eta) - term2)

z, prof = read_profile(f"{RUN}/vhc01/vhc01_profile_01_01_TEMPERATURE.csv")
times = [30.0, 60.0, 90.0, 180.0]
with open(f"{OUT}/vhc01_result.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["t_s", "x_m", "T_K"])
    for tt in times:
        t, T = nearest(prof, tt)
        err = max(abs(Ti - T_hc01(xi, t)) for xi, Ti in zip(z, T))
        report.append(f"V-HC-01  t={t:7.2f} s: Linf |Gpyro - exact erfc| = {err:.3f} K "
                      f"(T_s gpyro {T[0]:.2f} K, exact {T_hc01(0,t):.2f} K)")
        for xi, Ti in zip(z, T):
            w.writerow([f"{t:.4f}", f"{xi:.6e}", f"{Ti:.5f}"])
hdr, s = read_summary(f"{RUN}/vhc01/vhc01_summary_01_0001.csv")
with open(f"{OUT}/vhc01_surfaceT.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["t_s", "T_surface_K", "T_back_K"])
    for r in s:
        w.writerow([f"{r[0]:.4f}", f"{r[1]+C2K:.5f}", f"{r[2]+C2K:.5f}"])

# ---------------------------------------------------------------- V-RAD-01
delta, kappa = 0.01, 1000.0
def rad_eigenvalues(n):
    # u*tan(u) = Bi, u = lam*delta, Bi = h*delta/k; roots in (m*pi, m*pi + pi/2)
    Bi = h * delta / k
    us = []
    for m in range(n):
        lo, hi = m * math.pi + 1e-12, m * math.pi + math.pi / 2 - 1e-12
        f = lambda u: u * math.tan(u - m * math.pi) - Bi  # tan periodic
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if f(lo) * f(mid) <= 0: hi = mid
            else: lo = mid
        us.append(0.5 * (lo + hi))
    return us

# Exact solution per test/verification/exact/indepth_absorption.jl:
# theta_ss + eigen-series with coefficients from Simpson projection.
US = rad_eigenvalues(120)
LAMS = [u / delta for u in US]
def theta_ss(x):
    return q * ((1 - math.exp(-kappa * delta)) / h + (1 - math.exp(-kappa * x)) / (k * kappa)) \
        - (q / k) * math.exp(-kappa * delta) * x
_n = 4096
_dx = delta / _n
_xs = [i * _dx for i in range(_n + 1)]
_w = [1.0 if i in (0, _n) else (4.0 if i % 2 == 1 else 2.0) for i in range(_n + 1)]
_th = [theta_ss(x) for x in _xs]
CS = []
for m, lamm in enumerate(LAMS):
    num = _dx / 3 * sum(_w[i] * _th[i] * math.cos(lamm * (delta - _xs[i])) for i in range(_n + 1))
    den = delta / 2 + math.sin(2 * US[m]) / (4 * lamm)
    CS.append(num / den)
def T_rad01(x, t):
    if t <= 0: return T0
    th = theta_ss(x)
    for cm, lamm in zip(CS, LAMS):
        th -= cm * math.exp(-alpha * lamm * lamm * t) * math.cos(lamm * (delta - x))
    return T0 + th

z, prof = read_profile(f"{RUN}/vrad01/vrad01_profile_01_01_TEMPERATURE.csv")
with open(f"{OUT}/vrad01_result.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["t_s", "x_m", "T_K"])
    for tt in times:
        t, T = nearest(prof, tt)
        err = max(abs(Ti - T_rad01(xi, t)) for xi, Ti in zip(z, T))
        ipk = max(range(len(T)), key=lambda i: T[i])
        report.append(f"V-RAD-01 t={t:7.2f} s: Linf |Gpyro - exact eigenfn| = {err:.3f} K "
                      f"(peak {T[ipk]:.2f} K at x={z[ipk]*1e3:.2f} mm; exact T(0)={T_rad01(0,t):.2f} K, gpyro {T[0]:.2f} K)")
        for xi, Ti in zip(z, T):
            w.writerow([f"{t:.4f}", f"{xi:.6e}", f"{Ti:.5f}"])
hdr, s = read_summary(f"{RUN}/vrad01/vrad01_summary_01_0001.csv")
with open(f"{OUT}/vrad01_surfaceT.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["t_s", "T_front_K", "T_back_K"])
    for r in s:
        w.writerow([f"{r[0]:.4f}", f"{r[1]+C2K:.5f}", f"{r[2]+C2K:.5f}"])

# ---------------------------------------------------------------- V-CPL-01
ks_, rho_s, cs_ = 0.2, 900.0, 1700.0
Tm, dHm, T0s, Tinf_s = 433.0, 50e3, 438.0, 428.0
alpha_s = ks_ / (rho_s * cs_)   # = alpha_A
def stefan_lambda():
    St_ratio = (T0s - Tm) / (Tm - Tinf_s)
    rhs_c = dHm * math.sqrt(math.pi) / (cs_ * (Tm - Tinf_s))
    f = lambda lam: math.exp(-lam * lam) / math.erf(lam) \
        - St_ratio * math.exp(-lam * lam) / math.erfc(lam) - rhs_c * lam
    lo, hi = 1e-6, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(lo) * f(mid) <= 0: hi = mid
        else: lo = mid
    return 0.5 * (lo + hi)
lam = stefan_lambda()
report.append(f"V-CPL-01 eigenvalue lambda = {lam:.4f} (dissertation: 0.2317)")
def T_stefan(x, t):
    if t <= 0: return T0s
    s = 2 * lam * math.sqrt(alpha_s * t)
    if x <= s:
        return Tinf_s + (Tm - Tinf_s) * math.erf(x / (2 * math.sqrt(alpha_s * t))) / math.erf(lam)
    return T0s + (Tm - T0s) * math.erfc(x / (2 * math.sqrt(alpha_s * t))) / math.erfc(lam)

z, prof = read_profile(f"{RUN}/vcpl01/vcpl01_profile_01_01_TEMPERATURE.csv")
depths = [0.5e-3, 1.0e-3, 1.5e-3, 2.0e-3]
didx = {d: min(range(len(z)), key=lambda i: abs(z[i] - d)) for d in depths}
for d, i in didx.items():
    assert abs(z[i] - d) < 1e-9, (d, z[i])
def front_pos(T):
    for i in range(len(T) - 1):
        if (T[i] - Tm) * (T[i + 1] - Tm) <= 0 and T[i] != T[i + 1]:
            return z[i] + (z[i + 1] - z[i]) * (Tm - T[i]) / (T[i + 1] - T[i])
    return float("nan")
with open(f"{OUT}/vcpl01_result.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["t_s", "x_m", "T_K"])
    for t, T in prof:
        for d in depths:
            w.writerow([f"{t:.4f}", f"{d:.6e}", f"{T[didx[d]]:.5f}"])
with open(f"{OUT}/vcpl01_front.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["t_s", "s_front_m"])
    for t, T in prof:
        s_ = front_pos(T)
        if not math.isnan(s_) and t > 0:
            w.writerow([f"{t:.4f}", f"{s_:.6e}"])
for tt in [100.0, 200.0, 299.0]:
    t, T = nearest(prof, tt)
    s_num = front_pos(T); s_ex = 2 * lam * math.sqrt(alpha_s * t)
    errs = [abs(T[didx[d]] - T_stefan(d, t)) for d in depths]
    report.append(f"V-CPL-01 t={t:7.2f} s: front {s_num*1e3:.3f} mm vs exact {s_ex*1e3:.3f} mm "
                  f"({100*(s_num-s_ex)/s_ex:+.2f}%); |dT| at 0.5/1/1.5/2 mm = "
                  + "/".join(f"{e:.3f}" for e in errs) + " K")


# ---------------------------------------------------------------- V-HC-03
# finite slab, constant net flux q3 at x=0, adiabatic back (pure Neumann).
# Exact (C&J 3.8): T = T0 + (q3/(rho3 c3 L3)) t
#   + (q3 L3/k3) [ (3(L3-x)^2 - L3^2)/(6 L3^2)
#                  - (2/pi^2) sum 1/n^2 exp(-n^2 pi^2 Fo) cos(n pi x / L3) ]
k3, rho3, c3, q3, L3 = 0.5, 1000.0, 1000.0, 5000.0, 0.02
alpha3 = k3 / (rho3 * c3)
def T_hc03(x, t):
    if t <= 0: return T0
    Fo = alpha3 * t / (L3 * L3)
    ser = sum(math.exp(-n * n * math.pi ** 2 * Fo) * math.cos(n * math.pi * x / L3) / (n * n)
              for n in range(1, 400))
    return T0 + q3 * t / (rho3 * c3 * L3) + (q3 * L3 / k3) * (
        (3 * (L3 - x) ** 2 - L3 * L3) / (6 * L3 * L3) - 2 * ser / math.pi ** 2)

z, prof = read_profile(f"{RUN}/vhc03/vhc03_profile_01_01_TEMPERATURE.csv")
with open(f"{OUT}/vhc03_result.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["t_s", "x_m", "T_K"])
    for tt in [80.0, 400.0, 1200.0]:
        t, T = nearest(prof, tt)
        err = max(abs(Ti - T_hc03(xi, t)) for xi, Ti in zip(z, T))
        meanT = sum(T) / len(T)   # uniform nodes incl. half-weight ends -> trapezoid
        meanT_trap = (0.5 * T[0] + sum(T[1:-1]) + 0.5 * T[-1]) / (len(T) - 1)
        report.append(f"V-HC-03  t={t:7.2f} s (Fo={alpha3*t/L3**2:.2f}): Linf |Gpyro - exact series| = {err:.3f} K; "
                      f"mean T (trap) {meanT_trap:.3f} K vs exact {T0 + q3*t/(rho3*c3*L3):.3f} K; "
                      f"front-back dT {T[0]-T[-1]:.3f} K (exact {T_hc03(0,t)-T_hc03(L3,t):.3f}, plateau 100)")
        for xi, Ti in zip(z, T):
            w.writerow([f"{t:.4f}", f"{xi:.6e}", f"{Ti:.5f}"])
hdr, s3 = read_summary(f"{RUN}/vhc03/vhc03_summary_01_0001.csv")
with open(f"{OUT}/vhc03_frontback.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["t_s", "T_front_K", "T_back_K"])
    for r in s3:
        w.writerow([f"{r[0]:.4f}", f"{r[1]+C2K:.5f}", f"{r[2]+C2K:.5f}"])

# ---------------------------------------------------------------- V-RX-01
k1, k2 = 0.389, 0.262
def bateman(t):
    ya = math.exp(-k1 * t)
    yb = k1 / (k2 - k1) * (math.exp(-k1 * t) - math.exp(-k2 * t))
    return ya, yb, 1.0 - ya - yb
hdr, s = read_summary(f"{RUN}/vrx01/vrx01_summary_01_0001.csv")
maxe = [0.0, 0.0, 0.0]; maxm = 0.0; maxTdrift = 0.0
with open(f"{OUT}/vrx01_result.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["t_s", "Y_A", "Y_B", "Y_C"])
    for r in s:
        t, ya, yb, yc, mm0, TC = r[0], r[1], r[2], r[3], r[4], r[5]
        ex = bateman(t)
        for i, y in enumerate((ya, yb, yc)):
            maxe[i] = max(maxe[i], abs(y - ex[i]))
        maxm = max(maxm, abs(mm0 - 1.0))
        maxTdrift = max(maxTdrift, abs(TC + C2K - 400.0))
        w.writerow([f"{t:.4f}", f"{ya:.7f}", f"{yb:.7f}", f"{yc:.7f}"])
report.append(f"V-RX-01  max |Y - Bateman| A/B/C = {maxe[0]:.2e}/{maxe[1]:.2e}/{maxe[2]:.2e}; "
              f"max |M/M0 - 1| = {maxm:.1e}; max |T - 400 K| = {maxTdrift:.1e} K")

# copy decks + run logs
for case in ["vhc01", "vrad01", "vcpl01", "vrx01", "vhc03"]:
    shutil.copy(f"{RUN}/{case}/{case}.data", f"{OUT}/{case}.data")
    shutil.copy(f"{RUN}/{case}/{case}_run.log", f"{OUT}/{case}_run.log")

print("\n".join(report))
