#!/usr/bin/env python3
"""Checks behind README_fds.md. Standard library only.

    python3 fds_checks.py order [SET ...]        non-unit-order approximation (0D TGA)
    python3 fds_checks.py nist                   FDS vs NIST's own MaCFP-2 submission
    python3 fds_checks.py submission [SET ...]   FDS vs each group's own MaCFP-2 submission
    python3 fds_checks.py refine SET VARIANT...  base vs refined runs (same set, cases present)
    python3 fds_checks.py pyrolysis [SET ...]    FDS vs Pyrolysis.jl (results/pyrolysis/<set>/base)
    python3 fds_checks.py runs                   run counts and wall/CPU times from runs_*.csv
    python3 fds_checks.py stale                  kept inputs that the current decks would change
"""

import csv
import math
import sys
from pathlib import Path

import fds

MDB = Path.home() / "Documents" / "git" / "matl-db" / "PMMA" / "Computational_Results" / "Gasification_Predictions"
PYRO = fds.SUITE / "results" / "pyrolysis"


# ---------------------------------------------------------------------------
# Series helpers
# ---------------------------------------------------------------------------

def read_series(path):
    with open(path) as f:
        rows = list(csv.reader(f))
    return {k: [float(r[i]) for r in rows[1:]] for i, k in enumerate(rows[0])}


def read_matl(path):
    """matl-db submission: numeric rows only (header rows vary by institution)."""
    rows = []
    for raw in open(path, encoding="utf-8-sig"):
        f = [x.strip().strip('"') for x in raw.replace(";", ",").split(",")]
        try:
            v = [float(x) for x in f if x != ""]
        except ValueError:
            continue
        if v:
            rows.append(v)
    w = min(len(r) for r in rows)
    cols = list(zip(*[r[:w] for r in rows]))
    return {"time_s": list(cols[0]), "mass_g": list(cols[1]), "T_back_K": list(cols[2]),
            "T_top_K": list(cols[3])}


def interp(t, x, y):
    if t <= x[0]:
        return y[0]
    if t >= x[-1]:
        return y[-1]
    lo, hi = 0, len(x) - 1
    while hi - lo > 1:
        m = (lo + hi) // 2
        if x[m] <= t:
            lo = m
        else:
            hi = m
    w = (t - x[lo]) / (x[hi] - x[lo])
    return y[lo] * (1 - w) + y[hi] * w


def burnout(s, frac=1e-3):
    """First time the mass is below `frac` of its initial value (end if never)."""
    m0 = s["mass_g"][0]
    for t, m in zip(s["time_s"], s["mass_g"]):
        if m < frac * m0:
            return t
    return s["time_s"][-1]


def virgin_end(s, frac=0.01):
    """First time the mass is within `frac` of the initial mass of its final
    value: the end of the burning phase (burnout, or only residue left)."""
    m0, mf = s["mass_g"][0], s["mass_g"][-1]
    for t, m in zip(s["time_s"], s["mass_g"]):
        if m < mf + frac * m0:
            return t
    return s["time_s"][-1]


def diff(a, b, key, tmax):
    """max |b - a| and RMS over a's times <= tmax (b interpolated)."""
    d = [interp(t, b["time_s"], b[key]) - y for t, y in zip(a["time_s"], a[key]) if t <= tmax]
    d = [x for x in d if not math.isnan(x)]
    if not d:
        return math.nan, math.nan
    return max(map(abs, d)), math.sqrt(sum(x * x for x in d) / len(d))


def submission(deck, q, L):
    for kw in ("kW", "kw"):
        p = MDB / f"{deck['predictions']}_Gasification_{q}{kw}_{L}mm.csv"
        if p.exists():
            return p
    return None


# ---------------------------------------------------------------------------
# Non-unit order: one shrinking FDS cell vs the deck's law, 0D TGA
# ---------------------------------------------------------------------------

def _rhs(deck, m, T, form, ref):
    comps = deck["components"]
    ids = [c["id"] for c in comps]
    V = sum(max(mi, 0.0) / c["rho"] for mi, c in zip(m, comps))
    dm = [0.0] * len(ids)
    for r in deck["reactions"]:
        i = ids.index(r["reactant"])
        n = float(r["n"])
        k = math.exp(-float(r["E"]) / (fds.R_JL * T))
        if form == "deck":   # Pyrolysis.jl / Gpyro "own": A M_ref (m/M_ref)^n
            rate = float(r["A"]) * k * ref[i] * (max(m[i], 0.0) / ref[i]) ** n
        else:                # FDS: A_F rho_s^n per cell volume, cell shrinks with sum(m_i/rho_i)
            rate = 0.0 if V <= 0 else fds.fds_A(deck, r) * (max(m[i], 0.0) / V) ** n * k * V
        dm[i] -= rate
        if r.get("residue"):
            dm[ids.index(r["residue"])] += float(r["residue_yield"]) * rate
    return dm


def tga(deck, beta_Kmin, form, dT=0.01):
    rho0 = fds.initial_density(deck)
    comps = deck["components"]
    m = [c["y0"] * rho0 for c in comps]
    ref = [fds.reference_concentration(deck, c["id"]) or 1.0 for c in comps]
    dt = dT / (beta_Kmin / 60.0)
    T = 300.0
    Ts, ms = [T], [1.0]
    f = lambda m, T: _rhs(deck, m, T, form, ref)  # noqa: E731
    while T < 1000.0:  # RK4 in time, T = 300 + beta t
        k1 = f(m, T)
        k2 = f([a + 0.5 * dt * b for a, b in zip(m, k1)], T + 0.5 * dT)
        k3 = f([a + 0.5 * dt * b for a, b in zip(m, k2)], T + 0.5 * dT)
        k4 = f([a + dt * b for a, b in zip(m, k3)], T + dT)
        m = [max(a + dt / 6 * (b + 2 * c + 2 * d + e), 0.0) for a, b, c, d, e in zip(m, k1, k2, k3, k4)]
        T += dT
        Ts.append(T)
        ms.append(sum(m) / rho0)
    return Ts, ms


def tga_summary(Ts, ms, beta):
    rate = [-(ms[i + 1] - ms[i - 1]) / (Ts[i + 1] - Ts[i - 1]) for i in range(1, len(ms) - 1)]
    k = max(range(len(rate)), key=lambda i: rate[i])
    T50 = next(T for T, m in zip(Ts, ms) if m <= 0.5 * (1 + ms[-1]))
    return Ts[k + 1], rate[k] * beta / 60 * 100, T50


def cmd_order(sets):
    decks = [d for d in fds.load_decks(sets, include_all=not sets)
             if any(float(r["n"]) != 1 for r in d["reactions"])]
    print("set       beta  T_peak deck/FDS (K)   T_50% deck/FDS (K)    peak rate deck/FDS (%/s)  max|dm/m0| at equal T")
    for d in decks:
        for beta in (10, 100):
            Td, md = tga(d, beta, "deck")
            Tf, mf = tga(d, beta, "fds")
            pd, rd, hd = tga_summary(Td, md, beta)
            pf, rf, hf = tga_summary(Tf, mf, beta)
            dmax = max(abs(a - b) for a, b in zip(md, mf))
            print(f"{d['id']:9s} {beta:4d}  {pd:6.1f} / {pf:6.1f} ({pf - pd:+.1f})   {hd:6.1f} / {hf:6.1f} ({hf - hd:+.1f})   "
                  f"{rd:.3f} / {rf:.3f} ({100 * (rf / rd - 1):+.0f} %)     {dmax:.3f}")


# ---------------------------------------------------------------------------
# Comparisons
# ---------------------------------------------------------------------------

def cmd_nist(sets=("NIST",)):
    """FDS vs the group's own MaCFP-2 submission. The submitted 'top' column is
    compared with T_1mm for decks whose top temperature is 1 mm deep."""
    for deck in fds.load_decks(sets):
        cmp_submission(deck)


def cmp_submission(deck):
    top = "T_top_K" if deck.get("top_temperature", "surface") == "surface" else "T_1mm_K"
    print(f"{deck['id']}: FDS vs submission '{deck['predictions']}' (top column vs FDS {top})")
    print("case       burnout sub/FDS (s)   mass max/rms (g)   T_top max/rms (K)   T_back max/rms (K)  [to the earlier burnout - 20 s]")
    for q, L, _ in fds.gasification_cases():
        p = submission(deck, q, L)
        f = fds.RESULTS / deck["id"] / "base" / f"gasification_{fds.case_name(q, L)}.csv"
        if not (p and f.exists()):
            continue
        s, r = read_matl(p), read_series(f)
        tb_s, tb_f = burnout(s), burnout(r)
        lim = min(tb_s, tb_f) - 20
        mm = diff(s, r, "mass_g", min(s["time_s"][-1], r["time_s"][-1]))
        r = dict(r, T_top_K=r[top])
        tt = diff(s, r, "T_top_K", lim)
        tb = diff(s, r, "T_back_K", lim)
        print(f"{fds.case_name(q, L):10s} {tb_s:7.0f} / {tb_f:7.0f}   {mm[0]:6.3f} / {mm[1]:6.3f}     "
              f"{tt[0]:6.2f} / {tt[1]:5.2f}      {tb[0]:6.2f} / {tb[1]:5.2f}")


def cmd_refine(set_id, variants):
    base = fds.RESULTS / set_id / "base"
    print(f"{set_id}: change from base (cell {fds.CELL_SIZE} m, dt {fds.DT} s); temperatures to the end of the burning phase")
    print("case       variant          burnout base/var (s)  mass max (g)  mlr max (g/m2s)  T_top max  T_back max  T_1mm max (K)")
    for q, L, _ in fds.gasification_cases():
        name = f"gasification_{fds.case_name(q, L)}.csv"
        for v in variants:
            fv = fds.RESULTS / set_id / v / name
            if not (fv.exists() and (base / name).exists()):
                continue
            a, b = read_series(base / name), read_series(fv)
            lim = min(virgin_end(a), virgin_end(b))
            row = [diff(a, b, "mass_g", a["time_s"][-1])[0]]
            row += [diff(a, b, k, lim)[0] for k in ("mlr_g_m2s", "T_top_K", "T_back_K", "T_1mm_K")]
            print(f"{fds.case_name(q, L):10s} {v:16s} {burnout(a):7.0f} / {burnout(b):7.0f}      {row[0]:8.3f}     "
                  f"{row[1]:8.3f}    {row[2]:7.2f}    {row[3]:7.2f}    {row[4]:7.2f}")


def cmd_pyrolysis(sets):
    print("set        case       burnout PJ/FDS (s)  mass max (g)  mlr max/rms (g/m2s)  T_top max/rms   T_back max/rms  T_1mm max/rms (K)")
    for d in fds.load_decks(sets):
        for q, L, _ in fds.gasification_cases():
            name = f"gasification_{fds.case_name(q, L)}.csv"
            fj, ff = PYRO / d["id"] / "base" / name, fds.RESULTS / d["id"] / "base" / name
            if not (fj.exists() and ff.exists()):
                continue
            j, f = read_series(fj), read_series(ff)
            tend = min(j["time_s"][-1], f["time_s"][-1])
            lim = min(burnout(j), burnout(f), tend) - 10
            m = diff(j, f, "mass_g", tend)
            row = [diff(j, f, k, lim) for k in ("mlr_g_m2s", "T_top_K", "T_back_K", "T_1mm_K")]
            print(f"{d['id']:10s} {fds.case_name(q, L):10s} {burnout(j):7.0f} / {burnout(f):7.0f}   {m[0]:8.3f}     "
                  + "   ".join(f"{a:6.2f} / {b:5.2f}" for a, b in row))


def cmd_runs():
    tot = {}
    for f in sorted(fds.RESULTS.glob("*/runs_*.csv")):
        for r in csv.DictReader(f.open()):
            v = r["variant"]
            t = tot.setdefault(v, [0, 0, 0.0, 0.0, 0.0])
            t[0] += 1
            t[1] += r["ok"] == "True"
            t[2] += float(r["wall_s"])
            t[3] += float(r["cpu_s"]) if r["cpu_s"] not in ("", "nan") else 0.0
            t[4] = max(t[4], float(r["wall_s"]))
    print("variant          runs  ok   sum wall (s)  sum CPU (s)  longest wall (s)")
    for v, (n, ok, w, c, mx) in sorted(tot.items()):
        print(f"{v:16s} {n:4d} {ok:4d}   {w:11.0f}  {c:11.0f}  {mx:10.0f}")


def cmd_stale():
    """Kept inputs that differ from what the current decks and fds.py generate."""
    decks = {d["id"]: d for d in fds.load_decks(include_all=True)}
    cases = {fds.case_name(q, L): (q, L, t) for q, L, t in fds.gasification_cases()}
    n = stale = 0
    for f in sorted(fds.RESULTS.glob("*/*/run_*/gasification.fds")):
        set_id, vdir, run = f.parts[-4], f.parts[-3], f.parts[-2]
        q, L, t_end = cases[run[4:]]
        text = f.read_text()
        cell = float(next(x.split("=")[1] for x in text.replace(",", "\n").splitlines() if "CELL_SIZE(1)" in x))
        dt = float(next(x.split("=")[1] for x in text.replace(",", "\n").splitlines() if x.strip().startswith("DT=")))
        n += 1
        if fds.fds_input(decks[set_id], q, L, t_end, cell=cell, dt=dt) != text:
            stale += 1
            print("STALE", f)
    print(f"{n} inputs checked, {stale} differ from the current decks/generator")


def main(argv):
    if not argv:
        print(__doc__)
        return 1
    cmd, rest = argv[0], argv[1:]
    if cmd == "order":
        cmd_order(rest)
    elif cmd == "nist":
        cmd_nist()
    elif cmd == "submission":
        cmd_nist(rest)
    elif cmd == "refine":
        cmd_refine(rest[0], rest[1:])
    elif cmd == "pyrolysis":
        cmd_pyrolysis(rest)
    elif cmd == "runs":
        cmd_runs()
    elif cmd == "stale":
        cmd_stale()
    else:
        print(__doc__)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
