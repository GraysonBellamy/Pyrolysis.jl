"""Agreement of extracted series with an institution's own MaCFP-2 submission.

    python3 sanity.py CODE SET VARIANT [CASE ...]      e.g. thermakin UMD rand 25_6
    python3 sanity.py --table                          the README tables

For each case: mass and back-temperature differences (max and RMS over the
submission's time span, burnout tail excluded: submission mass > 1 % of its
initial), time to 50 % mass loss, and the submission's top-temperature column
against our T_top and, for ThermaKin, the candidates in aux.csv.
"""

from pathlib import Path
import csv
import math
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402


def interp(ts, ys, t):
    if not ts or t < ts[0] or t > ts[-1]:
        return math.nan
    lo, hi = 0, len(ts) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ts[mid] <= t:
            lo = mid
        else:
            hi = mid
    if ts[hi] == ts[lo]:
        return ys[lo]
    w = (t - ts[lo]) / (ts[hi] - ts[lo])
    return (1 - w) * ys[lo] + w * ys[hi]


def t_half(ts, ms):
    m0 = ms[0]
    for t, m in zip(ts, ms):
        if m < 0.5 * m0:
            return t
    return math.nan


def stats(diffs):
    d = [x for x in diffs if not math.isnan(x)]
    if not d:
        return math.nan, math.nan
    return max(abs(x) for x in d), math.sqrt(sum(x * x for x in d) / len(d))


def compare(code, set_name, variant, q, L):
    deck = C.load_deck(set_name)
    sp = C.prediction_path(deck, q, L)
    path = C.RESULTS / code / set_name / variant / f"gasification_{q}kW_{L}mm.csv"
    if sp is None or not path.exists():
        return None
    sub = [r for r in C.read_matl_csv(sp) if len(r) >= 4]
    X = C.read_series(path)
    st = [r[0] for r in sub]
    sm, sb, stop = [r[1] for r in sub], [r[2] for r in sub], [r[3] for r in sub]
    m0 = sm[0]
    aux = {}
    ap = path.parent / f"run_{q}kW_{L}mm" / "aux.csv"
    if code == "thermakin" and ap.exists():
        with open(ap) as f:
            rows = list(csv.reader(f))
        names = rows[0]
        cols = list(zip(*[[float(x) for x in r] for r in rows[1:]]))
        aux = {n: list(c) for n, c in zip(names, cols)}
    dm, db, dt_ = [], [], []
    dcand = {k: [] for k in ("T_mean_top_1mm_K", "T_0p5mm_K", "T_fixed_height_L_minus_1mm_K")}
    for t, m, b, tp in zip(st, sm, sb, stop):
        if m < 0.01 * m0 or t > X["time_s"][-1]:
            continue
        dm.append(interp(X["time_s"], X["mass_g"], t) - m)
        db.append(interp(X["time_s"], X["T_back_K"], t) - b)
        dt_.append(interp(X["time_s"], X["T_top_K"], t) - tp)
        for k in dcand:
            if k in aux:
                dcand[k].append(interp(aux["time_s"], aux[k], t) - tp)
    out = dict(case=f"{q}_{L}", m0_sub=m0, m0=X["mass_g"][0],
               t50_sub=t_half(st, sm), t50=t_half(X["time_s"], X["mass_g"]),
               mass=stats(dm), back=stats(db), top=stats(dt_),
               end_sub=st[-1], end=X["time_s"][-1])
    for k, v in dcand.items():
        if v:
            out[k] = stats(v)
    return out


def fmt(r):
    s = (f"{r['case']:>6}  m0 {r['m0_sub']:.2f}/{r['m0']:.2f} g  t50 {r['t50_sub']:.0f}/{r['t50']:.0f} s  "
         f"mass max/rms {r['mass'][0]:.2f}/{r['mass'][1]:.2f} g  T_back {r['back'][0]:.1f}/{r['back'][1]:.1f} K  "
         f"T_top {r['top'][0]:.1f}/{r['top'][1]:.1f} K")
    for k, lab in (("T_mean_top_1mm_K", "mean1mm"), ("T_0p5mm_K", "T0.5mm"), ("T_fixed_height_L_minus_1mm_K", "fixedH")):
        if k in r:
            s += f"  {lab} {r[k][0]:.1f}/{r[k][1]:.1f}"
    return s + f"  end {r['end_sub']:.0f}/{r['end']:.0f} s"


def runs_table():
    """Every run's meta.json: set, variant, case, wall time, end, flags."""
    import json
    tot = {}
    for code in ("gpyro", "thermakin"):
        metas = sorted((C.RESULTS / code).glob("*/*/run_*/meta.json"))
        print(f"== {code}: {len(metas)} runs")
        for m in metas:
            d = json.loads(m.read_text())
            flags = []
            if d.get("diverged"):
                flags.append("diverged")
            if d.get("stopped_at_burnout"):
                flags.append("stopped at burnout")
            if d.get("zero_mass_shutdown_at_burnout") is not None:
                flags.append("zero-mass shutdown")
            if d.get("error"):
                flags.append("ERROR")
            v = d["variant"]
            print(f"  {d['set']:10s} {v:14s} {d['q']:>2d} kW {d['L_mm']:>2d} mm  wall {d['wall_s']:8.1f} s  "
                  f"t_last {d.get('t_last', float('nan')):7.0f}/{d['t_end']:.0f} s  "
                  f"final {d.get('final_mass_g', float('nan')):.3f} g  {', '.join(flags)}")
            tot[(code, v)] = tot.get((code, v), 0.0) + d["wall_s"]
    for k, v in sorted(tot.items()):
        print(f"  total wall {k[0]} {k[1]}: {v:.0f} s")


def main(argv):
    if argv and argv[0] == "--runs":
        runs_table()
        return
    if argv and argv[0] == "--table":
        pairs = [("thermakin", "UMD", v) for v in ("base", "rand")] + \
                [("thermakin", "UMET_TK", v) for v in ("base", "rand")] + \
                [("gpyro", "GIDAZE+", "base"), ("gpyro", "UMET_GP", "base")]
        for code, s, v in pairs:
            print(f"== {code} {s} {v}  (submission/ours; differences ours - submission)")
            for q, L, _ in C.gasification_cases():
                r = compare(code, s, v, q, L)
                if r:
                    print("  " + fmt(r))
        return
    code, s, v = argv[:3]
    cases = argv[3:] or [f"{q}_{L}" for q, L, _ in C.gasification_cases()]
    for c in cases:
        q, L = (int(x) for x in c.split("_"))
        r = compare(code, s, v, q, L)
        print(fmt(r) if r else f"{c}: missing")


if __name__ == "__main__":
    main(sys.argv[1:])
