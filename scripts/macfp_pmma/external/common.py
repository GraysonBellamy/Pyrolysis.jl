"""Shared pieces of the Gpyro and ThermaKin2Ds MaCFP-2 PMMA drivers.

Python 3.12+, standard library only. The decks (../decks/*.toml) and the
scenario (GASIFICATION_CASES in ../MacfpPMMA.jl) are the single source of
truth; nothing here restates a deck value.
"""

from pathlib import Path
import concurrent.futures
import csv
import math
import os
import re
import subprocess
import time
import tomllib

HERE = Path(__file__).resolve().parent
SUITE = HERE.parent
REPO = SUITE.parent.parent
DECKS = SUITE / "decks"
RESULTS = SUITE / "results"

T_AMBIENT = 293.0
SAMPLE_AREA = 0.01  # m², 10 cm × 10 cm (MaCFP-2 output scaling)
R_CODATA = 8.314462618  # J/(mol K): the reference reading of E
HEADER = "time_s,mass_g,mlr_g_m2s,T_top_K,T_back_K,T_1mm_K,thickness_mm"

# The calibrating group's temperature-dependent density, which the canonical
# decks hold at its 293 K value (deck notes). Used only by the extra "rhoT"
# variant in the group's own code, to explain differences from its submission.
RHO_T_LAWS = {
    "UMET_TK": {"form": "linear", "a": 1380.0, "b": -0.6},
    "UMET_GP": {"form": "power", "c0": 1200.0, "Tr": 300.0, "n": -0.236},
}

GASIFICATION_SETS = ["NIST", "GIDAZE+", "Aalto_I", "BUW-FZJ_A", "BUW-FZJ_B",
                     "DBI_1", "UMD", "UMET_GP", "UMET_TK"]


def load_deck(name):
    with open(DECKS / f"{name}.toml", "rb") as f:
        return tomllib.load(f)


def gasification_cases():
    """(q kW/m², L mm, t_end s) from MacfpPMMA.jl's GASIFICATION_CASES."""
    text = (SUITE / "MacfpPMMA.jl").read_text()
    block = re.search(r"const GASIFICATION_CASES = \[(.*?)\]", text, re.S)
    if not block:
        raise RuntimeError("GASIFICATION_CASES not found in MacfpPMMA.jl")
    cases = [(int(q), int(L), float(t)) for q, L, t in
             re.findall(r"\(\s*(\d+)\s*,\s*(\d+)\s*,\s*([\d.eE+]+)\s*\)", block[1])]
    if not cases:
        raise RuntimeError("GASIFICATION_CASES is empty")
    return cases


def law(spec):
    """A deck property form as a function of T [K]."""
    form = spec["form"]
    if form == "constant":
        v = float(spec["value"])
        return lambda T: v
    if form == "linear":
        a, b = float(spec["a"]), float(spec["b"])
        return lambda T: a + b * T
    if form == "power":
        c0, Tr, n = float(spec["c0"]), float(spec["Tr"]), float(spec["n"])
        return lambda T: c0 * (T / Tr) ** n
    if form == "table":
        Ts, vs = [float(x) for x in spec["T"]], [float(x) for x in spec["v"]]

        def table(T):
            if T <= Ts[0]:
                return vs[0]
            if T >= Ts[-1]:
                return vs[-1]
            i = max(j for j in range(len(Ts)) if Ts[j] <= T)
            w = (T - Ts[i]) / (Ts[i + 1] - Ts[i])
            return (1 - w) * vs[i] + w * vs[i + 1]
        return table
    if form == "switch":
        Tsw = float(spec["T"])
        lo, hi = law(spec["below"]), law(spec["above"])
        return lambda T: lo(T) if T < Tsw else hi(T)
    raise ValueError(f"unknown property form {form}")


def initial_condensed_density(deck):
    return 1.0 / sum(c["y0"] / c["rho"] for c in deck["components"] if c["y0"] > 0)


def centred_derivative(t, y):
    n = len(y)
    if n < 2:
        return [math.nan] * n
    d = [0.0] * n
    d[0] = (y[1] - y[0]) / (t[1] - t[0])
    d[-1] = (y[-1] - y[-2]) / (t[-1] - t[-2])
    for i in range(1, n - 1):
        d[i] = (y[i + 1] - y[i - 1]) / (t[i + 1] - t[i - 1])
    return d


def value_at_depth(depths, values, depth):
    """Linear interpolation in depth below the current surface; the outermost
    value for depths shallower than the first point, the innermost beyond the
    last (the object is thinner than the requested depth)."""
    pts = sorted(zip(depths, values))
    if not pts:
        return math.nan
    if depth <= pts[0][0]:
        return pts[0][1]
    if depth >= pts[-1][0]:
        return pts[-1][1]
    for (d0, v0), (d1, v1) in zip(pts, pts[1:]):
        if d0 <= depth <= d1:
            w = (depth - d0) / (d1 - d0) if d1 > d0 else 0.0
            return (1 - w) * v0 + w * v1
    return math.nan


def write_series(path, rows):
    """rows: tuples (time_s, mass_g, mlr_g_m2s, T_top_K, T_back_K, T_1mm_K, thickness_mm)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write(HEADER + "\n")
        for r in rows:
            f.write(",".join("NaN" if (isinstance(x, float) and math.isnan(x)) else f"{x:.9g}"
                             for x in r) + "\n")
    return path


def series_from(t, mass_kg_m2, T_top, T_back, T_1mm, thickness_m):
    """Common post-processing: mass for a 0.01 m² sample in g, MLR in g/(m² s)
    as minus the centred difference of the areal condensed mass."""
    mlr = [-d * 1e3 for d in centred_derivative(t, mass_kg_m2)]
    return [(t[i], mass_kg_m2[i] * SAMPLE_AREA * 1e3, mlr[i], T_top[i], T_back[i], T_1mm[i],
             thickness_m[i] * 1e3) for i in range(len(t))]


def read_series(path):
    with open(path) as f:
        rows = list(csv.reader(f))
    names = rows[0]
    cols = {n: [] for n in names}
    for r in rows[1:]:
        for n, v in zip(names, r):
            cols[n].append(float(v))
    return cols


def matl_db():
    return Path(os.environ.get("MATL_DB", Path.home() / "Documents" / "git" / "matl-db"))


def read_matl_csv(path):
    """Numeric rows of a matl-db prediction file (header rows skipped)."""
    out = []
    with open(path, encoding="utf-8-sig", errors="replace") as f:
        for line in f:
            fields = [x.strip().strip('"') for x in re.split(r"[,;\t]", line.strip())]
            if not fields or not fields[0]:
                continue
            try:
                first = float(fields[0])
            except ValueError:
                continue
            row = [first]
            for x in fields[1:]:
                try:
                    row.append(float(x))
                except ValueError:
                    if x:
                        row.append(math.nan)
            out.append(row)
    return out


def prediction_path(deck, q, L):
    d = matl_db() / "PMMA" / "Computational_Results" / "Gasification_Predictions"
    for kw in ("kW", "kw"):
        p = d / f"{deck['predictions']}_Gasification_{q}{kw}_{L}mm.csv"
        if p.is_file():
            return p
    return None


def run_parallel(jobs, worker, max_workers=4):
    """Run worker(job) for every job with at most max_workers processes alive."""
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {pool.submit(worker, j): j for j in jobs}
        for fut in concurrent.futures.as_completed(futures):
            job = futures[fut]
            try:
                results.append((job, fut.result(), None))
            except Exception as err:  # recorded, the campaign continues
                results.append((job, None, repr(err)[:500]))
                print(f"FAILED {job}: {err!r}"[:600], flush=True)
    return results


def timed_run(cmd, cwd, log_path, stdin_text=None, env=None, timeout=None):
    t0 = time.monotonic()
    with open(log_path, "w") as log:
        proc = subprocess.run(cmd, cwd=cwd, input=stdin_text, stdout=log, stderr=subprocess.STDOUT,
                              text=True, env=env, timeout=timeout)
    return proc.returncode, time.monotonic() - t0
