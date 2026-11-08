#!/usr/bin/env python3
"""FDS 6.11 on the MaCFP-2 PMMA 1D gasification targets.

Reads the canonical decks (../decks/*.toml), writes one FDS SOLID_PHASE_ONLY
input per set x case, runs FDS, and extracts the MaCFP quantities in the
common CSV format of results/<code>/<set>/<variant>/gasification_<q>kW_<L>mm.csv.
Standard library only (Python >= 3.11). See README_fds.md.

    python3 fds.py generate [SET ...] [--cell=5e-5] [--dt=0.01]
    python3 fds.py run      [SET ...] [--cases=25x6,65x12] [--workers=6] [...]
    python3 fds.py extract  [SET ...]

`run` generates, runs and extracts. Sets default to the nine MaCFP-2
gasification sets; `--all` adds the 2023 decks.
"""

import argparse
import csv
import hashlib
import math
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import tomllib
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
SUITE = HERE.parent
DECKS = SUITE / "decks"
RESULTS = SUITE / "results" / "fds"
FDS_EXE = Path(os.environ.get("FDS_EXE", Path.home() / "FDS" / "FDS6" / "bin" / "fds"))

T_AMBIENT = 293.0          # K, initial and far-field temperature (MaCFP-2 v1.3)
SAMPLE_AREA = 0.01         # m², 10 cm x 10 cm (MaCFP-2 output scaling)
TMPM = 273.15              # FDS °C -> K offset
SIGMA_FDS = 5.670373e-8    # FDS 6.11 cons.f90:372
R_FDS = 8.314472           # FDS 6.11 cons.f90:371 (8314.472 J/kmol/K)
R_JL = 8.31446261815324    # Pyrolysis.jl src/Materials/constants.jl

# Property tables: FDS stores k(T) and cp(T) of a MATL at integer kelvin
# (read.f90:7587-7604, ML%K_S(0:I_MAX_TEMP)) and looks them up at INT(T)
# (wall.f90:2820). Each RAMP below has a breakpoint at every integer kelvin of
# [RAMP_T_LO, RAMP_T_HI] and NUMBER_INTERPOLATION_POINTS equal to the span, so
# FDS's own resampling (func.f90:844, nearest node) lands on those breakpoints.
# The entry at I holds the law at I + 0.5 K, the centre of the bin [I, I+1)
# that FDS maps to it.
RAMP_T_LO = 250
RAMP_T_HI = 1500

# Solid grid and time step of the base runs (see README_fds.md, convergence).
CELL_SIZE = 5e-5           # m
DT = 0.01                  # s
DT_DEVC = 1.0              # s
# The results directory of the deck as written; refinement runs append a tag.
VARIANT = "base"

# Remove after a run: everything but the input, the device CSV and the .out
# (trimmed to its header and timing summary).
KEEP_SUFFIXES = (".fds", "_devc.csv", ".out")

# Depths of the 1 mm samples, in cells from 1 mm.
SAMPLE_OFFSETS = (-1.5, -0.5, 0.5, 1.5)

HEADER = ["time_s", "mass_g", "mlr_g_m2s", "T_top_K", "T_back_K", "T_1mm_K", "thickness_mm"]


# ---------------------------------------------------------------------------
# Decks and cases
# ---------------------------------------------------------------------------

def load_deck(path):
    with open(path, "rb") as f:
        return tomllib.load(f)


def load_decks(wanted=(), include_all=False):
    """Gasification decks: by default the MaCFP-2 sets (those with submitted
    predictions); `include_all` adds the 2023 sets; names select explicitly."""
    decks = []
    for p in sorted(DECKS.glob("*.toml")):
        d = load_deck(p)
        d["_path"] = p
        decks.append(d)
    decks = [d for d in decks if not d.get("kinetics_only", False)]
    if not wanted and not include_all:
        decks = [d for d in decks if "predictions" in d]
    if wanted:
        known = {d["id"] for d in decks}
        missing = [w for w in wanted if w not in known]
        if missing:
            raise SystemExit(f"unknown or kinetics-only set(s): {missing}")
        decks = [d for d in decks if d["id"] in wanted]
    return decks


def gasification_cases():
    """(q kW/m², L mm, t_end s) from MacfpPMMA.jl, the one source of truth."""
    text = (SUITE / "MacfpPMMA.jl").read_text()
    block = re.search(r"const GASIFICATION_CASES = \[(.*?)\]", text, re.S).group(1)
    cases = [(int(q), int(L), float(t)) for q, L, t in
             re.findall(r"\(\s*(\d+)\s*,\s*(\d+)\s*,\s*([\d.eE+]+)\s*\)", block)]
    if len(cases) != 6:
        raise RuntimeError(f"expected 6 gasification cases in MacfpPMMA.jl, found {cases}")
    return cases


def case_name(q, L):
    return f"{q}kW_{L}mm"


# ---------------------------------------------------------------------------
# Property forms (same reading as MacfpPMMA.law / MacfpPMMA.positive)
# ---------------------------------------------------------------------------

def law(spec):
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
        return lambda T: _table(Ts, vs, T)
    if form == "switch":
        Tsw = float(spec["T"])
        lo, hi = law(spec["below"]), law(spec["above"])
        return lambda T: lo(T) if T < Tsw else hi(T)
    raise ValueError(f"unknown property form {form}")


def _table(Ts, vs, T):
    """Piecewise linear, held constant outside (the FDS RAMP convention)."""
    if T <= Ts[0]:
        return vs[0]
    if T >= Ts[-1]:
        return vs[-1]
    i = max(j for j in range(len(Ts)) if Ts[j] <= T)
    w = (T - Ts[i]) / (Ts[i + 1] - Ts[i])
    return (1 - w) * vs[i] + w * vs[i + 1]


def positive(f):
    """Floor a law that turns non-positive in [200, 3000] K at 1 % of its
    300 K value (MacfpPMMA.positive; UMD's k = 0.34 - 4.2e-4 T reaches zero at
    810 K). FDS would replace a non-positive k by 10000 W/m/K (wall.f90:2831)."""
    floor_value = 0.01 * abs(f(300.0))
    if any(f(float(T)) <= 0 for T in range(200, 3001)):
        return (lambda T: max(f(T), floor_value)), True
    return f, False


def is_constant(spec):
    return spec["form"] == "constant"


# ---------------------------------------------------------------------------
# FDS input
# ---------------------------------------------------------------------------

def num(x):
    return f"{x:.12g}"


def initial_density(deck):
    """Bulk density of the virgin mixture, 1/sum(Y_i/rho_i) (volume additive,
    as MacfpPMMA._initial_condensed_density and FDS's layer density)."""
    return 1.0 / sum(c["y0"] / c["rho"] for c in deck["components"] if c["y0"] > 0)


def reference_concentration(deck, comp_id):
    """xi_ref of a reactant in Pyrolysis.jl's normalised rate law, its lineage
    inventory (docs/src/guide/user_guide/04_defining-reactions.md): a virgin
    component's initial concentration Y0 * rho0; an intermediate's, what its
    parents can yield (residue yield x the parent's inventory)."""
    rho0 = initial_density(deck)
    comps = {c["id"]: c for c in deck["components"]}
    parents = {}
    for r in deck["reactions"]:
        if r.get("residue"):
            parents.setdefault(r["residue"], []).append((r["reactant"], float(r["residue_yield"])))

    def inventory(cid, seen=()):
        y = comps[cid]["y0"]
        if y > 0:
            return y * rho0
        return sum(nu * inventory(p, seen + (cid,)) for p, nu in parents.get(cid, []) if p not in seen)

    return inventory(comp_id)


def fds_A(deck, r):
    """FDS >= 6.7.8 rate per unit volume: A_F rho_s^n exp(-E/RT) (wall.f90:3463).
    The deck's normalised law r = A xi_ref (xi/xi_ref)^n exp(-E/RT) is the same
    at a material point with A_F = A / xi_ref^(n-1)."""
    n = float(r["n"])
    A = float(r["A"])
    if n == 1.0:
        return A
    basis = deck.get("order_basis", "own")
    if basis != "own":
        raise NotImplementedError(f"{deck['id']}: order_basis {basis}")
    return A / reference_concentration(deck, r["reactant"]) ** (n - 1.0)


def fds_E(E):
    """Preserve E/R of Pyrolysis.jl under FDS's gas constant."""
    return E * R_FDS / R_JL


class RampBook:
    """One RAMP per distinct T-dependent law (identical laws share one)."""

    def __init__(self):
        self.ramps = {}      # key -> (id, lines)
        self.notes = []

    def get(self, kind, spec, scale, label):
        key = (kind, repr(spec))
        if key in self.ramps:
            return self.ramps[key][0]
        rid = f"{kind}_{sum(1 for k in self.ramps if k[0] == kind) + 1}"
        f, floored = positive(law(spec))
        if floored:
            self.notes.append(f"{label} floored at 1 % of its 300 K value above its zero crossing")
        span = RAMP_T_HI - RAMP_T_LO
        lines = []
        for I in range(RAMP_T_LO, RAMP_T_HI + 1):
            Tc = I - TMPM
            lines.append(f"&RAMP ID='{rid}', T={Tc:.2f}, F={num(scale * f(I + 0.5))}, "
                         f"NUMBER_INTERPOLATION_POINTS={span} /")
        self.ramps[key] = (rid, lines)
        return rid

    def lines(self):
        out = []
        for rid, lines in self.ramps.values():
            out += lines
        return out


def fds_input(deck, q, L_mm, t_end, *, cell=CELL_SIZE, dt=DT, dt_devc=DT_DEVC, chid="gasification"):
    """FDS 6.11 input for one set x case, the deck as written."""
    L = L_mm * 1e-3
    comps = deck["components"]
    in_depth = deck["absorption"] == "in-depth"
    book = RampBook()
    out = [
        f"&HEAD CHID='{chid}', TITLE='MaCFP-2 PMMA {deck['id']} {q} kW/m2 {L_mm} mm ({VARIANT})' /",
        "",
        "! One coarse gas cell column: SOLID_PHASE_ONLY solves only the 1D solid.",
        "&MESH IJK=3,3,3, XB=-0.15,0.15,-0.15,0.15,0.0,0.3 /",
        f"&TIME T_END={num(t_end)}, DT={num(dt)}, WALL_INCREMENT=1, LOCK_TIME_STEP=T /",
        f"&MISC SOLID_PHASE_ONLY=T, TMPA={T_AMBIENT - TMPM:.2f} /",
        "&RADI RADIATION=T /",
        f"&DUMP DT_DEVC={num(dt_devc)}, DT_HRR={num(t_end)}, DT_CPU={num(t_end)}, DT_RESTART=1E9, "
        "SMOKE3D=F, MASS_FILE=F, WRITE_XYZ=F /",
        "",
        "! Surrogate pyrolyzate. FDS needs a tracked species for the gas yield; with",
        "! ADJUST_H=F and SOLID_PHASE_ONLY its properties do not enter the solid.",
        "&SPEC ID='MMA', FORMULA='C5H8O2' /",
        "",
    ]
    for c in comps:
        props = [f"ID='{c['id']}'", f"DENSITY={num(c['rho'])}", f"EMISSIVITY={num(c['emissivity'])}"]
        if is_constant(c["k"]):
            props.append(f"CONDUCTIVITY={num(float(c['k']['value']))}")
        else:
            props.append(f"CONDUCTIVITY_RAMP='{book.get('k', c['k'], 1.0, c['id'] + '.k')}'")
        if is_constant(c["cp"]):
            props.append(f"SPECIFIC_HEAT={num(float(c['cp']['value']) / 1000)}")
        else:
            props.append(f"SPECIFIC_HEAT_RAMP='{book.get('cp', c['cp'], 1e-3, c['id'] + '.cp')}'")
        if in_depth:
            props.append(f"ABSORPTION_COEFFICIENT={num(float(c['kappa']))}")
        props.append("ADJUST_H=F")
        reacs = [r for r in deck["reactions"] if r["reactant"] == c["id"]]
        if reacs:
            props.append(f"N_REACTIONS={len(reacs)}")
        for j, r in enumerate(reacs, start=1):
            y = float(r.get("residue_yield", 0.0)) if r.get("residue") else 0.0
            props += [
                f"A({j})={num(fds_A(deck, r))}",
                f"E({j})={num(fds_E(float(r['E'])))}",
                f"N_S({j})={num(float(r['n']))}",
                f"HEAT_OF_REACTION({j})={num(float(r['dH']) / 1000)}",
            ]
            if y > 0:
                props += [f"MATL_ID(1,{j})='{r['residue']}'", f"NU_MATL(1,{j})={num(y)}"]
            props += [f"SPEC_ID(1,{j})='MMA'", f"NU_SPEC(1,{j})={num(1.0 - y)}"]
        out.append("&MATL " + ",\n      ".join(props) + " /")
    out.append("")
    virgin = [c for c in comps if c["y0"] > 0]
    nv = len(virgin)
    surf = [
        "ID='SAMPLE'",
        f"MATL_ID(1,1:{nv})=" + ",".join(f"'{c['id']}'" for c in virgin),
        f"MATL_MASS_FRACTION(1,1:{nv})=" + ",".join(num(c["y0"]) for c in virgin),
        f"THICKNESS(1)={num(L)}",
        f"CELL_SIZE(1)={num(cell)}",
        "STRETCH_FACTOR(1)=1",
        f"N_LAYER_CELLS_MAX(1)={max(999, int(math.ceil(L / cell)) + 10)}",
        f"TMP_INNER={T_AMBIENT - TMPM:.2f}",
        # Absorbed: EMISSIVITY*EXTERNAL_FLUX (radi.f90:4714, SKIP_INRAD drops the
        # gas-phase incoming radiation). Adding sigma*T_inf^4 makes the net
        # radiative flux eps*q + eps*sigma*(T_inf^4 - T_s^4): incident q absorbed
        # at the surface emissivity, re-radiation to 293 K.
        f"EXTERNAL_FLUX={num((q * 1e3 + SIGMA_FDS * T_AMBIENT ** 4) / 1e3)}",
        "SKIP_INRAD=T",
        "HEAT_TRANSFER_COEFFICIENT=0",
        "BACKING='INSULATED'",
    ]
    out.append("&SURF " + ",\n      ".join(surf) + " /")
    out.append("&VENT XB=-0.05,0.05,-0.05,0.05,0.0,0.0, SURF_ID='SAMPLE' /")
    out.append("")
    dev = "XYZ=0,0,0.0, IOR=3, TIME_AVERAGED=F"
    out += [
        f"&DEVC ID='mass', QUANTITY='SURFACE DENSITY', {dev} /",
        f"&DEVC ID='T_top', QUANTITY='WALL TEMPERATURE', {dev} /",
        f"&DEVC ID='T_back', QUANTITY='BACK WALL TEMPERATURE', {dev} /",
        f"&DEVC ID='thickness', QUANTITY='WALL THICKNESS', {dev} /",
        f"&DEVC ID='mass_flux', QUANTITY='MASS FLUX', SPEC_ID='MMA', {dev} /",
    ]
    # Temperature 1 mm below the moving surface (FDS DEPTH is measured from the
    # current front surface, dump.f90:8931-8944): the cells around it, sampled
    # at their initial centres (a depth on a face picks either side), with
    # their actual centre depths; the extractor interpolates to 1 mm.
    for k, d in enumerate(SAMPLE_OFFSETS, start=1):
        d = 1e-3 + d * cell
        out += [
            f"&DEVC ID='T_d{k}', QUANTITY='INSIDE WALL TEMPERATURE', DEPTH={num(d)}, {dev} /",
            f"&DEVC ID='z_d{k}', QUANTITY='INSIDE WALL DEPTH', DEPTH={num(d)}, {dev} /",
        ]
    out.append("")
    if book.ramps:
        out.append(f"! Property RAMPs: integer-kelvin breakpoints {RAMP_T_LO}-{RAMP_T_HI} K (T in degC), "
                   "each holding the law at I + 0.5 K.")
        out += book.lines()
        out.append("")
    out.append("&TAIL /")
    for note in book.notes:
        out.insert(2, f"! NOTE: {note}")
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------------------
# Running
# ---------------------------------------------------------------------------

def fds_env():
    env = {**os.environ, "OMP_NUM_THREADS": "1", "I_MPI_FABRICS": "shm"}
    mpi = FDS_EXE.parent / "intelmpi"
    if mpi.is_dir():
        env["LD_LIBRARY_PATH"] = str(mpi / "lib") + ":" + env.get("LD_LIBRARY_PATH", "")
        env["FI_PROVIDER_PATH"] = str(mpi / "prov")
    return env


def run_dir(set_id, variant, q, L, tag=""):
    return RESULTS / set_id / (variant + tag) / f"run_{case_name(q, L)}"


def parse_out(path):
    """FDS's own timing from the .out: time-stepping and total wall clock, s."""
    text = path.read_text(errors="replace")
    res = {}
    for key, pat in (("fds_time_stepping_s", r"Time Stepping Wall Clock Time \(s\):\s*([\d.Ee+-]+)"),
                     ("fds_total_s", r"Total Elapsed Wall Clock Time \(s\):\s*([\d.Ee+-]+)")):
        m = re.search(pat, text)
        res[key] = float(m.group(1)) if m else float("nan")
    m = re.search(r"Revision\s*:\s*(\S+)", text)
    res["revision"] = m.group(1) if m else "unknown"
    return res


def parse_cpu(path):
    """Total CPU seconds of the single MPI process from <chid>_cpu.csv."""
    if not path.exists():
        return float("nan")
    rows = list(csv.reader(path.open()))
    if len(rows) < 2:
        return float("nan")
    head = [h.strip() for h in rows[0]]
    for name in ("Total T_USED (s)", "Total T_USED"):
        if name in head:
            return float(rows[-1][head.index(name)])
    return float("nan")


def trim_out(path):
    """Drop the per-step diagnostics from an FDS .out (several MB for a long
    run), keeping the header (version, solid grid, devices) and the last step
    with the timing summary."""
    lines = path.read_text(errors="replace").splitlines(keepends=True)
    steps = [i for i, l in enumerate(lines) if re.match(r"\s+Time Step\s+\d+", l)]
    if len(steps) < 2:
        return
    path.write_text("".join(lines[:steps[0]]) + "       [per-step diagnostics removed by fds.py]\n\n"
                    + "".join(lines[steps[-1]:]))


def run_one(deck, q, L, t_end, *, variant=VARIANT, cell=CELL_SIZE, dt=DT, tag="", cpus=None,
            timeout=6 * 3600):
    path = run_dir(deck["id"], variant, q, L, tag)
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True)
    chid = "gasification"
    (path / f"{chid}.fds").write_text(fds_input(deck, q, L, t_end, cell=cell, dt=dt,
                                                chid=chid))
    cmd = [str(FDS_EXE), f"{chid}.fds"]
    if cpus and shutil.which("taskset"):
        cmd = ["taskset", "-c", cpus] + cmd
    t0 = time.perf_counter()
    with (path / "run.log").open("w") as log:
        proc = subprocess.run(cmd, cwd=path, env=fds_env(), stdout=log, stderr=subprocess.STDOUT,
                              timeout=timeout)
    wall = time.perf_counter() - t0
    log_text = (path / "run.log").read_text(errors="replace")
    ok = proc.returncode == 0 and "FDS completed successfully" in log_text
    info = parse_out(path / f"{chid}.out") if (path / f"{chid}.out").exists() else {}
    info["cpu_s"] = parse_cpu(path / f"{chid}_cpu.csv")
    info.update(wall_s=wall, returncode=proc.returncode, ok=ok)
    if not ok:
        info["error"] = " ".join(log_text.strip().splitlines()[-3:])[:300]
    else:
        (path / "run.log").unlink()
        trim_out(path / f"{chid}.out")
        for f in path.iterdir():
            if f.is_file() and not f.name.endswith(KEEP_SUFFIXES):
                f.unlink()
    return path, info


# ---------------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------------

def read_devc(path):
    with path.open() as f:
        rows = list(csv.reader(f))
    names = [h.strip().strip('"') for h in rows[1]]
    data = {n: [] for n in names}
    for row in rows[2:]:
        if not row or not row[0].strip():
            continue
        for n, v in zip(names, row):
            data[n].append(float(v))
    return data


def derivative(t, y):
    """Centred differences, one-sided at the ends (MacfpPMMA.derivative)."""
    n = len(y)
    if n < 2:
        return [math.nan] * n
    d = [0.0] * n
    d[0] = (y[1] - y[0]) / (t[1] - t[0])
    d[-1] = (y[-1] - y[-2]) / (t[-1] - t[-2])
    for i in range(1, n - 1):
        d[i] = (y[i + 1] - y[i - 1]) / (t[i + 1] - t[i - 1])
    return d


def t_at_depth(depth, thickness, samples, T_back):
    """Linear interpolation to `depth` below the moving surface from the
    (cell-centre depth, T) samples. A wall thinner than `depth` gives its back
    temperature (Pyrolysis.jl: the deepest live cell centre)."""
    if thickness <= depth:
        return T_back
    pts = sorted({(z, T) for z, T in samples if 0 < z <= thickness})
    if not pts:
        return math.nan
    below = [p for p in pts if p[0] <= depth]
    above = [p for p in pts if p[0] >= depth]
    if below and above:
        (z0, T0), (z1, T1) = below[-1], above[0]
        if z1 == z0:
            return T0
        return T0 + (T1 - T0) * (depth - z0) / (z1 - z0)
    z, T = min(pts, key=lambda p: abs(p[0] - depth))
    return T


def extract(path):
    """devc CSV -> rows of the common format."""
    d = read_devc(path / "gasification_devc.csv")
    t = d["Time"]
    sd = d["mass"]                                        # kg/m²
    mlr = [-x * 1e3 for x in derivative(t, sd)]           # g/(m² s)
    rows = []
    for i in range(len(t)):
        thick = d["thickness"][i]
        burnt = thick <= 0.0
        T_top = math.nan if burnt else d["T_top"][i] + TMPM
        T_back = math.nan if burnt else d["T_back"][i] + TMPM
        samples = [(d[f"z_d{k}"][i], d[f"T_d{k}"][i] + TMPM) for k in range(1, len(SAMPLE_OFFSETS) + 1)]
        T1 = math.nan if burnt else t_at_depth(1e-3, thick, samples, T_back)
        rows.append([t[i], sd[i] * SAMPLE_AREA * 1e3, mlr[i], T_top, T_back, T1, thick * 1e3])
    return rows, d


def write_series(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write(",".join(HEADER) + "\n")
        for r in rows:
            f.write(",".join("NaN" if (isinstance(x, float) and math.isnan(x)) else f"{x:.9g}" for x in r)
                    + "\n")
    return path


def extract_case(set_id, variant, q, L, tag=""):
    rd = run_dir(set_id, variant, q, L, tag)
    rows, _ = extract(rd)
    return write_series(rd.parent / f"gasification_{case_name(q, L)}.csv", rows)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_cases(spec):
    allc = gasification_cases()
    if not spec:
        return allc
    want = {tuple(int(x) for x in s.split("x")) for s in spec.split(",")}
    return [c for c in allc if (c[0], c[1]) in want]


RUNLOG_FIELDS = ["set", "variant", "case", "cell_m", "dt_s", "ok", "wall_s", "fds_time_stepping_s",
                 "fds_total_s", "cpu_s", "revision", "deck_sha256", "error"]


_lock = threading.Lock()


def deck_hash(deck):
    return hashlib.sha256(deck["_path"].read_bytes()).hexdigest()[:16]


def append_runlog(path, row):
    new = not path.exists()
    with _lock, path.open("a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=RUNLOG_FIELDS, extrasaction="ignore")
        if new:
            w.writeheader()
        w.writerow(row)


def main(argv):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("command", choices=["generate", "run", "extract"])
    p.add_argument("sets", nargs="*")
    p.add_argument("--cases", default="", help="e.g. 25x6,65x12 (default: all six)")
    p.add_argument("--cell", type=float, default=CELL_SIZE)
    p.add_argument("--dt", type=float, default=DT)
    p.add_argument("--tag", default="", help="suffix of the variant directory (refinement studies)")
    p.add_argument("--workers", type=int, default=6)
    p.add_argument("--cpus", default="0-19", help="taskset CPU list ('' to disable)")
    p.add_argument("--out", default="", help="generate: write inputs here instead of the run dirs")
    p.add_argument("--all", action="store_true", help="include the 2023 (MaCFP-3) gasification decks")
    a = p.parse_args(argv)
    decks = load_decks(a.sets, include_all=a.all)
    cases = parse_cases(a.cases)

    if a.command == "generate":
        for d in decks:
            for q, L, t_end in cases:
                dest = Path(a.out) / d["id"] if a.out else run_dir(d["id"], VARIANT, q, L, a.tag)
                dest.mkdir(parents=True, exist_ok=True)
                f = dest / (f"{case_name(q, L)}.fds" if a.out else "gasification.fds")
                f.write_text(fds_input(d, q, L, t_end, cell=a.cell, dt=a.dt))
                print(f)
        return 0

    if a.command == "extract":
        for d in decks:
            for q, L, _ in cases:
                if (run_dir(d["id"], VARIANT, q, L, a.tag) / "gasification_devc.csv").exists():
                    print(extract_case(d["id"], VARIANT, q, L, a.tag))
        return 0

    jobs = [(d, q, L, t) for d in decks for q, L, t in cases]
    # Longest first, so the pool drains evenly.
    jobs.sort(key=lambda j: -j[3] * j[2])
    print(f"{len(jobs)} FDS runs, {a.workers} at a time, cell {a.cell} m, dt {a.dt} s", flush=True)

    def task(job):
        d, q, L, t_end = job
        path, info = run_one(d, q, L, t_end, cell=a.cell, dt=a.dt, tag=a.tag,
                             cpus=a.cpus or None)
        row = {"set": d["id"], "variant": VARIANT + a.tag, "case": case_name(q, L), "cell_m": a.cell,
               "dt_s": a.dt, "deck_sha256": deck_hash(d), **info}
        if info["ok"]:
            extract_case(d["id"], VARIANT, q, L, a.tag)
        append_runlog(RESULTS / d["id"] / f"runs_{VARIANT}{a.tag}.csv", row)
        return row

    failures = 0
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        futs = {pool.submit(task, j): j for j in jobs}
        for fut in as_completed(futs):
            d, q, L, _ = futs[fut]
            try:
                r = fut.result()
                status = "ok" if r["ok"] else "FAILED " + r.get("error", "")
                failures += 0 if r["ok"] else 1
                print(f"{d['id']:10s} {case_name(q, L):10s} {status:8s} wall {r['wall_s']:8.1f} s "
                      f"cpu {r['cpu_s']:8.1f} s", flush=True)
            except Exception as e:  # noqa: BLE001 - report and continue
                failures += 1
                print(f"{d['id']:10s} {case_name(q, L):10s} ERROR {e!r}", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
