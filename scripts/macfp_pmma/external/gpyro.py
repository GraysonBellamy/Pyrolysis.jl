"""Gpyro V2 on the MaCFP-2 PMMA 1D gasification targets.

Python 3.12+, standard library only. Generates native Gpyro V2 inputs from the
canonical decks, runs them (at most --jobs concurrent processes, one OpenMP
thread each) and extracts the common series. See README_gpyro_thermakin.md
for every mapping and approximation.

    python3 gpyro.py run [--sets NIST,GIDAZE+] [--cases 25_6,65_12] [--dz 5e-5] [--dt 0.02]
                         [--suffix _half] [--jobs 4]
    python3 gpyro.py extract [same selectors]
    python3 gpyro.py deck --sets GIDAZE+ --cases 25_6     (print one input)
"""

from pathlib import Path
import argparse
import gzip
import json
import math
import os
import re
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

GPYRO = C.REPO / "docs" / "refs" / "gpyro_build" / "gpyro"
OUT = C.RESULTS / "gpyro"
VERSION = "bff0126c68f2a35c084ef04490ae327185f281a3"

# Gpyro's reaction rate uses R0 = 8.314E-3 kJ/(mol K) (a single-precision
# literal in REACTION_RATE_T); E is rescaled so E/R matches the deck's E read
# with the CODATA R.
R_GPYRO = 0.00831400015205145
TREF = 300.0  # Gpyro's single reference temperature of every power law
# KAPPA >= 1e6 1/m selects Gpyro's surface deposition of the net flux.
SURFACE_KAPPA = 9e9
# Equal reactant and product densities give no gas in Gpyro (solid fraction
# 1 + (ρB/ρA - 1)χ); the product is made this much lighter and χ scaled up.
EQUAL_DENSITY_SPLIT = 1e-4
# A linear property is used only if a power law reproduces it within this
# relative error over FIT_RANGE (near-exact); otherwise the set is skipped.
FIT_TOLERANCE = 1e-3
FIT_RANGE = (293.0, 900.0)

SETS = ["NIST", "GIDAZE+", "Aalto_I", "BUW-FZJ_A", "BUW-FZJ_B", "UMET_GP", "UMET_TK"]
DUMP = 1.0 - 1e-8  # point and profile dump interval (s), just under 1 s
BURNOUT_FRACTION = 5e-3  # a zero-mass shutdown below this mass fraction is burnout,
BURNOUT_SECONDS = 3.0     # or with less mass left than this many seconds at the last MLR


def num(x):
    return f"{x:.12g}"


def power_fit(f, lo, hi):
    """Minimax-in-log fit f(T) ≈ c0 (T/TREF)^n over [lo, hi]: (c0, n, max rel err)."""
    Ts = [lo + (hi - lo) * i / 4000 for i in range(4001)]
    lf = [math.log(f(T)) for T in Ts]
    lt = [math.log(T / TREF) for T in Ts]

    def spread(n):
        r = [a - n * b for a, b in zip(lf, lt)]
        return (max(r) - min(r)) / 2, (max(r) + min(r)) / 2
    a, b = -5.0, 5.0
    for _ in range(200):  # golden-section on the (convex) spread
        m1, m2 = a + 0.382 * (b - a), a + 0.618 * (b - a)
        if spread(m1)[0] < spread(m2)[0]:
            b = m2
        else:
            a = m1
    n = (a + b) / 2
    err, lc0 = spread(n)
    return math.exp(lc0), n, math.expm1(err)


def gp_law(spec, notes, what):
    """Deck property form -> Gpyro (φ0, n) of φ0 (T/TREF)^n."""
    form = spec["form"]
    if form == "constant":
        return float(spec["value"]), 0.0
    if form == "power":
        if float(spec["Tr"]) != TREF:
            raise ValueError(f"{what}: power-law Tr {spec['Tr']} != Gpyro TREF {TREF}")
        return float(spec["c0"]), float(spec["n"])
    if form == "linear":
        f = C.law(spec)
        c0, n, err = power_fit(f, *FIT_RANGE)
        if err > FIT_TOLERANCE:
            raise ValueError(f"{what}: linear law has no near-exact power law (error {err:.2%})")
        notes.append(f"{what}: {spec['a']} + {spec['b']}·T fitted as {c0:.6g}(T/300)^{n:.6g}, "
                     f"max relative error {err:.3%} over {FIT_RANGE[0]:.0f}-{FIT_RANGE[1]:.0f} K")
        return c0, n
    raise ValueError(f"{what}: Gpyro cannot express property form {form!r}")


def section(name, lines):
    return "&GPYRO_" + name + "\n" + "\n".join(lines) + "\n/\n\n"


def deck_text(deck, q, L_mm, t_end, dz, dt, notes=None):
    notes = [] if notes is None else notes
    L = L_mm * 1e-3
    in_depth = deck["absorption"] == "in-depth"
    comps = deck["components"]
    rho = {c["id"]: float(c["rho"]) for c in comps}
    # Products of a reaction with the reactant's own density (UMET's A -> B).
    for r in deck["reactions"]:
        res = r.get("residue", "")
        if res and float(r.get("residue_yield", 0)) > 0 and rho[res] == rho[r["reactant"]]:
            rho[res] = rho[r["reactant"]] * (1 - EQUAL_DENSITY_SPLIT)
            notes.append(f"{res}: density {rho[res]:.6g} (reactant's ×(1-{EQUAL_DENSITY_SPLIT:g})) "
                         "so the reaction can release gas; χ carries the yield")
    nr = 0.0
    if deck.get("rho_law"):  # the "rhoT" variant: the group's ρ(T) power law
        law_ = deck["rho_law"]
        if law_["form"] != "power" or float(law_["Tr"]) != TREF:
            raise ValueError("Gpyro density law must be a power law at TREF")
        scale = float(law_["c0"]) / C.law(law_)(C.T_AMBIENT)
        rho = {k: v * scale for k, v in rho.items()}
        nr = float(law_["n"])
        notes.append(f"density law {law_['c0']}(T/300)^{nr} (rhoT variant)")
    props = [f"NSSPEC={len(comps)},"]
    for i, c in enumerate(comps, 1):
        c0, nc = gp_law(c["cp"], notes, f"{c['id']}.cp")
        k0, nk = gp_law(c["k"], notes, f"{c['id']}.k")
        kappa = float(c["kappa"]) if in_depth else SURFACE_KAPPA
        props.append(
            f"NAME({i})='{c['id']}', R0({i})={num(rho[c['id']])}, NR({i})={num(nr)}, RS0({i})={num(rho[c['id']] * (1.2 if nr else 1.0))}, "
            f"C0({i})={num(c0)}, NC({i})={num(nc)}, K0Z({i})={num(k0)}, NKZ({i})={num(nk)}, "
            f"EMIS({i})={num(float(c['emissivity']))}, KAPPA({i})={num(kappa)}, "
            f"TMELT({i})=5000, DHMELT({i})=0, SIGMA2MELT({i})=1, GAMMA({i})=0, PERMZ({i})=1D-10,")
    rx = [f"NRXNS={len(deck['reactions'])},"]
    orders = set()
    for k, r in enumerate(deck["reactions"], 1):
        y = float(r.get("residue_yield", 0.0))
        res = r.get("residue", "") if y > 0 else ""
        if res:
            # Solid fraction SF = 1 + (ρB/ρA - 1)χ = y; the product's volume
            # follows from its own density.
            chi = (1 - y) / (1 - rho[res] / rho[r["reactant"]])
        else:
            chi = 1.0
        n = float(r["n"])
        orders.add(n)
        E = float(r["E"]) * 1e-3 * (R_GPYRO * 1e3) / C.R_CODATA
        # Gpyro's ΔH is endothermic-positive; ΔH_sol = ΔH_vol is per kg of reactant.
        rx.append(f"CFROM({k})='{r['reactant']}', CTO({k})='{res or 'gases'}', Z({k})={num(float(r['A']))}, "
                  f"E({k})={num(E)}, ORDER({k})={num(n)}, CHI({k})={num(chi)}, "
                  f"DHS({k})={num(float(r['dH']))}, DHV({k})={num(float(r['dH']))},")
    # Gpyro's own normalisation of an n-th order rate (cumulative formed mass)
    # unless the deck reads A against the initial mass.
    basis = deck.get("order_basis", "own")
    conventional = ".TRUE." if basis == "init" else ".FALSE."
    if orders != {1.0} and basis not in ("own", "init"):
        raise ValueError(f"order basis {basis!r} not expressible in Gpyro")
    general = [
        f"DT0={num(dt)}, TAMB={num(C.T_AMBIENT)}, TREF={num(TREF)}, P0=101325,",
        "GX=0, GY=0, GZ=0, NTDMA_ITERATIONS=100, TMPTOL=1D-6, HTOL=1D-10, YITOL=1D-8, EPS=1D-7,",
        "SOLVE_PRESSURE=.FALSE., SOLVE_GAS_YJ=.FALSE., SOLVE_GAS_ENERGY=.FALSE., THERMAL_EQUILIBRIUM=.TRUE.,",
        f"CONVENTIONAL_RXN_ORDER={conventional}, CONSTANT_DHVOL=.TRUE., SHYI_CORRECTION=.TRUE.,",
        "SOLVER_DEFORMATION_MODE=2, FULL_QSG=.FALSE., BLOWING=.FALSE., DTMIN_KILL=1D-9,",
    ]
    points = [("TEMPERATURE", 0, 0.0), ("TEMPERATURE", 0, L), ("THICKNESS", 0, 0.0),
              ("TOTAL_MASS", 0, 0.0), ("MLR", 0, 0.0)]
    output = [
        f"CASENAME='case', N_POINT_QUANTITIES={len(points)}, N_PROFILE_QUANTITIES=2,",
        f"DTDUMP_POINT={num(DUMP)}, DTDUMP_PROFILE={num(DUMP)},",
        "PROFILE_QUANTITY(1)='TEMPERATURE', PROFILE_QUANTITY_INDEX(1)=0, PROFILE_DIRECTION(1)='z', PROFILE_ISKIP(1)=1,",
        "PROFILE_QUANTITY(2)='DLTZN', PROFILE_QUANTITY_INDEX(2)=0, PROFILE_DIRECTION(2)='z', PROFILE_ISKIP(2)=1,",
    ]
    for i, (qty, ix, z) in enumerate(points, 1):
        output.append(f"POINT_QUANTITY({i})='{qty}', POINT_QUANTITY_INDEX({i})={ix}, POINT_Z({i})={num(z)},")
    ic = [f"NIC=1, TMP_INITIAL(1)={num(C.T_AMBIENT)}, TMPG_INITIAL(1)={num(C.T_AMBIENT)}, P_INITIAL(1)=101325, YJ0(1,1)=1,"]
    ic += [f"YI0(1,{i})={num(float(c['y0']))}," for i, c in enumerate(comps, 1)]
    T0 = num(C.T_AMBIENT)
    bc = ["NSURF_IDX=2,",
          # Front: QE absorbed at the surface emissivity, re-radiation to T∞, no convection.
          f"SURF_IDX(1)=1, QE(1)={num(q * 1e3)}, HC(1)=0, RERADIATION(1)=.TRUE., TINF(1)={T0}, TINFG(1)={T0}, "
          "TFIXED(1)=-1, PRES(1)=-1, MDOTPP(1)=0, HM(1)=0, YJINF(1,1)=1,",
          # Back: insulated and impermeable.
          f"SURF_IDX(2)=2, QE(2)=0, HC(2)=0, RERADIATION(2)=.FALSE., TINF(2)={T0}, TINFG(2)={T0}, "
          "TFIXED(2)=-1, PRES(2)=-1, MDOTPP(2)=0, HM(2)=0, YJINF(2,1)=1,"]
    nodes = round(L / dz) + 1
    geom = [f"NMESH=1, ZDIM(1)={num(L)}, NCELLZ(1)={nodes}, XDIM(1)=1, YDIM(1)=1, NCELLX(1)=1, NCELLY(1)=1, "
            "HALF_CELLS_AT_BC(1)=.TRUE., DEFAULT_IC(1)=1, DEFAULT_SURF_IDX(1,5)=1, DEFAULT_SURF_IDX(1,6)=2,"]
    sections = [
        ("GENERAL", general),
        ("OUTPUT", output),
        ("SPROPS", props),
        ("GPROPS", ["NGSPEC=1, IBG=1, IO2=1, CPG=1000, NAME(1)='pyrolysate', M(1)=100,"]),
        ("RXNS", rx),
        ("HGRXNS", ["NHGRXNS=0,"]),
        ("GYIELDS", [f"GYIELDS(1,{k})=1," for k in range(1, len(deck["reactions"]) + 1)]),
        ("HGYIELDS", []),
        # Dumps skip TSTOP itself: stop two steps past the last wanted output.
        ("CASES", [f"NCASES=1, IMESH(1)=1, TSTOP(1)={num(t_end + 2 * dt)}, ZEROD(1)=.FALSE.,"]),
        ("IC", ic),
        ("ALLBC", bc),
        ("GEOM", geom),
    ]
    head = (f"! MaCFP-2 PMMA {deck['id']} (scripts/macfp_pmma/decks/{deck['id']}.toml): "
            f"{q} kW/m2, {L_mm} mm. Generated by scripts/macfp_pmma/external/gpyro.py.\n")
    return head + "".join(section(*s) for s in sections)


def read_csv_rows(path, skip=0):
    rows = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        for i, line in enumerate(f):
            if i < skip:
                continue
            vals = [x for x in line.strip().split(",") if x.strip()]
            if not vals:
                continue
            try:
                rows.append([float(x) for x in vals])
            except ValueError:
                continue
    return rows


def extract(rdir, t_end=math.inf):
    """Common series from the point and profile dumps.

    Temperatures are °C natively. Node 1 sits on the exposed face and the last
    node on the back face (half cells at both). The depth of node i below the
    current surface follows from the dumped cell sizes DLTZN; a consumed cell
    keeps a residual mass (EPS) and a near-zero size.
    """
    summ = rdir / "case_summary_01_0001.csv"
    pts = read_csv_rows(summ, 1)
    prof_T = sorted(rdir.glob("case_profile_01_*_TEMPERATURE.csv*"))[0]
    prof_D = sorted(rdir.glob("case_profile_01_*_DLTZN.csv*"))[0]
    PT = {round(r[0], 4): r[1:] for r in read_csv_rows(prof_T, 1)}
    PD = {round(r[0], 4): r[1:] for r in read_csv_rows(prof_D, 1)}
    t, mass, Ttop, Tback, T1, thick, aux = [], [], [], [], [], [], []
    m0 = None
    for r in pts:
        ti, T_front, T_back_, h, mtot, mlr = r[:6]
        if ti > t_end + 1e-3:
            # Gpyro's own burnout stop jumps the clock past TSTOP for a last dump
            continue
        key = round(ti, 4)
        if key not in PT or key not in PD:
            continue
        T = [x + 273.15 for x in PT[key]]
        dzs = PD[key]
        # node depths: node 1 at 0; node i at sum of preceding sizes + half its own
        depth = [0.0]
        for i in range(1, len(dzs)):
            if i == len(dzs) - 1:
                depth.append(depth[-1] + (dzs[i - 1] if i == 1 else 0.5 * dzs[i - 1]) + dzs[i])
            else:
                depth.append(depth[-1] + (dzs[i - 1] if i == 1 else 0.5 * dzs[i - 1]) + 0.5 * dzs[i])
        m = mtot * 1e-3  # g/m² -> kg/m²
        m0 = m if m0 is None else m0
        burnt = m < 1e-5 * m0
        t.append(ti)
        mass.append(m)
        thick.append(h)
        Ttop.append(math.nan if burnt else T_front + 273.15)
        Tback.append(math.nan if burnt else T_back_ + 273.15)
        T1.append(math.nan if burnt else C.value_at_depth(depth, T, min(1e-3, depth[-1])))
        aux.append((ti, mlr, depth[-1], sum(dzs)))
    rows = C.series_from(t, mass, Ttop, Tback, T1, thick)
    return rows, aux


def run_one(job):
    deck_name, variant, q, L_mm, t_end, dz, dt = job
    deck = C.load_deck(deck_name)
    if variant.startswith("rhoT"):
        deck["rho_law"] = C.RHO_T_LAWS[deck_name]
    vdir = OUT / deck_name / variant
    case = f"run_{q}kW_{L_mm}mm"
    rdir = vdir / case
    rdir.mkdir(parents=True, exist_ok=True)
    for old in rdir.glob("case_*"):
        old.unlink()
    notes = []
    (rdir / "case.data").write_text(deck_text(deck, q, L_mm, t_end, dz, dt, notes))
    env = {**os.environ, "OMP_NUM_THREADS": "1"}
    rc, wall = C.timed_run([str(GPYRO), "case.data"], rdir, rdir / "run.log", env=env)
    log = (rdir / "run.log").read_text(errors="replace")
    meta = dict(set=deck_name, variant=variant, q=q, L_mm=L_mm, t_end=t_end, dz=dz, dt=dt,
                returncode=rc, wall_s=round(wall, 2), notes=notes,
                executable=str(GPYRO.relative_to(C.REPO)), version=VERSION)
    ok = "End of simulation reached successfully." in log
    if VERSION not in log:
        meta["warning"] = "pinned version hash not found in the log"
    if not ok and "Mass in cell has dropped to zero" in log:
        # The last cell emptied within one step: V2's zero-mass shutdown. It is
        # accepted as burnout when less than BURNOUT_FRACTION of the mass is
        # left or the rest would be gone in BURNOUT_SECONDS at the last MLR.
        pts = read_csv_rows(rdir / "case_summary_01_0001.csv", 1)
        if pts and (pts[-1][4] < BURNOUT_FRACTION * pts[0][4]
                    or pts[-1][4] < BURNOUT_SECONDS * pts[-1][5]):
            ok = True
            meta["zero_mass_shutdown_at_burnout"] = pts[-1][0]
    if not ok:
        meta["error"] = log[-600:]
        (rdir / "meta.json").write_text(json.dumps(meta, indent=1))
        raise RuntimeError(f"{deck_name} {variant} {case}: unsuccessful Gpyro termination: {log[-300:]}")
    rows, aux = extract(rdir, t_end)
    C.write_series(vdir / f"gasification_{q}kW_{L_mm}mm.csv", rows)
    with open(rdir / "aux.csv", "w") as f:
        f.write("time_s,native_MLR_g_m2s,back_node_depth_m,sum_DLTZN_m\n")
        for r in aux:
            f.write(",".join(f"{x:.9g}" for x in r) + "\n")
    for p in rdir.glob("case_profile_*.csv"):
        with open(p, "rb") as fi, gzip.open(str(p) + ".gz", "wb", compresslevel=6) as fo:
            shutil.copyfileobj(fi, fo)
        p.unlink()
    meta["t_last"] = rows[-1][0]
    meta["final_mass_g"] = rows[-1][1]
    meta["stopped_at_burnout"] = rows[-1][0] < t_end - 1.5
    (rdir / "meta.json").write_text(json.dumps(meta, indent=1))
    print(f"done {deck_name:10s} {variant:12s} {case:16s} wall {wall:8.1f} s  t_last {rows[-1][0]:.0f}", flush=True)
    return meta


def jobs_from(args):
    sets = args.sets.split(",") if args.sets else SETS
    wanted = set(args.cases.split(",")) if args.cases else None
    jobs = []
    for s in sets:
        for q, L, t_end in C.gasification_cases():
            if wanted and f"{q}_{L}" not in wanted:
                continue
            for v in args.variants.split(","):
                if v == "rhoT" and s not in C.RHO_T_LAWS:
                    continue
                jobs.append((s, v + args.suffix, q, L, t_end, args.dz, args.dt))
    jobs.sort(key=lambda j: -(j[4] * j[3] / j[6]))
    return jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["run", "extract", "deck"])
    ap.add_argument("--sets", default="")
    ap.add_argument("--cases", default="")
    ap.add_argument("--dz", type=float, default=5e-5)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--suffix", default="")
    ap.add_argument("--variants", default="base")
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()
    jobs = jobs_from(args)
    if args.action == "deck":
        s, v, q, L, t_end, dz, dt = jobs[0]
        notes = []
        d = C.load_deck(s)
        if v.startswith("rhoT"):
            d["rho_law"] = C.RHO_T_LAWS[s]
        print(deck_text(d, q, L, t_end, dz, dt, notes))
        print("\n".join("! " + n for n in notes))
        return
    if args.action == "extract":
        for s, v, q, L, t_end, dz, dt in jobs:
            vdir = OUT / s / v
            rdir = vdir / f"run_{q}kW_{L}mm"
            if not (rdir / "case_summary_01_0001.csv").exists():
                continue
            rows, aux = extract(rdir, t_end)
            C.write_series(vdir / f"gasification_{q}kW_{L}mm.csv", rows)
            with open(rdir / "aux.csv", "w") as f:
                f.write("time_s,native_MLR_g_m2s,back_node_depth_m,sum_DLTZN_m\n")
                for r in aux:
                    f.write(",".join(f"{x:.9g}" for x in r) + "\n")
            mp = rdir / "meta.json"
            if mp.exists():
                meta = json.loads(mp.read_text())
                meta["t_last"], meta["final_mass_g"] = rows[-1][0], rows[-1][1]
                meta["stopped_at_burnout"] = rows[-1][0] < t_end - 1.5
                mp.write_text(json.dumps(meta, indent=1))
        return
    results = C.run_parallel(jobs, run_one, args.jobs)
    log = OUT / "runs.jsonl"
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "a") as f:
        for job, meta, err in results:
            f.write(json.dumps(dict(job=list(job), meta=meta, error=err)) + "\n")


if __name__ == "__main__":
    main()
