"""ThermaKin2Ds on the MaCFP-2 PMMA 1D gasification targets.

Python 3.12+, standard library only. Generates native components (.cmp) and
conditions (.cnd) files from the canonical decks, runs them (at most --jobs
concurrent processes) and extracts the common series. See
README_gpyro_thermakin.md for every mapping and approximation.

    python3 thermakin.py run [--sets UMD,NIST] [--variants base,rand] [--cases 25_6,65_12]
                             [--dx 5e-5] [--dt 0.01] [--suffix _half] [--seed 1] [--jobs 4]
    python3 thermakin.py extract [same selectors]

Variants: base = MAX absorption (ThermaKin's default), rand = RAND absorption
with a fixed seed (THERMAKIN_SEED, patched build from build_thermakin_seeded.sh).
--suffix names refinement or seed studies (variant directory <variant><suffix>).
The ThermaKin source and executables stay in ignored docs/refs/ directories.
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

TK = C.REPO / "docs" / "refs" / "thermakin_build" / "thermakin"
TK_SEEDED = C.REPO / "docs" / "refs" / "thermakin_build_seeded" / "thermakin"
OUT = C.RESULTS / "thermakin"

# ThermaKin's own constants (Constants.h). Activation energies are rescaled so
# E/R matches the deck's reading with the CODATA R; σ stays native.
TK_R = 8.31447
TK_SIGMA = 5.6704e-08
# Mass absorption coefficient [m²/kg] for surface-absorbing decks: the top
# element absorbs the whole net flux.
SURFACE_ALPHA = 1e6
# Gas: the pyrolyzate leaves through the top boundary, J = a·ξ_top (LIN a 0),
# a = 0.05 m/s (the reference manual's rapid-removal value); λ for every
# condensed component (deck [gas] D when given).
GAS_LIN_A = 0.05
GAS_LAMBDA = 2e-5
GAS_CP_SMALL = 1.0  # J/(kg K) when the deck has no gas cp (its sensible heat is then negligible)
# Fast glass transition (UMD's pseudo-reaction): rate A·ξ above T_switch.
GLASS_A = 100.0
BACKING = "backing"
BACKING_RHO_FACTOR = 0.1
BACKING_ELEMENTS = 2

SETS = ["NIST", "GIDAZE+", "Aalto_I", "BUW-FZJ_A", "BUW-FZJ_B", "UMD", "UMET_TK"]
IN_DEPTH = {"GIDAZE+", "Aalto_I", "UMD", "UMET_TK"}


def num(x):
    return f"{x:.10g}"


def tk_law(spec, T_ref=None):
    """Deck property form -> ThermaKin (p0, p1, p2, n) of p0 + p1 T + p2 T^n."""
    form = spec["form"]
    if form == "constant":
        return (float(spec["value"]), 0.0, 0.0, 0.0)
    if form == "linear":
        return (float(spec["a"]), float(spec["b"]), 0.0, 0.0)
    if form == "power":
        c0, Tr, n = float(spec["c0"]), float(spec["Tr"]), float(spec["n"])
        return (0.0, 0.0, c0 / Tr**n, n)
    raise ValueError(f"ThermaKin cannot express property form {form!r} exactly")


def split_switch(spec):
    """(T_switch, below_spec, above_spec) or (None, spec, spec)."""
    if spec["form"] == "switch":
        return float(spec["T"]), spec["below"], spec["above"]
    return None, spec, spec


def material(deck, dx):
    """Components (name -> dict), transition reactions, decomposition reactions,
    initial mass fractions of the sample and of the backing element."""
    in_depth = deck["absorption"] == "in-depth"
    lam = float(deck["gas"].get("D", GAS_LAMBDA))
    comps, transitions = {}, []
    names = {}  # deck id -> ThermaKin name of the high-temperature (reacting) form
    initial = {}
    kappa_max = max(float(c.get("kappa", 0.0)) for c in deck["components"])

    rho_law = deck.get("rho_law")

    def add(name, c, cp_spec, k_spec, eps, alpha, transport):
        rho = tk_law(rho_law) if (rho_law and name != BACKING) else (float(c["rho"]), 0, 0, 0)
        comps[name] = dict(state="S", rho=rho, cp=tk_law(cp_spec),
                           k=tk_law(k_spec), transport=(transport, 0, 0, 0), eps=eps, alpha=alpha)

    for c in deck["components"]:
        rho = float(c["rho"])
        alpha = float(c.get("kappa", 0.0)) / rho if in_depth else SURFACE_ALPHA
        Tc, cp_lo, cp_hi = split_switch(c["cp"])
        Tk, k_lo, k_hi = split_switch(c["k"])
        Tsw = Tc if Tc is not None else Tk
        if Tc is not None and Tk is not None and Tc != Tk:
            raise ValueError("switch temperatures of cp and k differ")
        names[c["id"]] = c["id"]
        add(c["id"], c, cp_hi, k_hi, float(c["emissivity"]), alpha, lam)
        if Tsw is not None and c["y0"] > 0:
            if C.T_AMBIENT >= Tsw:
                raise ValueError("initial temperature above the switch")
            glass = c["id"] + "_glass"
            add(glass, c, cp_lo, k_lo, float(c["emissivity"]), alpha, lam)
            transitions.append((glass, c["id"], Tsw))
            initial[glass] = float(c["y0"])
        elif c["y0"] > 0:
            initial[c["id"]] = float(c["y0"])

    # Backing: two non-reacting, non-emitting elements under the sample, a
    # thermal copy of the first virgin component at a tenth of its density
    # (so a tenth of its heat capacity per volume). ThermaKin's bottom boundary
    # always re-radiates from the element its absorption search selects; the
    # bottom backing element is that element (ε = 0), so the back face is
    # adiabatic. Their absorption just wins the bottom search against the
    # sample element above, so from the top they are selected only once the
    # sample is about one element thick. Two elements, because ThermaKin merges
    # a top element thinner than 0.2·Δx into the one below: a last residue
    # layer then mixes into the upper backing element, never the bottom one.
    first = next(c for c in deck["components"] if c["y0"] > 0)
    backing = dict(first, rho=float(first["rho"]) * BACKING_RHO_FACTOR)
    s = 1.0 - math.exp(-kappa_max * dx) if in_depth else 1.0
    if in_depth:
        a_b = 1.05 * s / (1.0 + s)
        backing_alpha = -math.log(1.0 - a_b) / (backing["rho"] * dx)
    else:
        backing_alpha = SURFACE_ALPHA
    Tc, cp_lo, cp_hi = split_switch(first["cp"])
    Tk, k_lo, k_hi = split_switch(first["k"])
    add(BACKING, backing, cp_hi, k_hi, 0.0, backing_alpha, 1e-30)
    backing_initial = {BACKING: 1.0}
    if (Tc or Tk) is not None:
        add(BACKING + "_glass", backing, cp_lo, k_lo, 0.0, backing_alpha, 1e-30)
        transitions.append((BACKING + "_glass", BACKING, Tc or Tk))
        backing_initial = {BACKING + "_glass": 1.0}

    gas_cp = float(deck["gas"].get("cp", 0.0)) or GAS_CP_SMALL
    rho_gas = C.initial_condensed_density(deck)
    k_gas = C.law(first["k"])(C.T_AMBIENT)
    comps["gas"] = dict(state="G", rho=(rho_gas, 0, 0, 0), cp=(gas_cp, 0, 0, 0), k=(k_gas, 0, 0, 0),
                        transport=(lam, 0, 0, 0), eps=0.0, alpha=0.0)

    reactions = []
    for r in deck["reactions"]:
        if float(r["n"]) != 1.0:
            raise ValueError(f"{deck['id']}: ThermaKin rates are first order (n = {r['n']})")
        y = float(r.get("residue_yield", 0.0))
        res = r.get("residue", "")
        E = float(r["E"]) * TK_R / C.R_CODATA
        # ThermaKin's positive heat is exothermic; the deck's dH is endothermic-positive,
        # both per kg of reactant consumed (stoichiometry 1 for the reactant).
        if res and y > 0:
            reactions.append((names[r["reactant"]], names[res], "gas", (1, 0, y, 1 - y),
                              float(r["A"]), E, -float(r["dH"])))
        else:
            reactions.append((names[r["reactant"]], "gas", "NOCOMP", (1, 0, 1, 0),
                              float(r["A"]), E, -float(r["dH"])))
    return comps, transitions, reactions, initial, backing_initial


def cmp_text(deck, dx):
    comps, transitions, reactions, _, _ = material(deck, dx)
    law4 = lambda v: " ".join(num(x) for x in v)
    out = [f"//MaCFP-2 PMMA property set {deck['id']} (scripts/macfp_pmma/decks/{deck['id']}.toml)\n\n"]
    for name, c in comps.items():
        out.append(f"COMPONENT: {name}\nSTATE: {c['state']}\nDENSITY: {law4(c['rho'])}\n"
                   f"HEAT CAPACITY: {law4(c['cp'])}\nCONDUCTIVITY: {law4(c['k'])}\n"
                   f"TRANSPORT: {law4(c['transport'])}\n"
                   f"EMISSIVITY & ABSORPTION: {num(c['eps'])} {num(c['alpha'])}\n\n")
    # Parallel (volume-weighted) conductivity and gas transfer, as FDS, Gpyro
    # and Pyrolysis.jl mix k; gases occupy no volume.
    out.append("MIXTURES\nS SWELLING: 0\nL SWELLING: 0\nG SWELLING LIMIT: 1e-30\n"
               "PARALL CONDUCTIVITY: 1\nPARALL TRANSPORT: 1\n\n")
    for frm, to, Tsw in transitions:
        out.append(f"REACTION: {frm} + NOCOMP -> {to} + NOCOMP\nSTOICHIOMETRY: 1 0 1 0\n"
                   f"ARRHENIUS: {num(GLASS_A)} 0\nHEAT: 0 0 0 0\nTEMP LIMIT: L {num(Tsw)}\n\n")
    for r1, p1, p2, st, A, E, h in reactions:
        out.append(f"REACTION: {r1} + NOCOMP -> {p1} + {p2}\n"
                   f"STOICHIOMETRY: {' '.join(num(s) for s in st)}\n"
                   f"ARRHENIUS: {num(A)} {num(E)}\nHEAT: {num(h)} 0 0 0\nTEMP LIMIT: L 0\n\n")
    return "".join(out)


def cnd_text(deck, q, L_mm, t_end, dx, dt, mode):
    _, _, _, initial, backing_initial = material(deck, dx)
    L = L_mm * 1e-3
    T0 = C.T_AMBIENT
    layers = [(L, initial), (BACKING_ELEMENTS * dx, backing_initial)]
    structure = "".join(
        f"THICKNESS: {num(t)}\nTEMPERATURE: {num(T0)}\nMASS FRACTIONS:\n"
        + "".join(f"{c} {num(v)}\n" for c, v in fr.items()) + "\n" for t, fr in layers)
    # Absorbed εq and re-radiation εσ(T⁴ - T∞⁴): ThermaKin deposits ε(I - σT⁴)
    # at the absorbing element, so I = q + σT∞⁴.
    flux = q * 1e3 + TK_SIGMA * T0**4
    top = (f"TOP BOUNDARY\n\nMASS TRANSPORT: YES\ngas LIN {num(GAS_LIN_A)} 0\n\n"
           f"OUTSIDE INIT TEMP: {num(T0)}\nOUTSIDE HEAT RATE: 0 0 0 0\nCONVECTION COEFF: 0\n\n"
           f"EXTERNAL RADIATION: YES\nTIME PROG1: {num(flux)} 0 1e30\nTIME PROG2: {num(flux)} 0 0\n"
           f"REPEAT: NO\nABSORPTION MODE: {mode}\n\nFLAME: NO\n\n")
    bottom = (f"BOTTOM BOUNDARY\n\nMASS TRANSPORT: NO\n\nOUTSIDE INIT TEMP: {num(T0)}\n"
              f"OUTSIDE HEAT RATE: 0 0 0 0\nCONVECTION COEFF: 0\n\nEXTERNAL RADIATION: NO\n\nFLAME: NO\n\n")
    steps = 1.0 / dt
    if abs(steps - round(steps)) > 1e-9:
        raise ValueError("output interval of 1 s is not a whole number of steps")
    return (f"//MaCFP-2 PMMA {deck['id']}: {q} kW/m2, {L_mm} mm, {mode} absorption\n\n"
            "OBJECT TYPE: 1D\n\nOBJECT STRUCTURE\n****************\n\n" + structure
            + "OBJECT BOUNDARIES\n*****************\n\n" + top + bottom
            + "INTEGRATION PARAMETERS\n**********************\n\n"
            f"ELEMENT SIZE: {num(dx)}\nTIME STEP: {num(dt)}\n"
            # Half a step short: the loop runs while t < DURATION with t = k·dt.
            f"DURATION: {num(t_end - dt / 2)}\n\n"
            f"OUTPUT FREQUENCY:\nELEMENTS: 1\nTIME STEPS: {round(steps)}\n")


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def snapshots(path, allow_incomplete=False):
    """Stream the output: per block (t, depths, T, {comp: ξ}, thickness, total mass, flows).
    allow_incomplete accepts a run that stopped early (ThermaKin's "integration
    diverged"); its last complete block is the last one yielded."""
    names = flow_names = None
    block = None
    complete = False
    with open_text(path) as f:
        for line in f:
            if line.startswith("Time [s] ="):
                if block is not None:
                    yield block
                block = dict(t=float(line.split("=")[1]), depth=[], T=[], xi=[], flow={})
                names = None
                continue
            if block is None:
                continue
            if line.startswith("BOUNDARY"):
                flow_names = [x for x in line.rstrip("\n").split("\t") if x.strip()][4:]
            elif line.startswith(("TOP", "BOTTOM")):
                fields = [x for x in line.split("\t") if x.strip()]
                vals = [float(x) for x in fields[1:]]
                block["flow"][fields[0].strip()] = dict(zip(flow_names, vals[2:]))
            elif line.startswith("FROM TOP"):
                names = [x.strip() for x in line.rstrip("\n").split("\t") if x.strip()][3:]
                block["names"] = names
            elif line.startswith("Total thickness"):
                block["thickness"] = float(line.split("=")[1])
            elif line.startswith("Total mass"):
                block["mass"] = float(line.split("=")[1])
            elif line.startswith("Calculations are complete"):
                complete = True
            elif names is not None and line[:1] in "0123456789":
                v = [float(x) for x in line.split("\t") if x.strip()]
                block["depth"].append(v[0])
                block["T"].append(v[1])
                block["xi"].append(v[2:])
    if block is not None and "mass" in block:
        yield block
    if not complete and not allow_incomplete:
        raise RuntimeError(f"{path}: ThermaKin did not report completion")


def extract(run_dir, case, L_mm, allow_incomplete=False):
    """Common series plus an auxiliary file (native outflow, UMD-style probes).

    Sample = every condensed component but the backing; its mass is the printed
    total minus the backing (constant: ρ_b·2Δx) and the gas held in the
    elements. Thickness = printed total minus the backing's 2Δx. T_top is the
    top element (its centre is half an element below the surface); T_back is
    interpolated at the sample/backing interface. After burnout (sample below
    1e-5 of its initial mass) no sample is left and the temperatures are NaN.
    """
    out = run_dir / f"{case}.out.gz"
    if not out.exists():
        out = run_dir / f"{case}.out"
    cnd = (run_dir / f"{case}.cnd").read_text()
    cmpt = (run_dir / f"{case}.cmp").read_text()
    dx = float(re.search(r"ELEMENT SIZE: (\S+)", cnd)[1])
    rho_b = float(re.search(r"COMPONENT: " + BACKING + r"\n.*?DENSITY: (\S+)", cmpt, re.S)[1])
    backing_mass = rho_b * BACKING_ELEMENTS * dx
    backing_thick = BACKING_ELEMENTS * dx
    t, mass, Ttop, Tback, T1, thick, aux = [], [], [], [], [], [], []
    L = L_mm * 1e-3
    m0 = None
    for b in snapshots(out, allow_incomplete):
        names = b["names"]
        ig = names.index("gas")
        n = len(b["depth"])
        # widths from the bottom up: centre_i = S_i - w_i/2, S_N = total thickness
        widths = [0.0] * n
        S = b["thickness"]
        for j in range(n - 1, -1, -1):
            widths[j] = 2.0 * (S - b["depth"][j])
            S -= widths[j]
        gas_mass = sum(b["xi"][j][ig] * widths[j] for j in range(n))
        m = b["mass"] - backing_mass - gas_mass
        m0 = m if m0 is None else m0
        h = b["thickness"] - backing_thick
        burnt = m < 1e-5 * m0
        depth, T = b["depth"], b["T"]
        t.append(b["t"])
        mass.append(m)
        thick.append(h)
        Ttop.append(math.nan if burnt else T[0])
        Tback.append(math.nan if burnt else C.value_at_depth(depth, T, h))
        T1.append(math.nan if burnt else C.value_at_depth(depth, T, min(1e-3, h)))
        # UMD's "Top Temperature (1mm)" candidates: mean over the top 1 mm of
        # sample (element-width weighted), T at 0.5 mm, and a probe at a fixed
        # height L - 1 mm above the back face (the top element once passed).
        top = min(1e-3, h)
        acc = wsum = 0.0
        S = 0.0
        for j in range(n):
            a, z = S, S + widths[j]
            S = z
            ov = max(0.0, min(z, top) - a)
            acc += ov * T[j]
            wsum += ov
            if z >= top:
                break
        mean1 = acc / wsum if wsum > 0 else math.nan
        target = h - (L - 1e-3)
        probe = T[0] if target < 0 else C.value_at_depth(depth, T, target)
        top_flow = b["flow"].get("TOP", {}).get("gas", math.nan)
        aux.append((b["t"], b["mass"], top_flow, gas_mass, mean1,
                    C.value_at_depth(depth, T, min(0.5e-3, h)), probe))
    rows = C.series_from(t, mass, Ttop, Tback, T1, thick)
    return rows, aux


AUX_HEADER = ("time_s,total_mass_kg_m2,top_gas_outflow_kg_m2s,gas_in_sample_kg_m2,"
              "T_mean_top_1mm_K,T_0p5mm_K,T_fixed_height_L_minus_1mm_K")


def run_one(job):
    deck_name, variant, suffix, q, L_mm, t_end, dx, dt, seed = job
    deck = C.load_deck(deck_name)
    if variant == "rhoT":
        deck["rho_law"] = C.RHO_T_LAWS[deck_name]
    mode = "RAND" if variant == "rand" else "MAX"
    vdir = OUT / deck_name / (variant + suffix)
    case = f"run_{q}kW_{L_mm}mm"
    rdir = vdir / case
    rdir.mkdir(parents=True, exist_ok=True)
    (rdir / "case.cmp").write_text(cmp_text(deck, dx))
    (rdir / "case.cnd").write_text(cnd_text(deck, q, L_mm, t_end, dx, dt, mode))
    for f in ("case.out", "case.out.gz"):
        (rdir / f).unlink(missing_ok=True)  # ThermaKin refuses to overwrite its output
    env = dict(os.environ)
    exe = TK
    if mode == "RAND":
        exe = TK_SEEDED
        env["THERMAKIN_SEED"] = str(seed)
    rc, wall = C.timed_run([str(exe)], rdir, rdir / "run.log", "case.cmp\ncase.cnd\ncase.out\n", env=env)
    log = (rdir / "run.log").read_text(errors="replace")
    meta = dict(set=deck_name, variant=variant + suffix, q=q, L_mm=L_mm, t_end=t_end, dx=dx, dt=dt,
                mode=mode, seed=seed if mode == "RAND" else None, returncode=rc, wall_s=round(wall, 2),
                executable=str(exe.relative_to(C.REPO)))
    diverged = "integration diverged" in log
    if rc != 0 or ("Error" in log and not diverged) or not (rdir / "case.out").exists():
        meta["error"] = log[-400:]
        (rdir / "meta.json").write_text(json.dumps(meta, indent=1))
        raise RuntimeError(f"{deck_name} {variant}{suffix} {case}: rc={rc} {log[-300:]}")
    # A run that diverges is kept up to its last complete output (recorded).
    rows, aux = extract(rdir, "case", L_mm, allow_incomplete=diverged)
    meta["diverged"] = diverged
    with open(rdir / "case.out", "rb") as fi, gzip.open(rdir / "case.out.gz", "wb", compresslevel=6) as fo:
        shutil.copyfileobj(fi, fo)
    (rdir / "case.out").unlink()
    C.write_series(vdir / f"gasification_{q}kW_{L_mm}mm.csv", rows)
    with open(rdir / "aux.csv", "w") as f:
        f.write(AUX_HEADER + "\n")
        for r in aux:
            f.write(",".join(f"{x:.9g}" for x in r) + "\n")
    meta["t_last"] = rows[-1][0]
    meta["final_mass_g"] = rows[-1][1]
    (rdir / "meta.json").write_text(json.dumps(meta, indent=1))
    print(f"done {deck_name:10s} {variant+suffix:12s} {case:16s} wall {wall:8.1f} s", flush=True)
    return meta


def jobs_from(args):
    sets = args.sets.split(",") if args.sets else SETS
    variants = args.variants.split(",")
    wanted = set(args.cases.split(",")) if args.cases else None
    jobs = []
    for s in sets:
        for v in variants:
            if v == "rand" and s not in IN_DEPTH:
                continue  # one absorbing element either way for surface decks
            if v == "rhoT" and s not in C.RHO_T_LAWS:
                continue
            for q, L, t_end in C.gasification_cases():
                if wanted and f"{q}_{L}" not in wanted:
                    continue
                jobs.append((s, v, args.suffix, q, L, t_end, args.dx, args.dt, args.seed))
    # longest first
    jobs.sort(key=lambda j: -(j[5] * j[4] / j[6]))
    return jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["run", "extract", "decks"])
    ap.add_argument("--sets", default="")
    ap.add_argument("--variants", default="base")
    ap.add_argument("--cases", default="")
    ap.add_argument("--dx", type=float, default=5e-5)
    ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--suffix", default="")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()
    jobs = jobs_from(args)
    if args.action == "decks":
        for s, v, suffix, q, L, t_end, dx, dt, seed in jobs:
            d = C.load_deck(s)
            if v == "rhoT":
                d["rho_law"] = C.RHO_T_LAWS[s]
            print(cmp_text(d, dx))
            print(cnd_text(d, q, L, t_end, dx, dt, "RAND" if v == "rand" else "MAX"))
            break
        return
    if args.action == "extract":
        for s, v, suffix, q, L, t_end, dx, dt, seed in jobs:
            vdir = OUT / s / (v + suffix)
            rdir = vdir / f"run_{q}kW_{L}mm"
            rows, aux = extract(rdir, "case", L, allow_incomplete=True)
            C.write_series(vdir / f"gasification_{q}kW_{L}mm.csv", rows)
            with open(rdir / "aux.csv", "w") as f:
                f.write(AUX_HEADER + "\n")
                for r in aux:
                    f.write(",".join(f"{x:.9g}" for x in r) + "\n")
        return
    results = C.run_parallel(jobs, run_one, args.jobs)
    log = OUT / "runs.jsonl"
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "a") as f:
        for job, meta, err in results:
            f.write(json.dumps(dict(job=list(job), meta=meta, error=err)) + "\n")


if __name__ == "__main__":
    main()
