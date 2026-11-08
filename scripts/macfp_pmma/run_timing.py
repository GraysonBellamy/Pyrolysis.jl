#!/usr/bin/env python3
"""Quiet-machine timing of the full MaCFP PMMA cases.

  python3 scripts/macfp_pmma/run_timing.py OUTDIR [--steps pyrolysis,external] [--wait 14400]

Follows scripts/timing/README.md: one timed process at a time, pinned to CPU 1
with its sibling idle, a quiet gate before every process, the disturbance of
every run recorded, warm-up then 5 / 3 / 1 repeats, the median reported.

- pyrolysis: time_pyrolysis.jl in one pinned process (it discards its own
  warm-up solve and repeats each rung of its mesh × tolerance ladder).
- external: FDS, Gpyro and ThermaKin2Ds on a representative subset (four sets ×
  three cases) at every level the external runners produced (base and refined
  cells or steps). Each repeat is a fresh pinned process in a fresh directory;
  its output must equal the original run's byte for byte.
Python 3.12+, standard library only. Run from the repository root.
"""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import shutil
import statistics
import sys

REPO = Path(__file__).resolve().parents[2]
HERE = REPO / "scripts" / "macfp_pmma"
sys.path.insert(0, str(REPO / "scripts" / "timing"))
sys.path.insert(0, str(HERE / "external"))
import quiet  # noqa: E402

SETS = ["NIST", "Aalto_I", "UMD", "GIDAZE+"]
CASES = ["10kW_6mm", "25kW_6mm", "65kW_12mm"]
EXE = {
    "fds": Path.home() / "FDS" / "FDS6" / "bin" / "fds",
    "gpyro": REPO / "docs" / "refs" / "gpyro_build" / "gpyro",
    "thermakin": REPO / "docs" / "refs" / "thermakin_build" / "thermakin",
}
# Inputs to copy, the command, stdin, and the output compared byte for byte.
SPEC = {
    "fds": dict(inputs=["gasification.fds"], args=["gasification.fds"], stdin=None,
                output="gasification_devc.csv"),
    "gpyro": dict(inputs=["case.data"], args=["case.data"], stdin=None, output="case_summary_01_0001.csv"),
    "thermakin": dict(inputs=["case.cmp", "case.cnd"], args=[], stdin="case.cmp\ncase.cnd\ncase.out\n",
                      output="case.out"),
}
LEVELS = {  # variant directories the runners wrote: base, and refined cell or step
    "fds": ["base", "base_c2.5e-5", "base_dt0.005"],
    "gpyro": ["base", "base_half"],
    "thermakin": ["base", "base_half"],
}


def repeats_for(seconds):
    return 5 if seconds < 1 else 3 if seconds < 10 else 1


def digest(path):
    p = Path(path)
    if not p.exists() and Path(str(p) + ".gz").exists():
        return hashlib.sha256(gzip.open(str(p) + ".gz", "rb").read()).hexdigest()
    return hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else None


def env_for(code):
    if code == "fds":
        import fds  # the runner's environment (Intel MPI paths)
        return fds.fds_env()
    return None


def time_external(out, wait):
    rows = []
    jl = out / "external_runs.jsonl"
    for code in ("fds", "gpyro", "thermakin"):
        spec = SPEC[code]
        for s in SETS:
            for level in LEVELS[code]:
                for case in CASES:
                    src = HERE / "results" / code / s / level / f"run_{case}"
                    if not all((src / f).exists() for f in spec["inputs"]):
                        continue
                    reference = digest(src / spec["output"])
                    times = []
                    k = 0
                    while True:
                        work = out / "work" / code / s / level / case / f"rep{k}"
                        shutil.rmtree(work, ignore_errors=True)
                        work.mkdir(parents=True)
                        for f in spec["inputs"]:
                            shutil.copy2(src / f, work / f)
                        quiet.wait_quiet(timeout=wait)
                        r = quiet.run_pinned([str(EXE[code]), *spec["args"]], cwd=work, log_path=work / "run.log",
                                             env=env_for(code), stdin_text=spec["stdin"], timeout=7200)
                        same = digest(work / spec["output"]) == reference
                        r.update(code=code, set=s, level=level, case=case, repeat=k, warmup=(k == 0),
                                 identical_output=same)
                        quiet.append_jsonl(jl, r)
                        print(f"{code:9s} {s:8s} {level:13s} {case:10s} rep{k} {r['wall_s']:8.2f} s "
                              f"identical={same} rc={r['returncode']}", flush=True)
                        if k > 0:
                            times.append(r["wall_s"])
                        if k >= repeats_for(r["wall_s"]):
                            break
                        k += 1
                        if k > 1:
                            shutil.rmtree(out / "work" / code / s / level / case / f"rep{k - 2}", ignore_errors=True)
                    rows.append(dict(code=code, set=s, level=level, case=case, median_s=statistics.median(times),
                                     min_s=min(times), max_s=max(times), repeats=len(times)))
    with open(out / "external_times.csv", "w") as f:
        f.write("code,set,level,case,median_s,min_s,max_s,repeats\n")
        for r in rows:
            f.write("{code},{set},{level},{case},{median_s:.4f},{min_s:.4f},{max_s:.4f},{repeats}\n".format(**r))


def time_pyrolysis(out, wait, sets=(), name="pyrolysis_timing.csv"):
    quiet.wait_quiet(timeout=wait)
    cmd = ["julia", "--project=scripts/macfp_pmma", "--threads=1", "--startup-file=no",
           "scripts/macfp_pmma/time_pyrolysis.jl", str(out / name), *sets]
    r = quiet.run_pinned(cmd, cwd=REPO, log_path=out / "pyrolysis_timing.log", timeout=6 * 3600)
    r.update(step="pyrolysis")
    quiet.append_jsonl(out / "pyrolysis_run.jsonl", r)
    print(f"pyrolysis ladder: rc={r['returncode']} wall {r['wall_s']:.0f} s, "
          f"sibling busy {r['sibling_busy']:.3f}, others {r['others_busy']:.3f}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--steps", default="pyrolysis,external")
    ap.add_argument("--wait", type=float, default=4 * 3600, help="longest wait for a quiet machine [s]")
    ap.add_argument("--sets", default="", help="comma-separated deck ids for the Pyrolysis.jl ladder (default: all)")
    ap.add_argument("--name", default="pyrolysis_timing.csv", help="CSV name for the Pyrolysis.jl ladder")
    a = ap.parse_args()
    out = Path(a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    quiet.keep_driver_off()
    env = quiet.environment()
    env["git"] = quiet.git_stamp(REPO)
    (out / "environment.json").write_text(json.dumps(env, indent=2))
    try:
        quiet.check_environment(env)
    except RuntimeError as e:      # record and go on: the owner judges the record
        print(f"WARNING {e}", flush=True)
        (out / "environment_warning.txt").write_text(str(e) + "\n")
    steps = a.steps.split(",")
    if "pyrolysis" in steps:
        time_pyrolysis(out, a.wait, [x for x in a.sets.split(",") if x], a.name)
    if "external" in steps:
        time_external(out, a.wait)


if __name__ == "__main__":
    main()
