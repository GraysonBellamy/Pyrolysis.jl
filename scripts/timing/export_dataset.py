"""Write a campaign's matched-case times as the committed cost dataset read by
report mode (test/verification/reference/matched/timing/).

    python3 scripts/timing/export_dataset.py CAMPAIGN_DIR

external.csv: per external configuration, the time compared (FDS's internal
elapsed time, which excludes its MPI start-up; otherwise the process time) and
the process-time median and range. pyrolysis.csv: per case, rung and scored
curve/reference, this code's median warm construction-and-solve time and its
analytical error, at the report tolerances (tolerance 1) and at looser ones
(their multiplier). provenance.json: machine, commits, executables, protocol and
disturbance record. Python 3.11+, standard library only.
"""

from pathlib import Path
import csv
import json
import statistics
import sys
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / "test" / "verification" / "reference" / "matched" / "timing"


def jsonl(path):
    return [json.loads(x) for x in Path(path).read_text().splitlines() if x.strip()]


def main():
    camp = Path(sys.argv[1]).resolve()
    DEST.mkdir(parents=True, exist_ok=True)

    runs = [r for r in jsonl(camp / "external" / "external_runs.jsonl")]
    gp = camp / "external_gpyro" / "external_runs.jsonl"
    if gp.exists():
        runs += [r for r in jsonl(gp) if r["level"] != "committed"]
    groups = defaultdict(list)
    for r in runs:
        if r["repeat"] > 0:
            groups[(r["solver"], r["case"], r["level"])].append(r)
    with (DEST / "external.csv").open("w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["solver", "case", "level", "repeats", "compared_s", "process_median_s",
                    "process_min_s", "process_max_s", "internal_median_s"])
        for (s, c, lv), rs in sorted(groups.items()):
            assert all(r["success"] and r["outputs_identical"] for r in rs), (s, c, lv)
            proc = [r["wall_s"] for r in rs]
            internal = [r["internal_s"] for r in rs if r["internal_s"] is not None]
            internal_med = statistics.median(internal) if s == "fds" and internal else None
            compared = internal_med if internal_med is not None else statistics.median(proc)
            w.writerow([s, c, lv, len(rs), f"{compared:.6g}", f"{statistics.median(proc):.6g}",
                        f"{min(proc):.6g}", f"{max(proc):.6g}",
                        "" if internal_med is None else f"{internal_med:.6g}"])

    times = defaultdict(list)
    for r in csv.DictReader((camp / "matched" / "matched_times.csv").open()):
        times[(r["case"], r["variant"], int(r["n_cells"]), int(r["tolerance"]))].append(
            float(r["seconds"]))
    with (DEST / "pyrolysis.csv").open("w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["case", "variant", "n_cells", "tolerance", "repeats", "median_s", "min_s",
                    "max_s", "curve", "reference", "samples", "Linf", "RMS"])
        for r in csv.DictReader((camp / "matched" / "matched_errors.csv").open()):
            t = times[(r["case"], r["variant"], int(r["n_cells"]), int(r["tolerance"]))]
            w.writerow([r["case"], r["variant"], r["n_cells"], r["tolerance"], len(t),
                        f"{statistics.median(t):.6g}",
                        f"{min(t):.6g}", f"{max(t):.6g}", r["curve"], r["reference"], r["samples"],
                        r["Linf"], r["RMS"]])

    sessions = jsonl(camp / "external" / "external_sessions.jsonl")
    if gp.exists():
        sessions += jsonl(camp / "external_gpyro" / "external_sessions.jsonl")
    matched = [r for r in jsonl(camp / "paired_runs.jsonl") if r["name"] == "matched"][-1]
    timed = [r for g in groups.values() for r in g]
    env = sessions[0]["environment"]
    executables = {}
    for s in sessions:
        executables.update(s["executables"])
    for e in executables.values():
        e.pop("path", None)
    prov = dict(
        protocol="scripts/timing/README.md",
        measured=env["time"][:10],
        machine={k: env[k] for k in ("cpu_model", "timed_cpu", "timed_cpu_sibling", "governor",
                                     "energy_performance_preference", "platform_profile", "kernel")},
        pyrolysis=dict(commit=matched["commit"], tree=matched["tree"],
                       uncommitted_changes=matched["dirty"], julia="1.13.0", threads=1),
        executables=executables,
        executable_notes=dict(
            gpyro="rebuilt from the pinned source with the documented recipe; its hash differs from "
                  "the accuracy provenance only by the embedded build date, and every timed run "
                  "reproduced the committed native outputs byte for byte",
        ),
        timers=dict(
            fds="internal Total Elapsed Wall Clock Time (excludes MPI start-up, about 0.26 s here)",
            gpyro="process time",
            thermakin="process time (its own timer ticks in whole seconds)",
            pyrolysis="warm problem construction and solve, after a discarded compiling solve",
        ),
        repeats="5 when a run takes under 1 s (30 s for external codes), 3 under 10 s (300 s), else 1; median",
        disturbance=dict(
            external_max_sibling_busy=max(r["sibling_busy"] for r in timed),
            external_max_other_cores_busy=max(r["others_busy"] for r in timed),
            pyrolysis_sibling_busy=matched["sibling_busy"],
            pyrolysis_other_cores_busy=matched["others_busy"],
        ),
        outputs_identical_to_committed=True,
    )
    (DEST / "provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
    print(f"wrote {DEST}")


if __name__ == "__main__":
    main()
