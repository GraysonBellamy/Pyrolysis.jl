"""Tables from a timing campaign's raw records (campaign.sh OUTDIR).

    python3 scripts/timing/analyze.py OUTDIR [--accuracy-dir docs/guide/verification/assets]

Writes OUTDIR/report.md: the paired comparisons (warm solves, retired
instructions, first-solve latency), the external codes' times with their
identity and disturbance checks, and the cross-code cost at matched accuracy.

Cost at matched accuracy: for an external run at one resolution, this code's
cost is that of its cheapest ladder rung that is at least as accurate (Linf
against the same analytical reference) on every curve and reference the
external run is scored on. "≤" marks a rung that is the coarsest timed (the
cost to match could be lower); "not reached" means no rung was as accurate.
Python 3.11+, standard library only.
"""

from pathlib import Path
import argparse
import csv
import json
import re
import statistics
import tomllib
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[2]
MATCHED = ROOT / "test" / "verification" / "reference" / "matched"
LEVELS = ["base", "time_half", "time_quarter", "space_half", "refined"]


def jsonl(path):
    return [json.loads(x) for x in Path(path).read_text().splitlines() if x.strip()] if Path(path).exists() else []


def fmt(x, digits=4):
    if x is None or x == "":
        return "—"
    if isinstance(x, str):
        return x
    return f"{x:.{digits}g}"


def disturbance_note(rows):
    sib = max(r["sibling_busy"] for r in rows)
    oth = max(r["others_busy"] for r in rows)
    load = max(r["load_before"] for r in rows)
    return f"max sibling busy {sib:.3f}, max other cores busy {oth:.3f}, max load before {load:.2f}"


def paired_recheck(out, lines):
    d = out / "recheck"
    if not d.exists():
        return
    runs = [r for r in jsonl(out / "paired_runs.jsonl") if r["name"] == "recheck"]
    names = {}
    samples = defaultdict(list)  # (label, config) -> [(round, seconds, bytes, steps)]
    for f in sorted(d.glob("*_performance_recheck.csv")):
        label, rnd = re.match(r"(.+)_r(\d+)_performance_recheck\.csv", f.name).groups()
        for row in csv.DictReader(f.open()):
            samples[(label, int(row["configuration"]))].append(
                (int(rnd), float(row["seconds"]), int(row["bytes"]), int(row["accepted_steps"])))
    first_md = next(iter(sorted(d.glob("*_performance_recheck.md"))), None)
    if first_md is not None:
        rows = [re.match(r"\| ([^|]+) \| [0-9.]+ \|", x) for x in first_md.read_text().splitlines()]
        names = {i + 1: m[1].strip() for i, m in enumerate(m for m in rows if m)}
    labels = sorted({k[0] for k in samples})
    stamps = {r["label"]: (r["commit"][:8], r["tree"][:8]) for r in runs}
    lines += ["## Warm solves, paired (test/gates/performance_recheck.jl)", ""]
    lines += [", ".join(f"{l}: commit {stamps[l][0]} (tree {stamps[l][1]})" for l in labels if l in stamps) + ".",
              f"{max(r['round'] for r in runs)} interleaved rounds per checkout, five solves each after one warm-up "
              f"(pooled below); {disturbance_note(runs)}.", ""]
    if len(labels) == 2 and "new" in labels and "base" in labels:
        lines += ["| configuration | base median [s] (range) | new median [s] (range) | new/base | per-round ratios | steps base, new | MiB base, new |",
                  "|---|---|---|---|---|---|---|"]
        for c in sorted({k[1] for k in samples}):
            b, n = samples[("base", c)], samples[("new", c)]
            bt, nt = [x[1] for x in b], [x[1] for x in n]
            rounds = sorted({x[0] for x in b} & {x[0] for x in n})
            per = [statistics.median(x[1] for x in n if x[0] == k) / statistics.median(x[1] for x in b if x[0] == k)
                   for k in rounds]
            lines.append(
                f"| {names.get(c, c)} | {fmt(statistics.median(bt))} ({fmt(min(bt))}–{fmt(max(bt))}) | "
                f"{fmt(statistics.median(nt))} ({fmt(min(nt))}–{fmt(max(nt))}) | "
                f"{statistics.median(nt) / statistics.median(bt):.4f} | {', '.join(f'{p:.3f}' for p in per)} | "
                f"{sorted({x[3] for x in b})}, {sorted({x[3] for x in n})} | "
                f"{statistics.median(x[2] for x in b) / 2**20:.2f}, {statistics.median(x[2] for x in n) / 2**20:.2f} |")
    lines.append("")


def paired_instructions(out, lines):
    d = out / "instructions"
    files = {re.match(r"(.+)_r\d+_instructions\.csv", f.name)[1]: f for f in d.glob("*_instructions.csv")} if d.exists() else {}
    if not {"new", "base"} <= set(files):
        return
    read = lambda f: {r["operation"]: r for r in csv.DictReader(f.open())}
    b, n = read(files["base"]), read(files["new"])
    lines += ["## Retired instructions, paired (benchmark/instructions)", "",
              "| operation | base instructions | new/base instructions | new/base cycles | bytes base, new |",
              "|---|---|---|---|---|"]
    for op in b:
        if op in n:
            lines.append(f"| {op} | {float(b[op]['instructions']):.4g} | "
                         f"{float(n[op]['instructions']) / float(b[op]['instructions']):.4f} | "
                         f"{float(n[op]['cycles']) / float(b[op]['cycles']):.3f} | {b[op]['bytes']}, {n[op]['bytes']} |")
    lines.append("")


def paired_ttfx(out, lines):
    d = out / "ttfx"
    if not d.exists():
        return
    vals = defaultdict(list)
    for f in d.glob("*_ttfx.md"):
        label = re.match(r"(.+)_r\d+_ttfx\.md", f.name)[1]
        for line in f.read_text().splitlines():
            m = re.match(r"\| ([^|]+) \| ([0-9.]+) \| ([0-9.]+) \|", line)
            if m:
                vals[(label, m[1].strip())].append((float(m[2]), float(m[3])))
    configs = sorted({k[1] for k in vals}, key=lambda x: (x.startswith("Douglas"), x))
    lines += ["## First-solve latency, paired (test/gates/ttfx.jl)", "",
              "Fresh process per row; medians over rounds; the harness rounds to 0.01 s.", "",
              "| configuration | load base, new [s] | first solve base, new [s] | new/base |", "|---|---|---|---|"]
    for c in configs:
        b, n = vals.get(("base", c), []), vals.get(("new", c), [])
        if not b or not n:
            continue
        bs, ns = statistics.median(x[1] for x in b), statistics.median(x[1] for x in n)
        lines.append(f"| {c} | {statistics.median(x[0] for x in b):.2f}, {statistics.median(x[0] for x in n):.2f} | "
                     f"{bs:.2f}, {ns:.2f} | {ns / bs:.3f} |" if bs > 0 else f"| {c} | | {bs}, {ns} | |")
    lines.append("")


def external(out, lines):
    dirs = sorted(d for d in out.glob("external*") if (d / "external_runs.jsonl").exists())
    runs = [r for d in dirs for r in jsonl(d / "external_runs.jsonl")]
    if not runs:
        return {}
    timed = [r for r in runs if r["repeat"] > 0]
    groups = defaultdict(list)
    for r in timed:
        groups[(r["solver"], r["case"], r["level"])].append(r)
    bad = [k for k, rs in groups.items() if not all(r["success"] and r["outputs_identical"] for r in rs)]
    lines += ["## External codes", "",
              f"{len(timed)} timed processes over {len(groups)} configurations after one first pass each; "
              f"{disturbance_note(timed)}. Outputs differing from the committed ones: "
              f"{', '.join(map(str, bad)) if bad else 'none'}.", ""]
    seen = set()
    for d in dirs:
        for s in jsonl(d / "external_sessions.jsonl"):
            for solver, e in s["executables"].items():
                if (solver, e["sha256"]) not in seen:
                    seen.add((solver, e["sha256"]))
                    lines.append(f"- {solver}: `{e['sha256'][:16]}…` ({e['path']})")
    lines.append("")
    times = {}
    for k, rs in groups.items():
        times[k] = dict(
            process=statistics.median(r["wall_s"] for r in rs),
            internal=statistics.median(r["internal_s"] for r in rs) if all(r["internal_s"] for r in rs) else None,
            user=statistics.median(r["user_s"] for r in rs),
            n=len(rs), spread=(max(r["wall_s"] for r in rs) - min(r["wall_s"] for r in rs)) / statistics.median(r["wall_s"] for r in rs))
    lines += ["| solver | decks | configurations | median process [s] range | largest spread / median |",
              "|---|---|---|---|---|"]
    for solver, kind in sorted({(k[0], k[2] == "committed") for k in times}):
        ks = [k for k in times if k[0] == solver and (k[2] == "committed") == kind]
        lines.append(f"| {solver} | {'committed' if kind else 'matched'} | {len(ks)} | {fmt(min(times[k]['process'] for k in ks))}–{fmt(max(times[k]['process'] for k in ks))} | "
                     f"{max(times[k]['spread'] for k in ks):.3f} |")
    lines.append("")
    return times


def crosscode(out, lines, times, accuracy_dir):
    ours_path = out / "matched" / "matched_errors.csv"
    if not ours_path.exists() or not times:
        return
    ours = defaultdict(dict)  # case -> n -> {"t": s, (curve, ref): Linf}
    for r in csv.DictReader(ours_path.open()):
        key = r["case"] if not r["variant"] else f"{r['case']}:{r['variant']}"
        e = ours[key].setdefault(int(r["n_cells"]), {"t": float(r["median_seconds"])})
        e[(r["curve"], r["reference"])] = float(r["Linf"])
    entries = []
    for solver in ("fds", "gpyro", "thermakin"):
        entries += [dict(e, solver=solver) for e in tomllib.load((MATCHED / f"{solver}.toml").open("rb"))["series"]]
    acc = {}
    for case in {e["case"] for e in entries}:
        f = accuracy_dir / case / f"{case}_comparison_accuracy.csv"
        for r in csv.DictReader(f.open()):
            acc[(case, r["code"], r["curve"], r["reference"], r["resolution"])] = float(r["Linf"])
    lines += ["## Cost at matched accuracy", "",
              "External time: FDS's internal elapsed time (it excludes MPI start-up), else the process time; "
              "process times beside. This code: median warm construction-and-solve time of its "
              "cheapest ladder rung at least as accurate on every scored curve and reference. "
              "Ratio > 1: this code is cheaper.", "",
              "| case | solver | level | external Linf (worst ref) | external [s] (process) | our rung | our [s] | ratio |",
              "|---|---|---|---|---|---|---|---|"]
    by_run = defaultdict(list)
    for e in entries:
        by_run[(e["case"], e["solver"], Path(e["path"]).stem)].append(e)
    for (case, solver, mcase), es in sorted(by_run.items()):
        for level in LEVELS:
            t = times.get((solver, mcase, level))
            if t is None:
                continue
            need = {}
            for e in es:
                v = acc.get((case, e["label"], e["curve"], e["reference"], level))
                if v is not None:
                    need[(e["curve"], e["reference"])] = v
            if not need:
                continue
            # V-CPL-03's paths are separate problems on our side: matching all
            # of them costs one solve of each.
            if case == "V-CPL-03":
                ladders = [ours.get(f"{case}:{ref}", {}) for (_, ref) in need]
                rung_sets = []
                for (curve, ref), lad in zip(need, ladders):
                    rung_sets.append({n: v["t"] for n, v in lad.items() if v.get((curve, ref), float("inf")) <= need[(curve, ref)]})
                common = set.intersection(*(set(s) for s in rung_sets)) if rung_sets else set()
                cand = [(sum(s[n] for s in rung_sets), n) for n in common]
                all_ns = set.intersection(*(set(l) for l in ladders)) if ladders else set()
            else:
                lad = ours.get(case, {})
                cand = [(v["t"], n) for n, v in lad.items()
                        if all(v.get(k, float("inf")) <= linf for k, linf in need.items())]
                all_ns = set(lad)
            # FDS's internal timer excludes its MPI start-up (about 0.26 s here),
            # as this code's warm time excludes Julia's start-up; ThermaKin's
            # timer ticks in whole seconds, and its launch costs milliseconds.
            ext = t["internal"] if solver == "fds" and t["internal"] else t["process"]
            worst = max(need.values())
            if cand:
                cost, n = min(cand)
                bound = "≤ " if n == min(all_ns) else ""
                ratio = f"{'≥ ' if bound else ''}{ext / cost:.3g}"
                lines.append(f"| {case} | {solver} | {level} | {worst:.3g} | {fmt(ext)} ({fmt(t['process'])}) | "
                             f"{n} | {bound}{fmt(cost)} | {ratio} |")
            else:
                lines.append(f"| {case} | {solver} | {level} | {worst:.3g} | {fmt(ext)} ({fmt(t['process'])}) | "
                             f"not reached | — | — |")
    lines.append("")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir", type=Path)
    ap.add_argument("--accuracy-dir", type=Path, default=ROOT / "docs" / "guide" / "verification" / "assets")
    a = ap.parse_args()
    out = a.outdir.resolve()
    lines = ["# Timing campaign", ""]
    env = jsonl(out / "paired_sessions.jsonl")
    if env:
        e = env[0]["environment_before"]
        lines += [f"{e['cpu_model']}, CPU {e['timed_cpu']} (sibling {e['timed_cpu_sibling']}), governor "
                  f"{e['governor']}, EPP {e['energy_performance_preference']}, platform profile "
                  f"{e['platform_profile']}, AC {e['ac_online']}, kernel {e['kernel']}. Protocol: scripts/timing/README.md.", ""]
    paired_recheck(out, lines)
    paired_instructions(out, lines)
    paired_ttfx(out, lines)
    times = external(out, lines)
    crosscode(out, lines, times, a.accuracy_dir)
    (out / "report.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
