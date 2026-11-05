"""Run one measurement command in several checkouts, interleaved, for a paired
comparison: round k runs every checkout once, in alternating order (A B, B A,
A B, ...), so slow drift of the machine falls on both sides equally. Each
process is pinned and recorded by quiet.run_pinned; the command's own result
files (`--artifact`, relative to the checkout) are copied after each run.

    python3 scripts/timing/paired.py OUTDIR NAME --rounds 3 \
        --tree new=/path/to/checkout --tree base=/path/to/other \
        [--artifact test/gates/measurements/performance_recheck.csv] \
        -- julia --project=test --startup-file=no --threads=1 test/gates/performance_recheck.jl

Writes OUTDIR/NAME/<label>_r<k>.log and .<artifact name>, and appends one row
per process to OUTDIR/paired_runs.jsonl. Python 3.11+, standard library only.
"""

from pathlib import Path
import argparse
import os
import shutil
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
import quiet  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir", type=Path)
    ap.add_argument("name")
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--tree", action="append", required=True, help="label=checkout")
    ap.add_argument("--artifact", action="append", default=[])
    ap.add_argument("--allow-dirty", action="store_true")
    argv = sys.argv[1:]
    if "--" not in argv:
        raise SystemExit("give the command after --")
    cut = argv.index("--")
    a = ap.parse_args(argv[:cut])
    command = argv[cut + 1:]
    if not command:
        raise SystemExit("no command")
    trees = [(t.split("=", 1)[0], Path(t.split("=", 1)[1]).resolve()) for t in a.tree]
    out = (a.outdir / a.name).resolve()
    out.mkdir(parents=True, exist_ok=True)
    quiet.keep_driver_off()
    env0 = quiet.environment()
    quiet.check_environment(env0)
    stamps = {label: quiet.git_stamp(path) for label, path in trees}
    for label, s in stamps.items():
        # The gates write their tables under test/gates/measurements, which
        # their own stamp ignores; anything else must be clean.
        dirty = [x for x in os.popen(f"git -C {dict(trees)[label]} status --porcelain").read().splitlines()
                 if "test/gates/measurements/" not in x]
        if dirty and not a.allow_dirty:
            raise SystemExit(f"{label}: checkout has changes: {dirty[:5]}")
    env = {**os.environ, "JULIA_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"}
    for k in range(1, a.rounds + 1):
        order = trees if k % 2 else list(reversed(trees))
        for label, path in order:
            quiet.wait_quiet()
            log = out / f"{label}_r{k}.log"
            rec = quiet.run_pinned(command, path, log, env=env)
            for art in a.artifact:
                src = path / art
                if src.exists():
                    shutil.copyfile(src, out / f"{label}_r{k}_{Path(art).name}")
            row = dict(name=a.name, label=label, round=k, checkout=str(path), command=command,
                       time=time.strftime("%Y-%m-%dT%H:%M:%S"), **stamps[label], **rec)
            quiet.append_jsonl(a.outdir.resolve() / "paired_runs.jsonl", row)
            print(f"{a.name} {label} r{k}: exit {rec['returncode']}, {rec['wall_s']:.1f} s, "
                  f"sibling {rec['sibling_busy']:.3f}, others {rec['others_busy']:.3f}", flush=True)
            if rec["returncode"] != 0:
                raise SystemExit(f"{label} round {k} failed; see {log}")
    quiet.append_jsonl(a.outdir.resolve() / "paired_sessions.jsonl",
                       dict(name=a.name, environment_before=env0, environment_after=quiet.environment(),
                            trees={l: str(p) for l, p in trees}, stamps=stamps, argv=sys.argv))


if __name__ == "__main__":
    main()
