# Timing protocol

Wall times on this repository's development machine (Intel Core Ultra 7 155H:
performance cores 0–11 as hyperthread pairs, efficiency cores 12–19,
low-power cores 20–21) move by tens of percent with load, core type and
thermal state. Every timing in the guides, the gate records and the
cross-code comparisons follows this protocol so that it can be repeated and
judged afterwards.

## Rules

1. **One timed process at a time, pinned to CPU 1, with its sibling CPU 2
   idle.** CPUs 1–2 are a 4.8 GHz performance core. The drivers move
   themselves onto the efficiency cores (`quiet.keep_driver_off`) and block
   in `wait4` while a timed process runs.
2. **A quiet gate before every timed process:** one-minute load below 2, and,
   over 0.3 s, the sibling at most 5 % and the other cores on average at most
   10 % busy (`quiet.wait_quiet`). Nothing is compiled or built while a
   campaign runs; every environment is instantiated and precompiled first.
3. **Every run records its disturbance:** busy fractions of the timed core,
   its sibling and the rest over the run (from `/proc/stat`), involuntary
   context switches, the load and the package temperature before it. Runs are
   kept and judged on that record, never dropped silently.
4. **The machine configuration is recorded and checked:** AC power, turbo on,
   platform profile `performance`, the timed CPU a performance core
   (`quiet.environment`, `quiet.check_environment`).
5. **Warm-up, then repeats, then the median.** Julia: one discarded solve per
   problem type (compilation), then 5 solves when one takes under 1 s, 3 under
   10 s, otherwise 1, each after `GC.gc()`. External codes: one first pass per
   configuration (binary and page cache warm-up, and the duration that sizes
   the repeats), then fresh processes with the same 5/3/1 rule. Report the
   median with the range.
6. **Comparisons between commits are paired and interleaved** (`paired.py`):
   rounds A B, B A, A B, ..., so drift falls on both sides. Retired
   instructions (`benchmark/instructions`) are the load-independent second
   reading; steps, residuals and allocations must be reported beside any time.
7. **Times between codes are compared only at matched accuracy.** Every timed
   external run must reproduce the committed outputs byte for byte, so its time
   belongs to the computation whose error the comparison pages report.
   `pyrolysis_matched.jl` scores this code on the same curves and analytical
   references over its resolution ladder and looser solver tolerances, the
   counterpart of the external codes' coarser steps. A cost comparison reads
   both at the same error, never at the same mesh or tolerance.
8. **Measure commits, not working trees.** Run from clean detached worktrees
   (`git worktree add --detach`), and record the commit and tree.

## Running a campaign

```sh
# Clean checkouts of the commits to compare, environments built first:
git worktree add --detach ../wt_new <commit>
git worktree add --detach ../wt_base <commit>
for wt in ../wt_new ../wt_base; do
  (cd $wt && julia --project=test -e 'using Pkg; Pkg.instantiate(); using Pyrolysis' &&
   julia --project=benchmark/instructions -e 'using Pkg; Pkg.instantiate()')
done
(cd ../wt_new && julia --project=scripts/verification -e 'using Pkg; Pkg.instantiate(); using Pyrolysis')

# Then, on a quiet machine on AC power (about 2–3 hours for everything):
FDS=~/FDS/FDS6/bin/fds GPYRO=/path/to/gpyro THERMAKIN=docs/refs/thermakin_build/thermakin \
  scripts/timing/campaign.sh OUTDIR ../wt_new ../wt_base
```

The steps, each serial and pinned:

| step | what | output |
|---|---|---|
| `recheck` | warm solves of `test/gates/performance_recheck.jl`, 3 paired rounds | `recheck/`, `paired_runs.jsonl` |
| `instructions` | retired instructions of `benchmark/instructions`, 1 paired round | `instructions/` |
| `ttfx` | first-solve latency of `test/gates/ttfx.jl`, 2 paired rounds | `ttfx/` |
| `external` | FDS, Gpyro V2 and ThermaKin2Ds on every matched deck at all seven levels | `external/external_summary.csv` |
| `matched` | this code on the matched cases' curves over its ladder, at the report tolerances and 10×, 100×, 1000× looser | `matched/matched_{times,errors}.csv` |

The other committed reference decks are timed by `external.py --deck-runs`,
for example `fds:test/verification/reference/fds/custom/*.fds`. The verification pages'
`wall (s)` columns come from
`scripts/verification/run_reports.jl --timings=true`, which applies the same
repeat rule. Run it pinned on the same quiet machine with `taskset -c 1`.

The executables: FDS is the official 6.11.1 Linux bundle. Gpyro V2 and
ThermaKin2Ds are built as described in
`test/verification/reference/matched/README.md` and
`test/verification/reference/matched/build_thermakin.sh`. Their hashes must
equal the provenance files beside the matched data. `external.py` records them.
