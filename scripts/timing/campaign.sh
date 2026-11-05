#!/bin/bash
# One quiet-machine timing campaign, every step serial and pinned (README.md).
#
#   scripts/timing/campaign.sh OUTDIR NEW_CHECKOUT BASE_CHECKOUT [steps...]
#
# NEW/BASE are clean checkouts (git worktree add --detach) whose test and
# benchmark/instructions environments are already instantiated and precompiled:
# compilation must not overlap a timed step. Steps (default: all, in order):
#   recheck       warm solves of test/gates/performance_recheck.jl, 3 paired rounds
#   instructions  retired instructions of benchmark/instructions, 1 paired round
#   ttfx          first-solve latency of test/gates/ttfx.jl, 2 paired rounds
#   external      FDS, Gpyro and ThermaKin on the matched decks (paths: FDS, GPYRO,
#                 THERMAKIN environment variables; a missing one is skipped)
#   matched       Pyrolysis.jl on the matched cases (scripts/timing/pyrolysis_matched.jl)
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
out=$(realpath -m "$1"); new=$(realpath "$2"); base=$(realpath "$3"); shift 3
steps=${*:-recheck instructions ttfx external matched}
mkdir -p "$out"
julia_run=(julia --project=test --startup-file=no --threads=1)

for step in $steps; do
  echo "=== $step $(date -Is)"
  case $step in
    recheck)
      python3 "$here/paired.py" "$out" recheck --rounds 3 --tree new="$new" --tree base="$base" \
        --artifact test/gates/measurements/performance_recheck.csv \
        --artifact test/gates/measurements/performance_recheck.md \
        -- "${julia_run[@]}" test/gates/performance_recheck.jl ;;
    instructions)
      python3 "$here/paired.py" "$out" instructions --rounds 1 --tree new="$new" --tree base="$base" \
        --artifact test/gates/measurements/instructions.csv \
        -- "${julia_run[@]}" benchmark/instructions/instructions.jl \
           --record test/gates/measurements/instructions.csv ;;
    ttfx)
      # Every row is a fresh child process; children inherit the pinning.
      python3 "$here/paired.py" "$out" ttfx --rounds 2 --tree new="$new" --tree base="$base" \
        --artifact test/gates/measurements/ttfx.md \
        -- julia --project=test --startup-file=no test/gates/ttfx.jl ;;
    external)
      solvers=()
      [[ -n "${FDS:-}" ]] && solvers+=(fds)
      [[ -n "${GPYRO:-}" ]] && solvers+=(gpyro)
      [[ -n "${THERMAKIN:-}" ]] && solvers+=(thermakin)
      if (( ${#solvers[@]} )); then
        args=()
        [[ -n "${FDS:-}" ]] && args+=(--fds "$FDS")
        [[ -n "${GPYRO:-}" ]] && args+=(--gpyro "$GPYRO")
        [[ -n "${THERMAKIN:-}" ]] && args+=(--thermakin "$THERMAKIN")
        python3 "$here/external.py" "$out/external" "${args[@]}" --solvers "$(IFS=,; echo "${solvers[*]}")"
      else
        echo "no external executables given; skipped"
      fi ;;
    matched)
      python3 "$here/paired.py" "$out" matched --rounds 1 --tree new="$new" \
        -- julia --project=scripts/verification --startup-file=no --threads=1 \
           "$here/pyrolysis_matched.jl" "$out/matched" ;;
    *) echo "unknown step $step"; exit 2 ;;
  esac
done
echo "=== done $(date -Is)"
