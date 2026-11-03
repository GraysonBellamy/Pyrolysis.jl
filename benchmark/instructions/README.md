# Instruction-count benchmark

Time is unusable on a loaded machine, and retired user-space instructions are
not: they repeat to within a fraction of a percent whatever else is running.
`instructions.jl` counts them, with cycles and bytes beside, for one residual
pass, one Jacobian pass and a short solve of four decks in both energy forms.

```sh
# Once: this directory holds the counters' package alone; it stacks on the test
# environment, with which it does not resolve.
julia --project=benchmark/instructions -e 'using Pkg; Pkg.instantiate()'

# A record, then a comparison in a later tree (same Julia, same manifest):
taskset -c 0-11 julia --project=test --startup-file=no --threads=1 \
    benchmark/instructions/instructions.jl --record before.csv
taskset -c 0-11 julia --project=test --startup-file=no --threads=1 \
    benchmark/instructions/instructions.jl --compare before.csv
```

**Pin to the cores with counters.** On a hybrid CPU the event may read exactly
zero on efficiency cores; the script stops when it does. `perf_event_paranoid`
must be 2 or lower; no root and no `perf` binary are needed.

**Read the columns in order.** Instructions are the figure of merit. Cycles are a
second opinion for a change to data layout, or one that trades arithmetic for
memory: fewer instructions with more cycles is not a win, and on a short run
under load cycles scatter by tens of percent. Fewer steps or residuals is not a
gain unless the error held; this benchmark does not measure error.
