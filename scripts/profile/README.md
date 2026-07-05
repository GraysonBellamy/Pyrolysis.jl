# Deep profiling harness

Answers "*where* is the time/allocation going?" — the complement to
[`benchmark/`](../../benchmark/) (which answers "*did* it get slower?").

Runs the canonical configurations from
[`test/perf/perf_helpers.jl`](../../test/perf/perf_helpers.jl) under:

- **`summary`** — wall clock, bytes allocated, and GC share per residual
  evaluation, per Jacobian assembly, and per end-to-end solve (plus
  integrator stats: `nf`, `njacs`, `naccept`, `nreject`).
- **`cpu`** — `Profile.@profile` flamegraphs exported as pprof `.pb.gz`
  files, with a flat top-frames teaser printed to the terminal.
- **`allocs`** — `Profile.Allocs.@profile` allocation-site profiles,
  aggregated by Pyrolysis `file:line` and by type in the terminal, full
  detail in the pprof file.

## Running

From the package root (first run instantiates the isolated env):

```sh
julia --project=scripts/profile -e 'using Pkg; Pkg.instantiate()'
julia --project=scripts/profile scripts/profile/profile.jl
```

Defaults profile the heaviest configuration (`ale_beer_lambert_darcy`) in
all modes and targets. Narrow it down with:

```sh
# Allocation hunt in the full solve of the PMMA cone case, 20 s horizon
julia --project=scripts/profile scripts/profile/profile.jl \
  --config=eulerian_surface_absorption --target=solve --mode=allocs --t-final=20

# CPU flamegraph of just the structured Jacobian, every config
julia --project=scripts/profile scripts/profile/profile.jl \
  --config=all --target=jacobian --mode=cpu
```

Options: `--config`, `--mode`, `--target` (each accepts `all`),
`--t-final` (solve horizon, s), `--seconds` (wall-clock budget per hot-op
profile window), `--alloc-rate` (allocation sample rate for the solve
target; hot ops always sample at 1.0), `--strategy` (Structured-backend
strategy: `direct` (default, full assembly — the verification anchor),
`approx` (the `:ale_reduced` profile that drops the dense ALE couplings —
the production target for ALE), `structured`, `matrixfree`). Non-default
strategies are suffixed into the output filenames, e.g.
`ale_no_radiation__jacobian__approx__cpu.pb.gz`, so direct/approx profiles
can sit side by side.

```sh
# Quantify what ApproxSolve buys on the ALE Jacobian
julia --project=scripts/profile scripts/profile/profile.jl \
  --config=ale_no_radiation --mode=summary --strategy=approx
```

## Inspecting profiles

Outputs land in `scripts/profile/output/` (gitignored). Open interactively:

```sh
julia --project=scripts/profile -e \
  'using PProf; PProf.refresh(file="scripts/profile/output/<file>.pb.gz"); readline()'
```

then browse the flamegraph at the printed URL (Flame Graph / Top / Graph
views; for allocation profiles switch the sample type between
`size` and `count`).

## Memory safety

Allocation profiling stores a full backtrace per sampled allocation — an
early unbounded run was OOM-killed at 22 GiB RSS while converting millions
of records to pprof format. The harness now bounds this: hot-op windows
are sized from a measured records-per-eval rate (~200k records max), the
pprof conversion only sees the capped slice, the solve target samples at
`--alloc-rate=0.01` by default, and profiling is skipped outright when
free RAM is under ~4 GiB. If you raise `--alloc-rate`, consider wrapping
the run in a cgroup, e.g.:

```sh
systemd-run --user --scope -p MemoryMax=12G -p MemorySwapMax=2G \
  julia --project=scripts/profile scripts/profile/profile.jl ...
```

## Caveats

- The wood/Darcy `solve` cases use a `ρ_intrinsic = 1500` variant of the
  hot-op fixture: the perf-gate fixture has no skeletal density (φ₀ = 0),
  which `solve` rejects at validation for permeable materials. The
  residual/jacobian targets keep the exact gate fixture.

- Residual/Jacobian hot ops are evaluated at the **initial state** — the
  same convention as the zero-allocation gates in `test/perf/`. Branches
  that only activate mid-burn (e.g. deep char) show up in the `solve`
  target instead.
- The end-to-end solve includes a short warm-up solve so compilation is
  excluded from profiles; first invocation of the script still pays the
  usual package-load latency.
