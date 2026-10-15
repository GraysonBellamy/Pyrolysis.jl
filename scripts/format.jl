#!/usr/bin/env julia
# Format (or check the formatting of) the Pyrolysis.jl codebase with JuliaFormatter.
#
# Usage:
#   julia --project=scripts scripts/format.jl            # format files in place
#   julia --project=scripts scripts/format.jl --check    # verify only; exit 1 if any file needs formatting
#
# Style is read from .JuliaFormatter.toml at the repo root.
using JuliaFormatter

const ROOT = dirname(@__DIR__)
const TARGETS = ("src", "test", "ext", "examples", "benchmark", "scripts")

check = "--check" in ARGS
paths = String[joinpath(ROOT, t) for t in TARGETS if isdir(joinpath(ROOT, t))]

# `format` returns true when every file was already correctly formatted.
# With overwrite=false (check mode) it reports without modifying files.
# `throw_on_error` makes a file the formatter cannot process fail the run;
# otherwise it is logged and counted as already formatted.
already_formatted = format(paths; overwrite = !check, verbose = true, throw_on_error = true)

if check && !already_formatted
    @error "Formatting check failed — run `julia --project=scripts scripts/format.jl` to fix."
    exit(1)
end

println(
    check ? "Format check passed: all files already formatted." : "Formatting complete.",
)
