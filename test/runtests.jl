using Test
using Pyrolysis

# Each test file is evaluated in its own throwaway module, so files that define
# top-level structs and constants cannot clash. The module name is derived from
# the path so new files are picked up automatically.
function _include_isolated(relpath::AbstractString)
    name = Symbol(replace(relpath, r"[/\\.\-]" => "_"))
    mod = Module(name)
    Core.eval(mod, :(include(path::AbstractString) = Base.include($mod, path)))
    Core.eval(mod, :(eval(x) = Core.eval($mod, x)))
    Base.include(mod, joinpath(@__DIR__, relpath))
    return mod
end

function _run_isolated_test_script(relpath::AbstractString)
    script = joinpath(@__DIR__, relpath)
    cmd = `$(Base.julia_cmd()) --project=$(dirname(Base.active_project())) --startup-file=no --check-bounds=yes $script`
    return success(cmd)
end

# The test files under `dir` and its subdirectories, as paths relative to this
# directory, in order.
function _jl_files_recursive(dir)
    out = String[]
    for (root, dirs, files) in walkdir(joinpath(@__DIR__, dir))
        for f in sort(files)
            startswith(f, "test_") &&
                endswith(f, ".jl") &&
                push!(out, relpath(joinpath(root, f), @__DIR__))
        end
    end
    return sort(out)
end

# The directories of module suites, each run file by file in its own module.
# `verification` holds the verification cases (`domains`); `gates` holds
# measurement scripts, `quality` the package-hygiene environment and
# `type_analysis` a suite that runs in a process of its own.
const _SUITE_DIRS = filter(sort(readdir(@__DIR__))) do d
    isdir(joinpath(@__DIR__, d)) &&
        !(d in ("verification", "gates", "quality", "type_analysis"))
end

# `PYROLYSIS_TEST_GROUP` runs one part of the tests, so CI can run the parts
# side by side: `jacobian` (the Jacobian, multi-D and type-analysis suites),
# `burnout` (the burnout path's suites), `other` (every other module suite),
# `domains` (the verification cases), or `all`, the default.
const _GROUP = get(ENV, "PYROLYSIS_TEST_GROUP", "all")
const _GROUPS = ("all", "jacobian", "burnout", "other", "domains")
_GROUP in _GROUPS ||
    error("PYROLYSIS_TEST_GROUP must be one of $(join(_GROUPS, ", ")), not $(_GROUP)")

const _BURNOUT = (
    joinpath("closures", "test_correction.jl"),
    joinpath("closures", "test_layer_time.jl"),
    joinpath("closures", "test_mesh_motion.jl"),
    joinpath("deactivation", "test_deactivation.jl"),
    joinpath("integration", "test_policies.jl"),
    joinpath("kernels", "test_relative_flux.jl"),
)
function _group(f)
    any(d -> startswith(f, d * "/"), ("jacobian", "multid", "type_analysis")) &&
        return "jacobian"
    return f in _BURNOUT ? "burnout" : "other"
end
_runs(f) = _GROUP == "all" || _GROUP == _group(f)

# The module suites run both energy forms themselves, so a run in the
# enthalpy form (`PYROLYSIS_ENERGY_STATE=enthalpy`) is the verification cases'.
const _ENTHALPY_RUN = get(ENV, "PYROLYSIS_ENERGY_STATE", "temperature") == "enthalpy"

# Each file's wall time, so a CI partition can be balanced from a run's log.
function _timed(run, label)
    t = @elapsed run()
    println("test file ", label, ": ", round(t; digits = 1), " s")
    flush(stdout)
    return nothing
end

@testset "Pyrolysis.jl" begin
    if !_ENTHALPY_RUN
        for dir in _SUITE_DIRS, f in filter(_runs, _jl_files_recursive(dir))
            _timed(() -> _include_isolated(f), f)
        end
        # The type analysis of the entry points, in a fresh process.
        type_analysis = joinpath("type_analysis", "runtests.jl")
        if _runs(type_analysis)
            @testset "type analysis (isolated)" begin
                _timed(type_analysis) do
                    @test _run_isolated_test_script(type_analysis)
                end
            end
        end
    end

    # The verification cases, in one namespace (`runcases.jl` includes the
    # case definitions once).
    if _GROUP in ("all", "domains")
        f = joinpath("verification", "runcases.jl")
        _timed(() -> include(joinpath(@__DIR__, f)), f)
    end
end
