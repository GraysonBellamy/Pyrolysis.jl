using Test
using Pyrolysis

# Each file under jacobian/ is evaluated in its own throwaway module: the
# operator/structured tests define top-level structs and consts that would
# clash if dumped into one shared namespace. The module name is derived from
# the path so new files are picked up automatically.
function _include_isolated(relpath::AbstractString)
    name = Symbol(replace(relpath, r"[/\\.\-]" => "_"))
    mod = Module(name)
    Core.eval(mod, :(include(path::AbstractString) = Base.include($mod, path)))
    Core.eval(mod, :(eval(x) = Core.eval($mod, x)))
    Base.include(mod, joinpath(@__DIR__, relpath))
    return mod
end

function _test_project_for_subprocess()
    project = Base.active_project()
    project_dir = dirname(project)
    if occursin("BenchmarkTools", read(project, String))
        return project_dir
    end
    return @__DIR__
end

function _run_isolated_test_script(relpath::AbstractString)
    project_dir = _test_project_for_subprocess()
    script = joinpath(@__DIR__, relpath)
    cmd = `$(Base.julia_cmd()) --project=$project_dir --startup-file=no --check-bounds=yes $script`
    return success(cmd)
end

_jl_files(dir) = sort(filter(f -> endswith(f, ".jl"), readdir(joinpath(@__DIR__, dir))))

# Test files under `dir` and its subdirectories, as paths relative to this
# directory; `gates/` holds measurement scripts, not tests.
function _jl_files_recursive(dir)
    out = String[]
    for (root, dirs, files) in walkdir(joinpath(@__DIR__, dir))
        filter!(d -> d != "gates", dirs)
        for f in sort(files)
            endswith(f, ".jl") && push!(out, relpath(joinpath(root, f), @__DIR__))
        end
    end
    return sort(out)
end

# `PYROLYSIS_CORE=new` routes the verification suite through the
# finite-volume core (see verification/casedefs.jl) and runs only that core's
# tests and the suite.
const _NEW_CORE = get(ENV, "PYROLYSIS_CORE", "old") == "new"

# `PYROLYSIS_TEST_GROUP` runs one part of a run, so CI can run the parts side
# by side: `newcore` the finite-volume core suites, or one of their three
# groups (`newcore-jacobian`, `newcore-burnout`, `newcore-other`); `domains`
# everything else; `all`, the default, both.
const _GROUP = get(ENV, "PYROLYSIS_TEST_GROUP", "all")
const _NEWCORE_GROUPS = ("newcore-jacobian", "newcore-burnout", "newcore-other")
_GROUP in ("all", "newcore", "domains", _NEWCORE_GROUPS...) || error(
    "PYROLYSIS_TEST_GROUP must be all, newcore, $(join(_NEWCORE_GROUPS, ", ")) " *
    "or domains, not $(_GROUP)",
)
const _RUN_DOMAINS = _GROUP in ("all", "domains")

# The finite-volume core suites' groups, of similar runtime on the CI runner:
# the Jacobian and multi-D suites, the burnout path's suites, and the rest.
const _NEWCORE_BURNOUT = joinpath.(
    "newcore",
    (
        "test_correction.jl",
        "test_deactivation.jl",
        "test_mesh_motion.jl",
        "test_policies.jl",
        "test_relative_flux.jl",
    ),
)
function _newcore_group(f)
    any(d -> startswith(f, joinpath("newcore", d)), ("jacobian", "multid")) &&
        return "newcore-jacobian"
    return f in _NEWCORE_BURNOUT ? "newcore-burnout" : "newcore-other"
end
_run_newcore(f) = _GROUP in ("all", "newcore") || _GROUP == _newcore_group(f)

# Each file's wall time, so a CI partition can be balanced from a run's log.
function _timed(run, label)
    t = @elapsed run()
    println("test file ", label, ": ", round(t; digits = 1), " s")
    flush(stdout)
    return nothing
end

@testset "Pyrolysis.jl" begin
    # Domain suites — shared namespace; shared fixtures come from helpers.jl
    # (idempotent under repeated include via its own load guard).
    domains =
        !_RUN_DOMAINS ? String[] :
        _NEW_CORE ? ["verification"] :
        [
            "boundary_conditions",
            "discretization",
            "solver",
            "verification",
            "ale",
            "adaptivity",
            "radiation",
            "geometry",
            "materials",
            "physics",
            "residual",
        ]
    for dir in domains
        for f in _jl_files(dir)
            # The analytical tests solve on the old core, which the new-core
            # runs leave to the old-core run.
            _NEW_CORE && f == "test_analytical.jl" && continue
            _timed(() -> include(joinpath(@__DIR__, dir, f)), joinpath(dir, f))
        end
    end

    # Jacobian suites — one module each.
    if !_NEW_CORE && _RUN_DOMAINS
        for f in _jl_files("jacobian")
            _timed(
                () -> _include_isolated(joinpath("jacobian", f)),
                joinpath("jacobian", f),
            )
        end
    end

    # Finite-volume core suites — one module each. They run both energy forms
    # themselves, so the enthalpy-form run of the verification suite skips them.
    if !(_NEW_CORE && get(ENV, "PYROLYSIS_ENERGY_STATE", "temperature") == "enthalpy")
        for f in filter(_run_newcore, _jl_files_recursive("newcore"))
            _timed(() -> _include_isolated(f), f)
        end
    end

    # Performance invariants are process-state sensitive, especially the
    # allocation probes. Run the perf driver in a fresh Julia process so
    # `Pkg.test()` exercises the same hermetic path as the CI perf job.
    if !_NEW_CORE && _RUN_DOMAINS
        @testset "performance invariants (isolated)" begin
            _timed("perf/runtests.jl") do
                @test _run_isolated_test_script("perf/runtests.jl")
            end
        end
    end
end
