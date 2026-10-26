"""
No source, test, example, script, extension, benchmark or workflow file names a
definition the package no longer has, or cites the planning records' phases,
gate items and plan sections: a renamed or deleted name left in a comment, a
docstring or a script that is not run would otherwise go unnoticed. The names
and markers are `removed_names.txt`.
"""

using Test

const _ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const _SCOPE = ("src", "test", "examples", "scripts", "ext", "benchmark", ".github")
const _EXTENSIONS = (".jl", ".md", ".toml", ".yml", ".yaml")
const _SELF = (
    joinpath("test", "code_scan", "test_removed_names.jl"),
    joinpath("test", "code_scan", "removed_names.txt"),
)

# The two lists of `removed_names.txt` as regular expressions.
function _patterns()
    lists = Dict("names" => Regex[], "markers" => Regex[])
    section = ""
    for line in eachline(joinpath(@__DIR__, "removed_names.txt"))
        line = strip(line)
        (isempty(line) || startswith(line, "#")) && continue
        if startswith(line, "[")
            section = strip(line, ['[', ']'])
            continue
        end
        body =
            startswith(line, "re:") ? line[4:end] :
            replace(line, r"[\\^$.|?*+()\[\]{}]" => s"\\\0")
        section == "names" && (body = "(?<![A-Za-z0-9_!])" * body * "(?![A-Za-z0-9_!])")
        push!(lists[section], Regex(body))
    end
    return lists
end

# The files in scope: in a git checkout the tracked files and the untracked
# ones git does not ignore (the local records are ignored); otherwise every
# file under the scope's directories. Manifests list dependencies, not names.
function _files()
    rel = if ispath(joinpath(_ROOT, ".git")) && Sys.which("git") !== nothing
        readlines(`git -C $_ROOT ls-files --cached --others --exclude-standard`)
    else
        [
            relpath(joinpath(r, f), _ROOT) for d in _SCOPE if isdir(joinpath(_ROOT, d))
            for (r, _, fs) in walkdir(joinpath(_ROOT, d)) for f in fs
        ]
    end
    return filter(rel) do f
        any(d -> startswith(f, d * "/"), _SCOPE) &&
            any(e -> endswith(f, e), _EXTENSIONS) &&
            basename(f) != "Manifest.toml" &&
            !(f in _SELF) &&
            isfile(joinpath(_ROOT, f))
    end
end

@testset "no removed name or planning marker in code" begin
    lists = _patterns()
    files = _files()
    @test length(files) > 100
    for (list, patterns) in lists
        hits = String[]
        for f in files, (n, line) in enumerate(eachline(joinpath(_ROOT, f)))
            for p in patterns
                m = match(p, line)
                m === nothing || push!(hits, "$f:$n: $(m.match)")
            end
        end
        isempty(hits) || @info "removed $list found" hits
        @test isempty(hits)
    end
end
