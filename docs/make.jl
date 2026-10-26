# Build the Pyrolysis.jl documentation site.
#
# Local build:
#     julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
#     julia --project=docs docs/make.jl
#
# The static site lands in docs/build/. Deployment to Cloudflare Pages is done
# by .github/workflows/docs.yml — there is deliberately no deploydocs() call;
# the site is a single unversioned build that tracks main.

using Documenter
using Pyrolysis

# Generates docs/src/guide/** from the canonical monoliths (also defines
# chapter_files via scripts/guide_split_core.jl).
include("generate_guide.jl")

guide_pages(name) =
    ["guide/$name/$f" for f in chapter_files(joinpath(@__DIR__, "src", "guide", name))]

# The documentation scan: no name the package no longer defines and no
# planning marker in what the site publishes (the guides, the hand-written
# pages, the verification report). The lists are the code scan's,
# test/code_scan/removed_names.txt.
function documentation_scan()
    lists = Dict("names" => Regex[], "markers" => Regex[])
    section = ""
    for line in eachline(joinpath(@__DIR__, "..", "test", "code_scan", "removed_names.txt"))
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
    files = vcat(
        joinpath.(@__DIR__, ["User_Guide.md", "Technical_Reference_Guide.md"]),
        [
            joinpath(@__DIR__, "src", f) for
            f in readdir(joinpath(@__DIR__, "src")) if endswith(f, ".md")
        ],
        [
            joinpath(@__DIR__, "src", "api", f) for
            f in readdir(joinpath(@__DIR__, "src", "api"))
        ],
        let d = joinpath(@__DIR__, "guide", "verification")
            isdir(d) ? [joinpath(d, f) for f in readdir(d) if endswith(f, ".md")] : String[]
        end,
    )
    hits = String[]
    for f in files,
        (n, line) in enumerate(eachline(f)),
        patterns in values(lists),
        p in patterns

        m = match(p, line)
        m === nothing || push!(hits, "$(relpath(f, @__DIR__)):$n: $(m.match)")
    end
    isempty(hits) ||
        error("documentation scan: removed names or planning markers\n" * join(hits, "\n"))
    return nothing
end

# The keyword check: the User Guide's `solve` keyword table names exactly the
# keywords `solve` accepts.
function keyword_check()
    lines = readlines(joinpath(@__DIR__, "User_Guide.md"))
    h = findfirst(
        l -> startswith(l, "#") && occursin("solve keywords", lowercase(l)),
        lines,
    )
    h === nothing && error("keyword check: no heading naming the `solve` keywords")
    documented = Symbol[]
    for l in lines[(h+1):end]
        startswith(l, "#") && break
        m = match(r"^\|\s*`([A-Za-z_][A-Za-z0-9_]*)`\s*\|", l)
        m === nothing || push!(documented, Symbol(m[1]))
    end
    accepted = collect(Pyrolysis.Integration._SOLVE_KEYWORDS)
    missing_ = setdiff(accepted, documented)
    extra = setdiff(documented, accepted)
    isempty(missing_) && isempty(extra) && allunique(documented) ||
        error("keyword check: undocumented $(missing_), not accepted $(extra)")
    return nothing
end

documentation_scan()
keyword_check()

# Verification chapters are generated OFFLINE by scripts/verification/
# run_reports.jl into the COMMITTED docs/guide/verification/ (heavy solver
# ladders + figures; the verification suite's report mode). The docs
# build only copies them into the gitignored Documenter tree.
let src = joinpath(@__DIR__, "guide", "verification"),
    dst = joinpath(@__DIR__, "src", "guide", "verification")

    if isdir(src)
        rm(dst; force = true, recursive = true)
        cp(src, dst)
    end
end

# index.md first, then one chapter per case family (V-HC-01.md, ...); the
# NN_slug.md convention of chapter_files doesn't apply here.
function verification_pages()
    dir = joinpath(@__DIR__, "src", "guide", "verification")
    isdir(dir) || return String[]
    chapters = sort(filter(f -> endswith(f, ".md") && f != "index.md", readdir(dir)))
    return ["guide/verification/$f" for f in vcat("index.md", chapters)]
end

pages = Any[
    "Home"=>"index.md",
    "User Guide"=>guide_pages("user_guide"),
    "Technical Reference"=>guide_pages("technical_reference"),
]
let vpages = verification_pages()
    isempty(vpages) || push!(pages, "Verification" => vpages)
end
append!(
    pages,
    Any[
        "Examples"=>"examples.md",
        "API Reference"=>[
            "api/index.md",
            "api/materials.md",
            "api/properties.md",
            "api/reactions.md",
            "api/boundary_conditions.md",
            "api/radiation.md",
            "api/meshes.md",
            "api/problem.md",
            "api/solve.md",
            "api/solution.md",
            "api/states.md",
            "api/kernels.md",
            "api/deactivation.md",
            "api/jacobian.md",
            "api/extensions.md",
            "api/internal.md",
        ],
        "Development"=>"development.md",
    ],
)

makedocs(
    sitename = "Pyrolysis.jl",
    authors = "Grayson Bellamy",
    modules = [Pyrolysis],
    checkdocs = :exports,
    doctest = false,   # docstring examples are ```julia blocks, not jldoctests (yet)
    warnonly = [:cross_references, :missing_docs],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://pyrolysis.graysonbellamy.dev",
        edit_link = nothing,   # most pages are generated; no useful edit target
        size_threshold = 2_000_000,        # guide chapters are large
        size_threshold_warn = 1_000_000,
    ),
    pages = pages,
)
