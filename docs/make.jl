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

# Verification chapters are generated OFFLINE by scripts/verification/
# run_reports.jl into the COMMITTED docs/guide/verification/ (heavy solver
# ladders + figures; report mode of the verification suite plan). The docs
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
            "api/geometry.md",
            "api/problem.md",
            "api/solve.md",
            "api/residual.md",
            "api/jacobian.md",
            "api/sensitivity.md",
            "api/diagnostics.md",
            "api/adaptivity.md",
            "api/experimental.md",
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
