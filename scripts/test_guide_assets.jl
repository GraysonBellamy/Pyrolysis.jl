# Focused generator regression checks; no solver or plotting dependencies.
# Run: julia scripts/test_guide_assets.jl
using Test
include("guide_split_core.jl")

@testset "canonical figure destinations" begin
    original = [
        "![mesh anatomy](assets/figures/01_mesh_anatomy.svg)",
        "[figure](assets/figures/01_mesh_anatomy.svg \"Mesh\")",
        "![external](https://example.org/assets/figures/x.svg)",
        "`![literal](assets/figures/example.svg)` and ![real](assets/figures/real.svg)",
        "````markdown",
        "```",
        "![example](assets/figures/example.svg)",
        "```",
        "````",
        "~~~markdown",
        "![example](assets/figures/example.svg)",
        "~~~",
        "![after](assets/figures/after.svg)",
    ]
    result = rebase_figure_paths(original, "../../assets/figures")
    @test result[1] == "![mesh anatomy](../../assets/figures/01_mesh_anatomy.svg)"
    @test result[2] == "[figure](../../assets/figures/01_mesh_anatomy.svg \"Mesh\")"
    @test result[3] == original[3]
    @test result[4] ==
          "`![literal](assets/figures/example.svg)` and ![real](../../assets/figures/real.svg)"
    @test result[5:12] == original[5:12]
    @test result[13] == "![after](../../assets/figures/after.svg)"
    @test rebase_figure_paths(result, "../../assets/figures") == result
end

@testset "Documenter asset copy and chapter destinations" begin
    root = normpath(joinpath(@__DIR__, "..", "docs"))
    run(`$(Base.julia_cmd()) $(joinpath(root, "generate_guide.jl"))`)
    source = joinpath(root, "assets", "figures")
    copied = joinpath(root, "src", "assets", "figures")
    exports = sort(filter(f -> endswith(f, ".svg"), readdir(source)))
    @test sort(readdir(copied)) == exports
    @test "08_quickstart_compact.svg" in exports
    for file in exports
        @test read(joinpath(source, file)) == read(joinpath(copied, file))
    end
    for name in ("user_guide", "technical_reference")
        dir = joinpath(root, "src", "guide", name)
        for file in chapter_files(dir)
            for m in
                eachmatch(r"!\[[^\]]+\]\(([^)]+\.svg)\)", read(joinpath(dir, file), String))
                @test isfile(normpath(joinpath(dir, m[1])))
            end
        end
    end
end

@testset "all guide image destinations exist" begin
    root = normpath(joinpath(@__DIR__, "..", "docs"))
    figures = Set{String}()
    for name in ("User_Guide.md", "Technical_Reference_Guide.md")
        for chapter in split_chapters(joinpath(root, name))
            canonical = join(chapter, "\n")
            derived = join(rebase_figure_paths(chapter, "../../assets/figures"), "\n")
            for m in eachmatch(r"!\[[^\]]+\]\((assets/figures/[^)]+)\)", canonical)
                push!(figures, m[1])
                @test isfile(joinpath(root, m[1]))
                @test occursin("](" * "../../" * m[1] * ")", derived)
            end
        end
    end
    @test length(figures) == 15
end
