# Regenerate the per-chapter sources in docs/guide/ from the monolithic guides.
#
# The monoliths (docs/Technical_Reference_Guide.md, docs/User_Guide.md) are the
# canonical, hand-edited documents; the per-chapter files are derived views.
# Run this after editing a monolith:
#
#     julia scripts/split_guides.jl
#
# Chapters are delimited by top-level `# N. Title` headings outside fenced code
# blocks (see scripts/guide_split_core.jl). Each chapter is written to the
# existing file whose name starts with the zero-padded chapter number; a count
# mismatch aborts so a renamed or added chapter is resolved deliberately rather
# than silently misfiled.

include(joinpath(@__DIR__, "guide_split_core.jl"))

function regenerate(monolith::AbstractString, chapterdir::AbstractString)
    chapters = split_chapters(monolith)
    files = chapter_files(chapterdir)
    length(files) == length(chapters) || error(
        "$monolith has $(length(chapters)) chapters but $chapterdir has " *
        "$(length(files)) files — rename/add chapter files first, then rerun",
    )
    for (k, chunk) in enumerate(chapters)
        path = joinpath(chapterdir, files[k])
        write(path, join(chunk, "\n") * "\n")
        println("  wrote $(files[k]) ($(length(chunk)) lines)")
    end
end

root = normpath(joinpath(@__DIR__, ".."))
println("Technical Reference Guide:")
regenerate(
    joinpath(root, "docs", "Technical_Reference_Guide.md"),
    joinpath(root, "docs", "guide", "technical_reference"),
)
println("User Guide:")
regenerate(
    joinpath(root, "docs", "User_Guide.md"),
    joinpath(root, "docs", "guide", "user_guide"),
)
