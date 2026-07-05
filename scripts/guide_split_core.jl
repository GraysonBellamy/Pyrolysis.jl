# Shared chapter-splitting logic for the guide monoliths.
#
# The monoliths (docs/Technical_Reference_Guide.md, docs/User_Guide.md) are the
# canonical, hand-edited documents. Two consumers derive per-chapter views:
#
#   * scripts/split_guides.jl   — regenerates docs/guide/** for GitHub browsing
#   * docs/generate_guide.jl    — generates docs/src/guide/** for the Documenter
#                                 site (with Markdown-dialect conversion)
#
# Chapters are delimited by top-level `# N. Title` headings outside fenced code
# blocks. The preamble (document title and Contents list) precedes the first
# chapter heading and is deliberately dropped.

function split_chapters(monolith::AbstractString)
    lines = readlines(monolith)
    starts = Int[]
    in_fence = false
    for (i, line) in enumerate(lines)
        startswith(line, "```") && (in_fence = !in_fence)
        if !in_fence && occursin(r"^# \d+\. ", line)
            push!(starts, i)
        end
    end
    isempty(starts) && error("no chapter headings found in $monolith")
    chapters = Vector{Vector{String}}()
    for (k, s) in enumerate(starts)
        e = k < length(starts) ? starts[k+1] - 1 : length(lines)
        chunk = lines[s:e]
        # Strip the inter-chapter separator: trailing blanks and a lone `---`.
        while !isempty(chunk) && isempty(strip(chunk[end]))
            pop!(chunk)
        end
        if !isempty(chunk) && strip(chunk[end]) == "---"
            pop!(chunk)
            while !isempty(chunk) && isempty(strip(chunk[end]))
                pop!(chunk)
            end
        end
        push!(chapters, chunk)
    end
    return chapters
end

"Sorted `NN_slug.md` chapter filenames present in `chapterdir`."
function chapter_files(chapterdir::AbstractString)
    return sort(filter(f -> occursin(r"^\d{2}_.*\.md$", f), readdir(chapterdir)))
end
