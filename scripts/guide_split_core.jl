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

"""
    rebase_figure_paths(lines, prefix)

Rebase canonical `assets/figures/` Markdown destinations for a derived view.
Only standalone image/link destinations outside fenced code are rewritten;
code examples and external URLs stay unchanged. `prefix` is the relative
path to that view's figure directory, with no trailing slash.
"""
function rebase_figure_paths(lines::Vector{String}, prefix::AbstractString)
    out = String[]
    fence = 0
    for line in lines
        m = match(r"^\s{0,3}(`{3,}|~{3,})", line)
        if m !== nothing
            n = length(m[1])
            # Positive = backticks, negative = tildes. A closing fence must
            # use the same character and contain no info string.
            kind = startswith(m[1], "`") ? 1 : -1
            if fence == 0
                fence = kind * n
            elseif sign(fence) == kind &&
                   n >= abs(fence) &&
                   occursin(r"^\s{0,3}(`{3,}|~{3,})\s*$", line)
                fence = 0
            end
            push!(out, line)
        elseif fence != 0
            push!(out, line)
        else
            # Do not touch literal Markdown links inside inline code spans.
            pieces = split(line, r"(`+[^`]*`+)", keepempty = true)
            spans = collect(eachmatch(r"`+[^`]*`+", line))
            rebased = IOBuffer()
            for (i, piece) in enumerate(pieces)
                print(rebased, replace(piece, r"\]\(assets/figures/" => "]($prefix/"))
                i <= length(spans) && print(rebased, spans[i].match)
            end
            push!(out, String(take!(rebased)))
        end
    end
    return out
end
