# Generate the Documenter guide sources (docs/src/guide/**) from the canonical
# monolithic guides. Chapters are split with the same logic as
# scripts/split_guides.jl (shared via scripts/guide_split_core.jl) and then
# converted from GitHub-flavored to Documenter-flavored Markdown:
#
#   * `$$` display-math blocks  ->  ```math fenced blocks
#   * `$...$` inline math       ->  ``...`` inline math (pairs may span a
#     hard-wrapped line break within a paragraph; the match is joined onto
#     one line, with blockquote `>` markers handled)
#
# Fenced code blocks and inline code spans are left untouched. Any `$` or
# `](#...)` anchor link remaining outside code after conversion aborts the
# build: unconverted math would render as mangled prose, and no intra-document
# anchors are expected to survive chapter splitting.
#
# Runs standalone (`julia docs/generate_guide.jl`) or via `docs/make.jl`.

include(joinpath(@__DIR__, "..", "scripts", "guide_split_core.jl"))

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

# --- fence tracking ---------------------------------------------------------

# Tracks fenced code blocks, honoring CommonMark's rule that a fence opened
# with N backticks is only closed by a bare run of >= N backticks (so ````-
# wrapped Markdown examples containing ``` lines are handled correctly).
mutable struct FenceTracker
    ticks::Int   # 0 = outside a fence
end

"Advance the tracker with `line`; return true if the line is a fence delimiter."
function fence_line!(ft::FenceTracker, line::AbstractString)
    m = match(r"^\s{0,3}(`{3,})", line)
    m === nothing && return false
    if ft.ticks == 0
        ft.ticks = length(m.captures[1])
        return true
    elseif length(m.captures[1]) >= ft.ticks && occursin(r"^\s{0,3}`{3,}\s*$", line)
        ft.ticks = 0
        return true
    end
    return false   # backtick run inside a longer-delimited fence: plain content
end

# --- inline-math conversion (paragraph level) -------------------------------

const _CODE_SPAN_RE = r"(`+)(.*?)\1"

# `$...$` where the content has no `$`, and any newline is a soft wrap (not a
# blank line / paragraph break).
const _INLINE_MATH_RE = r"\$((?:[^\$\n]|\n(?!\s*\n))+?)\$"

# Convert inline math in a multi-line text chunk (no code fences, no display
# math, no blockquote markers). Inline code spans are protected first so `$`
# inside code survives; math pairs that would swallow a code span are left
# untouched for the validator to flag.
function convert_inline(text::String)
    saved = String[]
    prot = replace(text, _CODE_SPAN_RE => s -> begin
        push!(saved, s)
        "\x01$(length(saved))\x02"
    end)
    prot = replace(prot, _INLINE_MATH_RE => s -> begin
        content = s[2:prevind(s, lastindex(s))]
        if isempty(strip(content)) || occursin('\x01', content)
            s   # degenerate or code-span-crossing match: leave for validation
        else
            "``" * replace(content, r"\s*\n\s*" => " ") * "``"
        end
    end)
    return replace(prot, r"\x01(\d+)\x02" => s -> saved[parse(Int, s[2:prevind(s, lastindex(s))])])
end

const _QUOTE_RE = r"^\s{0,3}>\s?"

# Convert a run of consecutive plain-text lines, handling blockquote paragraphs
# by stripping the `>` markers, converting, and re-prefixing (so math wrapped
# across quoted lines doesn't capture the `>` marker as an operator). A quote
# paragraph only starts at a paragraph boundary; a `>`-leading line directly
# under a text line is hard-wrapped prose/math (e.g. `$N\n> 1$`), not a quote.
function convert_text_block(blk::Vector{String})
    res = String[]
    i = 1
    while i <= length(blk)
        quote_start(k) = occursin(_QUOTE_RE, blk[k]) && (k == 1 || isempty(strip(blk[k-1])))
        local j
        if quote_start(i)
            j = i
            while j <= length(blk) && occursin(_QUOTE_RE, blk[j])
                j += 1
            end
            inner = join([replace(l, _QUOTE_RE => "") for l in blk[i:(j-1)]], "\n")
            for l in split(convert_inline(inner), "\n")
                push!(res, isempty(l) ? ">" : "> " * l)
            end
        else
            j = i + 1
            while j <= length(blk) && !quote_start(j)
                j += 1
            end
            append!(res, split(convert_inline(join(blk[i:(j-1)], "\n")), "\n"))
        end
        i = j
    end
    return res
end

# --- chapter transform ------------------------------------------------------

function transform_chapter(lines::Vector{String}, label::AbstractString)
    # Pass 1 (line-based): convert display math, classify lines.
    out = String[]
    text = Bool[]
    ft = FenceTracker(0)
    in_math = false
    math_indent = ""
    for line in lines
        if in_math
            # Closer: bare `$$` or trailing content followed by `$$`.
            m = match(r"^(\s*)(.*?)\$\$\s*$", line)
            if m === nothing
                push!(out, line)
                push!(text, false)
            else
                if !isempty(m.captures[2])
                    push!(out, m.captures[1] * m.captures[2])
                    push!(text, false)
                end
                push!(out, math_indent * "```")
                push!(text, false)
                in_math = false
            end
            continue
        end
        if fence_line!(ft, line) || ft.ticks != 0
            push!(out, line)
            push!(text, false)
            continue
        end
        m_single = match(r"^(\s*)\$\$(.+?)\$\$\s*$", line)
        m_open = match(r"^(\s*)\$\$(.*)$", line)   # bare `$$` or `$$content`
        if m_single !== nothing
            ind = m_single.captures[1]
            append!(out, [ind * "```math", ind * strip(m_single.captures[2]), ind * "```"])
            append!(text, [false, false, false])
        elseif m_open !== nothing
            math_indent = m_open.captures[1]
            push!(out, math_indent * "```math")
            push!(text, false)
            if !isempty(strip(m_open.captures[2]))
                push!(out, math_indent * strip(m_open.captures[2]))
                push!(text, false)
            end
            in_math = true
        else
            push!(out, line)
            push!(text, true)
        end
    end
    in_math && error("$label: unterminated \$\$ display-math block")

    # Pass 2: convert inline math over runs of consecutive text lines.
    final = String[]
    i = 1
    while i <= length(out)
        if !text[i]
            push!(final, out[i])
            i += 1
            continue
        end
        j = i
        while j <= length(out) && text[j]
            j += 1
        end
        append!(final, convert_text_block(out[i:(j-1)]))
        i = j
    end
    return final
end

# --- validation -------------------------------------------------------------

function validate_chapter!(problems::Vector{String}, lines::Vector{String}, label::AbstractString)
    ft = FenceTracker(0)
    for (i, line) in enumerate(lines)
        (fence_line!(ft, line) || ft.ticks != 0) && continue
        bare = replace(line, _CODE_SPAN_RE => "")
        occursin('$', bare) && push!(problems, "$label:$i: unconverted \$: $line")
        occursin("](#", bare) &&
            push!(problems, "$label:$i: intra-document anchor link: $line")
    end
end

# --- driver -----------------------------------------------------------------

function generate_guide_sources(; verbose::Bool = true)
    specs = (
        (monolith = "User_Guide.md", name = "user_guide"),
        (monolith = "Technical_Reference_Guide.md", name = "technical_reference"),
    )
    problems = String[]
    for spec in specs
        monolith = joinpath(REPO_ROOT, "docs", spec.monolith)
        registry = joinpath(REPO_ROOT, "docs", "guide", spec.name)
        outdir = joinpath(@__DIR__, "src", "guide", spec.name)
        chapters = split_chapters(monolith)
        files = chapter_files(registry)
        length(files) == length(chapters) || error(
            "$monolith has $(length(chapters)) chapters but $registry has " *
            "$(length(files)) files — run scripts/split_guides.jl bookkeeping first",
        )
        rm(outdir; force = true, recursive = true)
        mkpath(outdir)
        for (k, chunk) in enumerate(chapters)
            label = "guide/$(spec.name)/$(files[k])"
            transformed = transform_chapter(chunk, label)
            validate_chapter!(problems, transformed, label)
            write(joinpath(outdir, files[k]), join(transformed, "\n") * "\n")
        end
        verbose && println("generated $(length(chapters)) pages in docs/src/guide/$(spec.name)/")
    end
    if !isempty(problems)
        error(
            "guide conversion left unconverted Markdown " *
            "($(length(problems)) findings):\n" * join(problems, "\n"),
        )
    end
    return nothing
end

generate_guide_sources()
