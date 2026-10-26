# Small dependency-free SVG drawing vocabulary shared by the guide figures.
# Coordinates are in SVG pixels; physical coordinates/units are labelled explicitly.
using Printf

const INK = "#202020"
const MUTED = "#555555"
const BLUE = "#295674"
const TEAL = "#365b58"
const ORANGE = "#92512f"
const PURPLE = "#675878"
const PALE = "#f7f7f7"
const GREEN = "#ededed"
const GREY = "#d9d9d9"
const OUT = normpath(joinpath(@__DIR__, "..", "..", "docs", "assets", "figures"))

xml(s) = replace(string(s), '&' => "&amp;", '<' => "&lt;", '>' => "&gt;", '"' => "&quot;")
num(x) = @sprintf("%.3f", x)

# Explicit subscripts/superscripts keep mathematical notation distinct from
# case-sensitive API identifiers (the latter are enclosed in backticks).
function typeset(s)
    text = replace(
        xml(s),
        r"[_^]\{[^{}]+\}" =>
            m -> begin
                shift = startswith(m, "_") ? "sub" : "super"
                body = m[3:prevind(m, lastindex(m))]
                "<tspan baseline-shift=\"$shift\" font-size=\"75%\">$body</tspan>"
            end,
    )
    return replace(
        text,
        r"`[^`]+`" =>
            m ->
                "<tspan font-family=\"Liberation Mono, monospace\" font-size=\"90%\">" *
                m[2:prevind(m, lastindex(m))] *
                "</tspan>",
    )
end

function svg(f, name, title, description; width = 1200, height = 600)
    mkpath(OUT)
    io = IOBuffer()
    println(
        io,
        """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height" role="img" aria-labelledby="title desc">
<title id="title">$(xml(title))</title><desc id="desc">$(xml(description))</desc>
<defs>""",
    )
    for (name, color) in (
        ("ink", INK),
        ("blue", BLUE),
        ("teal", TEAL),
        ("orange", ORANGE),
        ("purple", PURPLE),
    )
        println(
            io,
            """<marker id="$name" viewBox="0 0 10 10" refX="9" refY="5" markerUnits="userSpaceOnUse" markerWidth="11" markerHeight="11" orient="auto-start-reverse"><path d="M 1 1 L 9 5 L 1 9" fill="none" stroke="$color" stroke-width="1.3"/></marker>""",
        )
    end
    println(
        io,
        """</defs><rect width="100%" height="100%" fill="white"/><g font-family="Liberation Sans, Arial, sans-serif" fill="$INK">""",
    )
    label(io, 32, 39, title; size = 24, weight = "bold")
    f(io)
    println(io, "</g></svg>")
    write(joinpath(OUT, name * ".svg"), String(take!(io)))
    println("  wrote ", name, ".svg")
end

function label(
    io,
    x,
    y,
    s;
    size = 20,
    color = INK,
    anchor = "start",
    weight = "normal",
    rotate = nothing,
)
    transform = rotate === nothing ? "" : "transform=\"rotate($rotate $x $y)\""
    println(
        io,
        """<text x="$(num(x))" y="$(num(y))" font-size="$size" fill="$color" text-anchor="$anchor" font-weight="$weight" $transform>$(typeset(s))</text>""",
    )
end
function lines(io, x, y, items; step = 27, kw...)
    for (k, s) in enumerate(items)
        label(io, x, y + (k - 1) * step, s; kw...)
    end
end
function rect(io, x, y, w, h; fill = PALE, stroke = INK, radius = 0, dash = "")
    println(
        io,
        """<rect x="$x" y="$y" width="$w" height="$h" rx="$radius" fill="$fill" stroke="$stroke" stroke-width="1.2" stroke-dasharray="$dash"/>""",
    )
end
function line(
    io,
    x1,
    y1,
    x2,
    y2;
    color = INK,
    width = 1.5,
    dash = "",
    arrow = nothing,
    both = false,
)
    markers =
        arrow === nothing ? "" :
        "marker-end=\"url(#$arrow)\"" * (both ? " marker-start=\"url(#$arrow)\"" : "")
    println(
        io,
        """<line x1="$x1" y1="$y1" x2="$x2" y2="$y2" stroke="$color" stroke-width="$width" stroke-dasharray="$dash" $markers/>""",
    )
end
function circle(io, x, y, r = 5; fill = INK, stroke = fill)
    println(
        io,
        """<circle cx="$x" cy="$y" r="$r" fill="$fill" stroke="$stroke" stroke-width="1.3"/>""",
    )
end
function path(io, d; color = BLUE, fill = "none", width = 1.5, dash = "")
    println(
        io,
        """<path d="$d" fill="$fill" stroke="$color" stroke-width="$width" stroke-dasharray="$dash"/>""",
    )
end
function polyline(io, xs, ys; color = BLUE, width = 2.2, dash = "")
    # Break at NaN (e.g. dead cells); do not bridge a missing interval.
    pen = false
    d = IOBuffer()
    for (x, y) in zip(xs, ys)
        if isfinite(x) && isfinite(y)
            print(d, pen ? " L " : " M ", num(x), " ", num(y))
            pen = true
        else
            pen = false
        end
    end
    path(io, String(take!(d)); color, width, dash)
end
function panel(io, x, y, w, h, title)
    line(io, x, y, x + w, y; color = MUTED, width = 1)
    label(io, x + 18, y + 30, title; weight = "bold", size = 20)
end
function box(io, x, y, w, h, title, body; fill = "white", color = INK)
    rect(io, x, y, w, h; fill, stroke = color)
    label(io, x + 16, y + 30, title; weight = "bold", size = 20)
    lines(io, x + 16, y + 59, body; size = 18, step = 25)
end
function dimension(io, x1, x2, y, text; color = MUTED)
    line(io, x1, y, x2, y; color)
    for x in (x1, x2)
        line(io, x, y - 6, x, y + 6; color)
    end
    label(io, (x1 + x2) / 2, y + 26, text; color, anchor = "middle", size = 18)
end
function meshstrip(
    io,
    nodes,
    y;
    h = 65,
    fills = fill(PALE, length(nodes) - 1),
    ids = true,
    centers = false,
)
    for i = 1:(length(nodes)-1)
        rect(io, nodes[i], y, nodes[i+1] - nodes[i], h; fill = fills[i], stroke = INK)
        ids && label(
            io,
            (nodes[i] + nodes[i+1]) / 2,
            y + h / 2 + 7,
            string(i);
            anchor = "middle",
            size = 18,
        )
        centers && circle(io, (nodes[i] + nodes[i+1]) / 2, y + h / 2; fill = TEAL)
    end
    for x in nodes
        circle(io, x, y + h, 3.5; fill = INK)
    end
end

# Draw a chart inside a specified rectangle, returning physical-to-SVG maps.
function chart_axes(io, x, y, w, h; xlim, ylim, xticks, yticks, xlabel, ylabel)
    X(v) = x + w * (v - xlim[1]) / (xlim[2] - xlim[1])
    Y(v) = y + h - h * (v - ylim[1]) / (ylim[2] - ylim[1])
    for (v, s) in yticks
        line(io, x, Y(v), x + w, Y(v); color = "#e5e5e5", width = 0.7)
        label(io, x - 12, Y(v) + 6, s; anchor = "end", size = 17, color = MUTED)
    end
    for (v, s) in xticks
        line(io, X(v), y + h, X(v), y + h + 6; color = MUTED)
        label(io, X(v), y + h + 27, s; anchor = "middle", size = 17, color = MUTED)
    end
    line(io, x, y, x, y + h; color = MUTED)
    line(io, x, y + h, x + w, y + h; color = MUTED)
    label(io, x + w / 2, y + h + 56, xlabel; anchor = "middle", size = 19)
    label(io, x - 56, y + h / 2, ylabel; anchor = "middle", size = 19, rotate = -90)
    return X, Y
end
function key(io, x, y, s, color; dash = "")
    line(io, x, y - 6, x + 32, y - 6; color, width = 2.2, dash)
    label(io, x + 42, y, s; size = 18)
end
