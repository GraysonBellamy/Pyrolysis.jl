# ============================================================================
# Colored sparse forward-mode Jacobian
# ============================================================================

"""
	ResidualFunction(def, ws)

The residual as `f!(du, u, t)`, the form the differentiation backend calls.
"""
struct ResidualFunction{D<:CoreDefinition,W<:CoreWorkspace}
    def::D
    ws::W
end
(f::ResidualFunction)(du, u, t) = residual!(du, u, f.def, t, f.ws)

"""
	LateralRows

The two dense rows of a variable cross section, assembled after the colored
pass. Every dependence of the area and lateral-rate rows on the state passes
through six column sums `s = (Σ G, Σ M, Σ Ġ, Σ Ṁ, θ_A, A)`
(`Closures.lateral_rows`), so each entry is `Σₖ (∂r/∂sₖ)(∂sₖ/∂u)`: `∂r/∂s`
by forward differentiation of the six-number function, `∂(Σ Ġ)/∂u` and
`∂(Σ Ṁ)/∂u` the column sums of the history rows of the colored Jacobian,
which the colored pass gets exactly, and the rest unit vectors.

`slots` lists, per dense row, its stored columns and their positions in the
Jacobian's stored values; `member` marks each history row (and column) as
`G` (1) or `M` (2), every other row 0; `colsum` holds the two column sums
and `du` a residual buffer.
"""
struct LateralRows{X}
    lateral::X
    slots::NTuple{2,Vector{Tuple{Int,Int}}}
    member::Vector{Int8}
    colsum::Matrix{Float64}
    du::Vector{Float64}
end

function LateralRows(def::CoreDefinition, pattern::SparseMatrixCSC{Bool,Int})
    layout = def.layout
    rows = (area_index(layout), lateral_rate_index(layout))
    slots = (Tuple{Int,Int}[], Tuple{Int,Int}[])
    for col in axes(pattern, 2), k in nzrange(pattern, col)
        r = pattern.rowval[k]
        r == rows[1] && push!(slots[1], (col, k))
        r == rows[2] && push!(slots[2], (col, k))
    end
    L = state_length(layout)
    member = zeros(Int8, L)
    for i = 1:n_cells(layout), k = 1:2
        member[history_index(layout, i, k)] = k
    end
    return LateralRows(def.closure.lateral, slots, member, zeros(2, L), zeros(L))
end

"""
	COLORING_ORDERS

The vertex orders the greedy column coloring tries, keeping the one with the
fewest colors. On a moving mesh the natural order takes one color more than
the largest row degree, the lower bound; smallest-last and incidence-degree
orders reach it on every configuration of the test suite.
"""
const COLORING_ORDERS = (NaturalOrder(), LargestFirst(), SmallestLast(), IncidenceDegree())

"""
	MAX_CHUNK

The widest batch of colors one dual-number residual carries. The colors are
split into the fewest batches of at most this width, of equal width, so up
to `MAX_CHUNK` colors cost one pass.
"""
const MAX_CHUNK = 24

# The batch width for `K` colors.
_chunk_width(K::Int) = cld(K, cld(K, MAX_CHUNK))

# The pattern the coloring sees: the dense rows reduced to their diagonal.
function _colored_pattern(pattern::SparseMatrixCSC{Bool,Int}, dense::Vector{Int})
    colored = copy(pattern)
    for r in dense, c in findall(colored[r, :])
        c == r || (colored[r, c] = false)
    end
    return dropzeros!(colored)
end

_color_count(colored::SparseMatrixCSC{Bool,Int}) =
    ncolors(coloring(colored, ColoringProblem(), GreedyColoringAlgorithm(COLORING_ORDERS)))

"""
	color_count(def) -> Int

The number of colors of the residual's Jacobian outside its dense rows.
"""
color_count(def::CoreDefinition) =
    _color_count(_colored_pattern(jacobian_pattern(def), dense_rows(def.layout)))

"""
	ColoredJacobianCache(def, ws, u, t)

The residual's Jacobian by column-colored forward-mode automatic
differentiation on the connectivity pattern of [`jacobian_pattern`](@ref):
SparseMatrixColorings partitions the columns greedily in the best of
[`COLORING_ORDERS`](@ref), and DifferentiationInterface evaluates one
dual-number residual per batch of at most [`MAX_CHUNK`](@ref) colors and
decompresses into a `SparseMatrixCSC`. Preparation happens once
per problem structure; `jacobian!` then reuses it at any state and time.

The [`dense_rows`](@ref) of a variable cross section are colored by their
diagonal alone, the colored result is scattered into the full pattern, and
[`LateralRows`](@ref) writes them afterwards.

The decompression writes stored values in order and checks only the matrix
size and stored-entry count, so construction asserts that the prepared
pattern and the colored pattern share `colptr` and `rowval`; the
integrator's Jacobian must be allocated from [`jac_prototype`](@ref).
"""
struct ColoredJacobianCache{F,B,P,D}
    f!::F
    backend::B
    prep::P
    pattern::SparseMatrixCSC{Bool,Int}
    colored::SparseMatrixCSC{Bool,Int}
    J_colored::SparseMatrixCSC{Float64,Int}
    scatter::Vector{Int}
    dense::D
    du::Vector{Float64}
    n_colors::Int
end

function ColoredJacobianCache(
    def::CoreDefinition,
    ws::CoreWorkspace,
    u::AbstractVector,
    t::Real,
)
    pattern = jacobian_pattern(def)
    dense = dense_rows(def.layout)
    colored = _colored_pattern(pattern, dense)
    f! = ResidualFunction(def, ws)
    backend = AutoSparse(
        AutoForwardDiff(; chunksize = _chunk_width(_color_count(colored)));
        sparsity_detector = KnownJacobianSparsityDetector(colored),
        coloring_algorithm = GreedyColoringAlgorithm(COLORING_ORDERS),
    )
    du = zeros(length(u))
    x = Vector{Float64}(u)
    prep = prepare_jacobian(f!, du, backend, x, Constant(Float64(t)))
    prepared = sparsity_pattern(prep)
    (
        size(prepared) == size(colored) &&
        SparseMatrixCSC(prepared).colptr == colored.colptr &&
        SparseMatrixCSC(prepared).rowval == colored.rowval
    ) || error("the prepared Jacobian's sparsity differs from the connectivity pattern")
    scatter = Int[]
    for col in axes(colored, 2), k in nzrange(colored, col)
        r = colored.rowval[k]
        range = nzrange(pattern, col)
        push!(scatter, range[searchsortedfirst(view(pattern.rowval, range), r)])
    end
    J_colored = SparseMatrixCSC(
        colored.m,
        colored.n,
        copy(colored.colptr),
        copy(colored.rowval),
        zeros(nnz(colored)),
    )
    lateral = isempty(dense) ? nothing : LateralRows(def, pattern)
    n_colors = length(column_groups(prep))
    return ColoredJacobianCache(
        f!,
        backend,
        prep,
        pattern,
        colored,
        J_colored,
        scatter,
        lateral,
        du,
        n_colors,
    )
end

"""
	jac_prototype(cj::ColoredJacobianCache) -> SparseMatrixCSC{Float64,Int}

A zero matrix with exactly the Jacobian's structure.
"""
function jac_prototype(cj::ColoredJacobianCache)
    p = cj.pattern
    return SparseMatrixCSC(p.m, p.n, copy(p.colptr), copy(p.rowval), zeros(nnz(p)))
end

"""
	jacobian!(J, cj, u, t) -> J

Write the Jacobian at `(u, t)` into `J`, whose structure must be the pattern's.
"""
function jacobian!(
    J::SparseMatrixCSC{Float64,Int},
    cj::ColoredJacobianCache,
    u::AbstractVector,
    t::Real,
)
    _jacobian!(J, cj, cj.dense, u, t)
    return J
end

function _jacobian!(J, cj::ColoredJacobianCache, ::Nothing, u, t)
    DifferentiationInterface.jacobian!(
        cj.f!,
        cj.du,
        J,
        cj.prep,
        cj.backend,
        u,
        Constant(Float64(t)),
    )
    return J
end

function _jacobian!(J, cj::ColoredJacobianCache, dense::LateralRows, u, t)
    Jc = cj.J_colored
    DifferentiationInterface.jacobian!(
        cj.f!,
        cj.du,
        Jc,
        cj.prep,
        cj.backend,
        u,
        Constant(Float64(t)),
    )
    fill!(J.nzval, 0.0)
    @inbounds for k in eachindex(cj.scatter)
        J.nzval[cj.scatter[k]] = Jc.nzval[k]
    end
    _lateral_rows!(J, Jc, dense, cj.f!.def, cj.f!.ws, u, t)
    return J
end

# Write the two dense rows from the column sums of the history rows.
function _lateral_rows!(J, Jc, d::LateralRows, def::CoreDefinition, ws, u, t)
    layout = def.layout
    residual!(d.du, u, def, t, ws)
    ΣG, ΣM, ΣĠ, ΣṀ = 0.0, 0.0, 0.0, 0.0
    @inbounds for i = 1:n_cells(layout)
        kG, kM = history_index(layout, i, 1), history_index(layout, i, 2)
        layout.mask[kG] || continue
        ΣG += u[kG]
        ΣM += u[kM]
        ΣĠ += d.du[kG]
        ΣṀ += d.du[kM]
    end
    s = SVector(ΣG, ΣM, ΣĠ, ΣṀ, u[lateral_rate_index(layout)], u[area_index(layout)])
    ∂r = ForwardDiff.jacobian(x -> lateral_rows(d.lateral, x), s)
    fill!(d.colsum, 0.0)
    @inbounds for col in axes(Jc, 2), k in nzrange(Jc, col)
        r = Jc.rowval[k]
        side = d.member[r]
        (side == 0 || !layout.mask[r]) && continue
        d.colsum[side, col] += Jc.nzval[k]
    end
    A, θ_A = area_index(layout), lateral_rate_index(layout)
    @inbounds for (row, slots) in enumerate(d.slots), (col, k) in slots
        side = layout.mask[col] ? d.member[col] : Int8(0)
        J.nzval[k] =
            (side == 1 ? ∂r[row, 1] : 0.0) +
            (side == 2 ? ∂r[row, 2] : 0.0) +
            ∂r[row, 3] * d.colsum[1, col] +
            ∂r[row, 4] * d.colsum[2, col] +
            (col == θ_A ? ∂r[row, 5] : 0.0) +
            (col == A ? ∂r[row, 6] : 0.0)
    end
    return J
end

"""
	same_structure(J, pattern) -> Bool

Whether `J` stores exactly the entries of `pattern`.
"""
same_structure(J::SparseMatrixCSC, pattern::SparseMatrixCSC) =
    size(J) == size(pattern) && J.colptr == pattern.colptr && J.rowval == pattern.rowval
