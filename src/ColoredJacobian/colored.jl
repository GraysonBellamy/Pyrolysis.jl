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
	ColoredJacobianCache(def, ws, u, t)

The residual's Jacobian by column-colored forward-mode automatic
differentiation on the connectivity pattern of [`jacobian_pattern`](@ref):
SparseMatrixColorings partitions the columns (greedy, natural order), and
DifferentiationInterface evaluates one dual-number residual per batch of
colors and decompresses into a `SparseMatrixCSC` whose structure is exactly
the pattern. Preparation happens once per problem structure; `jacobian!`
then reuses it at any state and time.

The decompression writes stored values in order and checks only the matrix
size and stored-entry count, so construction asserts that the prepared
pattern and [`jac_prototype`](@ref) share `colptr` and `rowval`; the
integrator's Jacobian must be allocated from that prototype.
"""
struct ColoredJacobianCache{F,B,P}
    f!::F
    backend::B
    prep::P
    pattern::SparseMatrixCSC{Bool,Int}
    du::Vector{Float64}
    n_colors::Int
end

function ColoredJacobianCache(
    def::CoreDefinition,
    ws::CoreWorkspace,
    u::AbstractVector,
    t::Real,
)
    pattern = jacobian_pattern(def.layout, def.setup.topology)
    f! = ResidualFunction(def, ws)
    backend = AutoSparse(
        AutoForwardDiff();
        sparsity_detector = KnownJacobianSparsityDetector(pattern),
        coloring_algorithm = GreedyColoringAlgorithm(),
    )
    du = zeros(length(u))
    x = Vector{Float64}(u)
    prep = prepare_jacobian(f!, du, backend, x, Constant(Float64(t)))
    prepared = sparsity_pattern(prep)
    (
        size(prepared) == size(pattern) &&
        SparseMatrixCSC(prepared).colptr == pattern.colptr &&
        SparseMatrixCSC(prepared).rowval == pattern.rowval
    ) || error("the prepared Jacobian's sparsity differs from the connectivity pattern")
    n_colors = length(column_groups(prep))
    return ColoredJacobianCache(f!, backend, prep, pattern, du, n_colors)
end

"""
	jac_prototype(cj::ColoredJacobianCache) -> SparseMatrixCSC{Float64,Int}

A zero matrix with exactly the colored Jacobian's structure.
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

"""
	same_structure(J, pattern) -> Bool

Whether `J` stores exactly the entries of `pattern`.
"""
same_structure(J::SparseMatrixCSC, pattern::SparseMatrixCSC) =
    size(J) == size(pattern) && J.colptr == pattern.colptr && J.rowval == pattern.rowval
