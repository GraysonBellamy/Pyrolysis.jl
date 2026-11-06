# =============================================================================
# Gauss–Legendre quadrature
# =============================================================================
#
# Nodes and weights on [−1, 1] by the Golub–Welsch eigenvalue method (Golub &
# Welsch, Math. Comp. 23, 1969): the nodes are the eigenvalues of the Jacobi
# matrix of the Legendre recurrence, the weights twice the squared first
# components of its normalised eigenvectors.

using LinearAlgebra: SymTridiagonal, eigen

function _gauss_legendre(n::Int)
    β = [k / sqrt(4.0 * k^2 - 1.0) for k = 1:(n-1)]
    x, V = eigen(SymTridiagonal(zeros(n), β))
    return x, 2 .* abs2.(V[1, :])
end
