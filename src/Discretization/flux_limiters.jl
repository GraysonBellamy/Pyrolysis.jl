"""
Flux limiters for TVD schemes.

TVD (Total Variation Diminishing) limiters prevent spurious oscillations
in second-order schemes while maintaining accuracy in smooth regions.

Available limiters:
- minmod: Most diffusive, but very robust
- vanleer: Good balance of accuracy and robustness
- superbee: Least diffusive, can be slightly oscillatory
"""

# ============================================================================
# Limiter Enum for dispatch
# ============================================================================

"""
	Limiter

Enumeration of available flux limiters.
"""
@enum Limiter begin
    MINMOD = 1
    VANLEER = 2
    SUPERBEE = 3
    NONE = 4      # No limiting (first-order upwind)
    CENTRAL = 5   # No limiting (second-order central)
end

# ============================================================================
# Limiter Functions (scalar r input)
# ============================================================================

"""
	minmod(r::Real) -> Float64

Minmod limiter function.

ψ(r) = max(0, min(1, r))

Most diffusive limiter, but guarantees TVD.
"""
@inline function minmod(r::Real)
    return max(0.0, min(1.0, r))
end

"""
	vanleer(r::Real) -> Float64

Van Leer limiter function.

ψ(r) = (r + |r|) / (1 + |r|)

Good balance between diffusion and accuracy.
"""
@inline function vanleer(r::Real)
    ar = abs(r)
    return (r + ar) / (1.0 + ar)
end

"""
	superbee(r::Real) -> Float64

Superbee limiter function.

ψ(r) = max(0, min(2r, 1), min(r, 2))

Least diffusive standard limiter.
"""
@inline function superbee(r::Real)
    return max(0.0, min(2.0 * r, 1.0), min(r, 2.0))
end

# ============================================================================
# Dispatch-based limiter application
# ============================================================================

"""
	apply_limiter(r::Float64, limiter::Limiter) -> Float64

Apply the specified limiter to slope ratio r.
"""
@inline function apply_limiter(r::Real, limiter::Limiter)
    # Promote once so every branch returns the same type — literal 0.0/1.0
    # here would make the return type Union{Float64, Dual} under AD, and the
    # union boxes in the limited-gradient arithmetic downstream.
    rf = float(r)
    if limiter == MINMOD
        return minmod(rf)
    elseif limiter == VANLEER
        return vanleer(rf)
    elseif limiter == SUPERBEE
        return superbee(rf)
    elseif limiter == NONE
        return zero(rf)
    elseif limiter == CENTRAL
        return one(rf)
    else
        return minmod(rf)  # Default to minmod
    end
end
