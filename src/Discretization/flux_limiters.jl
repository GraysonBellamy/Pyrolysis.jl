"""
Flux limiters for TVD schemes.

Reference: Design Document Section 8.2

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

"""
	no_limiter(r::Real) -> Float64

No limiting (returns 0 for first-order upwind).
"""
@inline no_limiter(r::Real) = 0.0

"""
	central_limiter(r::Real) -> Float64

Central scheme (returns 1 for linear interpolation).
"""
@inline central_limiter(r::Real) = 1.0

# ============================================================================
# Dispatch-based limiter application
# ============================================================================

"""
	apply_limiter(r::Float64, limiter::Limiter) -> Float64

Apply the specified limiter to slope ratio r.
"""
@inline function apply_limiter(r::Float64, limiter::Limiter)
	if limiter == MINMOD
		return minmod(r)
	elseif limiter == VANLEER
		return vanleer(r)
	elseif limiter == SUPERBEE
		return superbee(r)
	elseif limiter == NONE
		return 0.0
	elseif limiter == CENTRAL
		return 1.0
	else
		return minmod(r)  # Default to minmod
	end
end

# ============================================================================
# Limited Face Value Computation
# ============================================================================

"""
	limited_face_value(φ_LL::Float64, φ_L::Float64, φ_R::Float64, 
					   limiter::Limiter) -> Float64

Compute second-order limited face value for upwind scheme.

For flow from L to R (positive velocity):
	φ_f = φ_L + 0.5 * ψ(r) * (φ_R - φ_L)

where r = (φ_L - φ_LL) / (φ_R - φ_L)

# Arguments
- `φ_LL`: Value at upwind-upwind cell (L's left neighbor)
- `φ_L`: Value at upwind cell (left of face)
- `φ_R`: Value at downwind cell (right of face)
- `limiter`: Limiter type to apply

# Returns
Limited face value
"""
@inline function limited_face_value(φ_LL::Float64, φ_L::Float64, φ_R::Float64,
	limiter::Limiter)
	# Compute slope ratio
	denom = φ_R - φ_L
	if abs(denom) < 1e-20
		return φ_L  # Uniform field
	end
	r = (φ_L - φ_LL) / denom

	# Apply limiter
	ψ = apply_limiter(r, limiter)

	return φ_L + 0.5 * ψ * denom
end

"""
	limited_face_value(φ_L::Float64, φ_R::Float64, ∇φ_L::Float64, 
					   d_Lf::Float64, limiter::Limiter) -> Float64

Compute limited face value using cell gradient.

φ_f = φ_L + ψ * ∇φ_L * d_Lf

where ψ is computed from comparing gradient-based and central reconstructions.
"""
@inline function limited_face_value(φ_L::Float64, φ_R::Float64, ∇φ_L::Float64,
	d_Lf::Float64, d_LR::Float64, limiter::Limiter)
	# Slope ratio based on gradient vs central difference
	central_slope = (φ_R - φ_L) / d_LR
	if abs(central_slope) < 1e-20
		return φ_L
	end
	r = ∇φ_L / central_slope

	ψ = apply_limiter(r, limiter)

	# Limited reconstruction
	return φ_L + ψ * ∇φ_L * d_Lf
end

# ============================================================================
# Symmetric Limiter (for central schemes)
# ============================================================================

"""
	symmetric_limiter(r::Float64, limiter::Limiter) -> Float64

Compute symmetric limiter value ψ(r) * r for use with central schemes.
"""
@inline function symmetric_limiter(r::Float64, limiter::Limiter)
	return apply_limiter(r, limiter) * r
end

# ============================================================================
# Two-sided limiting for convection
# ============================================================================

"""
	convective_face_value(φ_LL::Float64, φ_L::Float64, φ_R::Float64, φ_RR::Float64,
						  v::Float64, limiter::Limiter) -> Float64

Compute limited face value for convective flux with automatic upwinding.

# Arguments
- `φ_LL`: Value at cell two to the left
- `φ_L`: Value at cell to the left
- `φ_R`: Value at cell to the right
- `φ_RR`: Value at cell two to the right
- `v`: Velocity at face (positive = flow from L to R)
- `limiter`: Limiter type

# Returns
Limited face value based on upwind direction.
"""
@inline function convective_face_value(φ_LL::Float64, φ_L::Float64,
	φ_R::Float64, φ_RR::Float64,
	v::Float64, limiter::Limiter)
	if v >= 0.0
		# Flow from L to R: upwind is L, need L's left neighbor (LL)
		return limited_face_value(φ_LL, φ_L, φ_R, limiter)
	else
		# Flow from R to L: upwind is R, need R's right neighbor (RR)
		return limited_face_value(φ_RR, φ_R, φ_L, limiter)
	end
end
