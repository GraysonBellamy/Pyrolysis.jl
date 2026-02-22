"""
Mixture property computations.

Reference: Design Document Section 2.4

PERFORMANCE NOTE: All functions that iterate over material.components use
recursive dispatch to avoid allocations from indexing heterogeneous tuples
with runtime indices. Each component access uses a compile-time constant index.
"""

"""
	mass_fractions(ξ::NTuple{NC,Float64}) -> NTuple{NC,Float64}

Compute mass fractions from mass concentrations.

Yⱼ = ξⱼ / Σξⱼ
"""
@inline function mass_fractions(ξ::NTuple{NC, Float64}) where NC
	ρ_total = sum(ξ)
	if ρ_total < 1e-20
		return _zero_ntuple(Val(NC))
	end
	inv_ρ = 1.0 / ρ_total
	return _scale_ntuple(ξ, inv_ρ, Val(NC))
end

# Helper functions for allocation-free ntuple operations
# Explicit specializations for common NC values to avoid ntuple() allocations
@inline _zero_ntuple(::Val{1}) = (0.0,)
@inline _zero_ntuple(::Val{2}) = (0.0, 0.0)
@inline _zero_ntuple(::Val{3}) = (0.0, 0.0, 0.0)
@inline _zero_ntuple(::Val{4}) = (0.0, 0.0, 0.0, 0.0)
@inline _zero_ntuple(::Val{5}) = (0.0, 0.0, 0.0, 0.0, 0.0)
@inline _zero_ntuple(::Val{6}) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
@inline _zero_ntuple(::Val{7}) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
@inline _zero_ntuple(::Val{8}) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
@inline _zero_ntuple(::Val{9}) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
@inline _zero_ntuple(::Val{10}) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
@inline _zero_ntuple(::Val{11}) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
@inline _zero_ntuple(::Val{12}) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
@inline _zero_ntuple(::Val{13}) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
@inline _zero_ntuple(::Val{14}) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
@inline _zero_ntuple(::Val{15}) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
@inline _zero_ntuple(::Val{16}) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
@inline _zero_ntuple(::Val{N}) where N = ntuple(_ -> 0.0, Val(N))

@inline _scale_ntuple(x::NTuple{1, Float64}, s::Float64, ::Val{1}) = (x[1]*s,)
@inline _scale_ntuple(x::NTuple{2, Float64}, s::Float64, ::Val{2}) = (x[1]*s, x[2]*s)
@inline _scale_ntuple(x::NTuple{3, Float64}, s::Float64, ::Val{3}) = (x[1]*s, x[2]*s, x[3]*s)
@inline _scale_ntuple(x::NTuple{4, Float64}, s::Float64, ::Val{4}) = (x[1]*s, x[2]*s, x[3]*s, x[4]*s)
@inline _scale_ntuple(x::NTuple{5, Float64}, s::Float64, ::Val{5}) = (x[1]*s, x[2]*s, x[3]*s, x[4]*s, x[5]*s)
@inline _scale_ntuple(x::NTuple{6, Float64}, s::Float64, ::Val{6}) = (x[1]*s, x[2]*s, x[3]*s, x[4]*s, x[5]*s, x[6]*s)
@inline _scale_ntuple(x::NTuple{7, Float64}, s::Float64, ::Val{7}) = (x[1]*s, x[2]*s, x[3]*s, x[4]*s, x[5]*s, x[6]*s, x[7]*s)
@inline _scale_ntuple(x::NTuple{8, Float64}, s::Float64, ::Val{8}) = (x[1]*s, x[2]*s, x[3]*s, x[4]*s, x[5]*s, x[6]*s, x[7]*s, x[8]*s)
@inline _scale_ntuple(x::NTuple{9, Float64}, s::Float64, ::Val{9}) = (x[1]*s, x[2]*s, x[3]*s, x[4]*s, x[5]*s, x[6]*s, x[7]*s, x[8]*s, x[9]*s)
@inline _scale_ntuple(x::NTuple{10, Float64}, s::Float64, ::Val{10}) = (x[1]*s, x[2]*s, x[3]*s, x[4]*s, x[5]*s, x[6]*s, x[7]*s, x[8]*s, x[9]*s, x[10]*s)
@inline _scale_ntuple(x::NTuple{11, Float64}, s::Float64, ::Val{11}) = (x[1]*s, x[2]*s, x[3]*s, x[4]*s, x[5]*s, x[6]*s, x[7]*s, x[8]*s, x[9]*s, x[10]*s, x[11]*s)
@inline _scale_ntuple(x::NTuple{12, Float64}, s::Float64, ::Val{12}) = (x[1]*s, x[2]*s, x[3]*s, x[4]*s, x[5]*s, x[6]*s, x[7]*s, x[8]*s, x[9]*s, x[10]*s, x[11]*s, x[12]*s)
@inline _scale_ntuple(x::NTuple{13, Float64}, s::Float64, ::Val{13}) = (x[1]*s, x[2]*s, x[3]*s, x[4]*s, x[5]*s, x[6]*s, x[7]*s, x[8]*s, x[9]*s, x[10]*s, x[11]*s, x[12]*s, x[13]*s)
@inline _scale_ntuple(x::NTuple{14, Float64}, s::Float64, ::Val{14}) = (x[1]*s, x[2]*s, x[3]*s, x[4]*s, x[5]*s, x[6]*s, x[7]*s, x[8]*s, x[9]*s, x[10]*s, x[11]*s, x[12]*s, x[13]*s, x[14]*s)
@inline _scale_ntuple(x::NTuple{15, Float64}, s::Float64, ::Val{15}) = (x[1]*s, x[2]*s, x[3]*s, x[4]*s, x[5]*s, x[6]*s, x[7]*s, x[8]*s, x[9]*s, x[10]*s, x[11]*s, x[12]*s, x[13]*s, x[14]*s, x[15]*s)
@inline _scale_ntuple(x::NTuple{16, Float64}, s::Float64, ::Val{16}) = (x[1]*s, x[2]*s, x[3]*s, x[4]*s, x[5]*s, x[6]*s, x[7]*s, x[8]*s, x[9]*s, x[10]*s, x[11]*s, x[12]*s, x[13]*s, x[14]*s, x[15]*s, x[16]*s)
@inline _scale_ntuple(x::NTuple{N, Float64}, s::Float64, ::Val{N}) where N = ntuple(j -> x[j]*s, Val(N))

"""
	total_density(ξ::NTuple{NC,Float64}) -> Float64

Compute total density (sum of mass concentrations).
"""
@inline total_density(ξ::NTuple{NC, Float64}) where NC = sum(ξ)

# ============================================================================
# Recursive helpers for heterogeneous tuple iteration
# These avoid runtime indexing by using compile-time dispatch
# ============================================================================

# Compute specific volume from components recursively
@inline function _sum_specific_vol(ξ, T, comps::Tuple{C}, ::Val{1}) where C
	ρ = comps[1].density(T)
	return ξ[1] > 0.0 && ρ > 0.0 ? ξ[1] / ρ : 0.0
end

@inline function _sum_specific_vol(ξ, T, comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	v = 0.0
	ρ1 = comps[1].density(T)
	if ξ[1] > 0.0 && ρ1 > 0.0
		v += ξ[1] / ρ1
	end
	ρ2 = comps[2].density(T)
	if ξ[2] > 0.0 && ρ2 > 0.0
		v += ξ[2] / ρ2
	end
	return v
end

@inline function _sum_specific_vol(ξ, T, comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	v = 0.0
	ρ1 = comps[1].density(T)
	if ξ[1] > 0.0 && ρ1 > 0.0
		;
		v += ξ[1] / ρ1;
	end
	ρ2 = comps[2].density(T)
	if ξ[2] > 0.0 && ρ2 > 0.0
		;
		v += ξ[2] / ρ2;
	end
	ρ3 = comps[3].density(T)
	if ξ[3] > 0.0 && ρ3 > 0.0
		;
		v += ξ[3] / ρ3;
	end
	return v
end

# Generic fallback for any number of components (NC > 3) using generated function
@inline function _sum_specific_vol(ξ, T, comps::Tuple, ::Val{NC}) where NC
	return _sum_specific_vol_gen(ξ, T, comps)
end

# Compute volume fractions tuple
@inline function _vol_fracs(ξ, T, inv_v_total, comps::Tuple{C}, ::Val{1}) where C
	ρ = comps[1].density(T)
	v = ξ[1] > 0.0 && ρ > 0.0 ? ξ[1] / ρ * inv_v_total : 0.0
	return (v,)
end

@inline function _vol_fracs(ξ, T, inv_v_total, comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	ρ1 = comps[1].density(T)
	v1 = ξ[1] > 0.0 && ρ1 > 0.0 ? ξ[1] / ρ1 * inv_v_total : 0.0
	ρ2 = comps[2].density(T)
	v2 = ξ[2] > 0.0 && ρ2 > 0.0 ? ξ[2] / ρ2 * inv_v_total : 0.0
	return (v1, v2)
end

@inline function _vol_fracs(ξ, T, inv_v_total, comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	ρ1 = comps[1].density(T)
	v1 = ξ[1] > 0.0 && ρ1 > 0.0 ? ξ[1] / ρ1 * inv_v_total : 0.0
	ρ2 = comps[2].density(T)
	v2 = ξ[2] > 0.0 && ρ2 > 0.0 ? ξ[2] / ρ2 * inv_v_total : 0.0
	ρ3 = comps[3].density(T)
	v3 = ξ[3] > 0.0 && ρ3 > 0.0 ? ξ[3] / ρ3 * inv_v_total : 0.0
	return (v1, v2, v3)
end

# Generic fallback for any number of components (NC > 3) using generated function
@inline function _vol_fracs(ξ, T, inv_v_total, comps::Tuple, ::Val{NC}) where NC
	return _vol_fracs_gen(ξ, T, inv_v_total, comps)
end

"""
	volume_fractions(ξ::NTuple{NC,Float64}, T::Float64, material::Material) -> NTuple{NC,Float64}

Compute volume fractions from mass concentrations and component densities.

vⱼ = (ξⱼ/ρⱼ) / Σᵢ(ξᵢ/ρᵢ)

Uses dispatch-based unrolling to avoid allocation from heterogeneous tuple indexing.
"""
@inline function volume_fractions(ξ::NTuple{NC, Float64}, T::Float64, material) where NC
	comps = material.components
	v_total = _sum_specific_vol(ξ, T, comps, Val(NC))

	if v_total < 1e-20
		return _zero_ntuple(Val(NC))
	end

	inv_v_total = 1.0 / v_total
	return _vol_fracs(ξ, T, inv_v_total, comps, Val(NC))
end

"""
	volume_fraction_derivatives(ξ, T, material) -> ∂vⱼ/∂ξᵢ

Compute analytical derivatives of volume fractions for Jacobian assembly.

For volume fraction vⱼ = (ξⱼ/ρⱼ) / V_total where V_total = Σₖ(ξₖ/ρₖ):

∂vⱼ/∂ξᵢ = δᵢⱼ/(ρⱼ × V_total) - vⱼ/(ρᵢ × V_total)

Returns an NC×NC matrix where element [j,i] = ∂vⱼ/∂ξᵢ
"""
@inline function volume_fraction_derivatives(ξ::NTuple{NC, Float64}, T::Float64, material) where NC
	comps = material.components
	v_total = _sum_specific_vol(ξ, T, comps, Val(NC))

	if v_total < 1e-20
		# Return zeros if no volume
		return ntuple(_ -> ntuple(_ -> 0.0, Val(NC)), Val(NC))
	end

	# Get current volume fractions
	v = volume_fractions(ξ, T, material)

	# Get component densities
	ρ = ntuple(i -> comps[i].density(T), Val(NC))

	# Compute derivatives: ∂vⱼ/∂ξᵢ stored as [j, i]
	# Returns tuple of tuples: outer index j (component whose v is differentiated),
	# inner index i (concentration being perturbed)
	return ntuple(Val(NC)) do j
		ntuple(Val(NC)) do i
			# Kronecker delta term
			δ_term = (i == j) ? 1.0 / (ρ[j] * v_total) : 0.0

			# Common term
			common_term = v[j] / (ρ[i] * v_total)

			δ_term - common_term
		end
	end
end

# ============================================================================
# Solid volume fractions (excluding components with swelling_factor ≈ 0)
# Used for conductivity and gas transfer calculations per ThermaKin convention:
# "when γ = 0, gaseous components do not contribute to the thermal conductivity"
# ============================================================================

# Sum specific volumes only for components with swelling_factor > 0
@inline function _sum_solid_specific_vol(ξ, T, comps::Tuple{C}, ::Val{1}) where C
	c = comps[1]
	if c.phase == GAS && c.swelling_factor < 1e-10
		return 0.0
	end
	ρ = c.density(T)
	return ξ[1] > 0.0 && ρ > 0.0 ? ξ[1] / ρ : 0.0
end

@inline function _sum_solid_specific_vol(ξ, T, comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	v = 0.0
	c1 = comps[1]
	if !(c1.phase == GAS && c1.swelling_factor < 1e-10)
		ρ1 = c1.density(T)
		if ξ[1] > 0.0 && ρ1 > 0.0
			v += ξ[1] / ρ1
		end
	end
	c2 = comps[2]
	if !(c2.phase == GAS && c2.swelling_factor < 1e-10)
		ρ2 = c2.density(T)
		if ξ[2] > 0.0 && ρ2 > 0.0
			v += ξ[2] / ρ2
		end
	end
	return v
end

@inline function _sum_solid_specific_vol(ξ, T, comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	v = 0.0
	c1 = comps[1]
	if !(c1.phase == GAS && c1.swelling_factor < 1e-10)
		ρ1 = c1.density(T)
		if ξ[1] > 0.0 && ρ1 > 0.0
			v += ξ[1] / ρ1
		end
	end
	c2 = comps[2]
	if !(c2.phase == GAS && c2.swelling_factor < 1e-10)
		ρ2 = c2.density(T)
		if ξ[2] > 0.0 && ρ2 > 0.0
			v += ξ[2] / ρ2
		end
	end
	c3 = comps[3]
	if !(c3.phase == GAS && c3.swelling_factor < 1e-10)
		ρ3 = c3.density(T)
		if ξ[3] > 0.0 && ρ3 > 0.0
			v += ξ[3] / ρ3
		end
	end
	return v
end

# Generic fallback for any number of components (NC > 3) using generated function
@inline function _sum_solid_specific_vol(ξ, T, comps::Tuple, ::Val{NC}) where NC
	return _sum_solid_specific_vol_gen(ξ, T, comps)
end

# Compute solid volume fractions tuple (0 for non-swelling gases)
@inline function _solid_vol_fracs(ξ, T, inv_v_total, comps::Tuple{C}, ::Val{1}) where C
	c = comps[1]
	if c.phase == GAS && c.swelling_factor < 1e-10
		return (0.0,)
	end
	ρ = c.density(T)
	v = ξ[1] > 0.0 && ρ > 0.0 ? ξ[1] / ρ * inv_v_total : 0.0
	return (v,)
end

@inline function _solid_vol_fracs(ξ, T, inv_v_total, comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	c1 = comps[1]
	if !(c1.phase == GAS && c1.swelling_factor < 1e-10)
		ρ1 = c1.density(T)
		v1 = ξ[1] > 0.0 && ρ1 > 0.0 ? ξ[1] / ρ1 * inv_v_total : 0.0
	else
		v1 = 0.0
	end
	c2 = comps[2]
	if !(c2.phase == GAS && c2.swelling_factor < 1e-10)
		ρ2 = c2.density(T)
		v2 = ξ[2] > 0.0 && ρ2 > 0.0 ? ξ[2] / ρ2 * inv_v_total : 0.0
	else
		v2 = 0.0
	end
	return (v1, v2)
end

@inline function _solid_vol_fracs(ξ, T, inv_v_total, comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	c1 = comps[1]
	if !(c1.phase == GAS && c1.swelling_factor < 1e-10)
		ρ1 = c1.density(T)
		v1 = ξ[1] > 0.0 && ρ1 > 0.0 ? ξ[1] / ρ1 * inv_v_total : 0.0
	else
		v1 = 0.0
	end
	c2 = comps[2]
	if !(c2.phase == GAS && c2.swelling_factor < 1e-10)
		ρ2 = c2.density(T)
		v2 = ξ[2] > 0.0 && ρ2 > 0.0 ? ξ[2] / ρ2 * inv_v_total : 0.0
	else
		v2 = 0.0
	end
	c3 = comps[3]
	if !(c3.phase == GAS && c3.swelling_factor < 1e-10)
		ρ3 = c3.density(T)
		v3 = ξ[3] > 0.0 && ρ3 > 0.0 ? ξ[3] / ρ3 * inv_v_total : 0.0
	else
		v3 = 0.0
	end
	return (v1, v2, v3)
end

# Generic fallback for any number of components (NC > 3) using generated function
@inline function _solid_vol_fracs(ξ, T, inv_v_total, comps::Tuple, ::Val{NC}) where NC
	return _solid_vol_fracs_gen(ξ, T, inv_v_total, comps)
end

"""
	solid_volume_fractions(ξ::NTuple{NC,Float64}, T::Float64, material::Material) -> NTuple{NC,Float64}

Compute volume fractions for solid/liquid components only (excludes components with swelling_factor ≈ 0).

This matches ThermaKin behavior: "when γ = 0, gaseous components do not contribute to the
thermal conductivity" (ThermaKin2Ds Reference, Equations 8-10).

For gases with swelling_factor = 0 (default), their volume fraction is returned as 0.
Volume fractions are normalized to sum to 1 over the contributing (non-zero swelling) components.

Used for:
- Effective thermal conductivity
- Effective gas transfer coefficient
"""
@inline function solid_volume_fractions(ξ::NTuple{NC, Float64}, T::Float64, material) where NC
	comps = material.components
	v_total = _sum_solid_specific_vol(ξ, T, comps, Val(NC))

	if v_total < 1e-20
		return _zero_ntuple(Val(NC))
	end

	inv_v_total = 1.0 / v_total
	return _solid_vol_fracs(ξ, T, inv_v_total, comps, Val(NC))
end

# ============================================================================
# mixture_density with swelling
# ============================================================================

@inline function _specific_vol_with_swelling(Y, T, comps::Tuple{C}, ::Val{1}) where C
	c = comps[1]
	ρ = c.density(T)
	if c.phase == GAS
		return c.swelling_factor * Y[1] / ρ
	else
		return Y[1] / ρ
	end
end

@inline function _specific_vol_with_swelling(Y, T, comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	v = 0.0
	c1 = comps[1]
	ρ1 = c1.density(T)
	if c1.phase == GAS
		v += c1.swelling_factor * Y[1] / ρ1
	else
		v += Y[1] / ρ1
	end
	c2 = comps[2]
	ρ2 = c2.density(T)
	if c2.phase == GAS
		v += c2.swelling_factor * Y[2] / ρ2
	else
		v += Y[2] / ρ2
	end
	return v
end

@inline function _specific_vol_with_swelling(Y, T, comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	v = 0.0
	c1 = comps[1]
	ρ1 = c1.density(T)
	v += c1.phase == GAS ? c1.swelling_factor * Y[1] / ρ1 : Y[1] / ρ1
	c2 = comps[2]
	ρ2 = c2.density(T)
	v += c2.phase == GAS ? c2.swelling_factor * Y[2] / ρ2 : Y[2] / ρ2
	c3 = comps[3]
	ρ3 = c3.density(T)
	v += c3.phase == GAS ? c3.swelling_factor * Y[3] / ρ3 : Y[3] / ρ3
	return v
end

# Generic fallback for any number of components (NC > 3) using generated function
@generated function _specific_vol_with_swelling(Y, T, comps::Tuple, ::Val{NC}) where NC
	expr = quote
		v = 0.0
	end
	for i in 1:NC
		push!(expr.args, quote
			c = comps[$i]
			ρ = c.density(T)
			if c.phase == GAS
				v += c.swelling_factor * Y[$i] / ρ
			else
				v += Y[$i] / ρ
			end
		end)
	end
	push!(expr.args, :(return v))
	return expr
end

# Optimized version that computes mass fractions inline to avoid intermediate NTuple allocation
@inline function _specific_vol_with_swelling_direct(ξ, inv_ρ, T, comps::Tuple{C}, ::Val{1}) where C
	c = comps[1]
	ρ_comp = c.density(T)
	Y_i = ξ[1] * inv_ρ
	if c.phase == GAS
		return c.swelling_factor * Y_i / ρ_comp
	else
		return Y_i / ρ_comp
	end
end

@inline function _specific_vol_with_swelling_direct(ξ, inv_ρ, T, comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	v = 0.0
	c1 = comps[1]
	ρ1 = c1.density(T)
	Y1 = ξ[1] * inv_ρ
	v += c1.phase == GAS ? c1.swelling_factor * Y1 / ρ1 : Y1 / ρ1
	c2 = comps[2]
	ρ2 = c2.density(T)
	Y2 = ξ[2] * inv_ρ
	v += c2.phase == GAS ? c2.swelling_factor * Y2 / ρ2 : Y2 / ρ2
	return v
end

@inline function _specific_vol_with_swelling_direct(ξ, inv_ρ, T, comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	v = 0.0
	c1 = comps[1]
	ρ1 = c1.density(T)
	Y1 = ξ[1] * inv_ρ
	v += c1.phase == GAS ? c1.swelling_factor * Y1 / ρ1 : Y1 / ρ1
	c2 = comps[2]
	ρ2 = c2.density(T)
	Y2 = ξ[2] * inv_ρ
	v += c2.phase == GAS ? c2.swelling_factor * Y2 / ρ2 : Y2 / ρ2
	c3 = comps[3]
	ρ3 = c3.density(T)
	Y3 = ξ[3] * inv_ρ
	v += c3.phase == GAS ? c3.swelling_factor * Y3 / ρ3 : Y3 / ρ3
	return v
end

# Generic fallback using generated function
@generated function _specific_vol_with_swelling_direct(ξ, inv_ρ, T, comps::Tuple, ::Val{NC}) where NC
	expr = quote
		v = 0.0
	end
	for i in 1:NC
		push!(expr.args, quote
			c = comps[$i]
			ρ = c.density(T)
			Y_i = ξ[$i] * inv_ρ  # Compute mass fraction inline
			if c.phase == GAS
				v += c.swelling_factor * Y_i / ρ
			else
				v += Y_i / ρ
			end
		end)
	end
	push!(expr.args, :(return v))
	return expr
end

"""
	mixture_density(ξ::NTuple{NC,Float64}, T::Float64, material::Material) -> Float64

Compute mixture density from mass concentrations.

For a mixture of solids and gases with swelling:
1/ρ = Σ_solid(Yⱼ/ρⱼ) + Σ_liquid(Yⱼ/ρⱼ) + γ·Σ_gas(Yⱼ/ρⱼ)

Where γ is the average swelling factor.

PERFORMANCE: Uses optimized version that computes mass fractions inline to avoid
allocating intermediate NTuple.
"""
@inline function mixture_density(ξ::NTuple{NC, Float64}, T::Float64, material) where NC
	ρ_total = total_density(ξ)
	if ρ_total < 1e-20
		return 0.0
	end

	inv_ρ = 1.0 / ρ_total
	# Use optimized version that computes mass fractions inline
	v_specific = _specific_vol_with_swelling_direct(ξ, inv_ρ, T, material.components, Val(NC))

	if v_specific < 1e-20
		return 0.0
	end

	return 1.0 / v_specific
end

# ============================================================================
# Porosity
# ============================================================================

@inline function _porosity_sum(v, comps::Tuple{C}, ::Val{1}) where C
	return comps[1].phase == GAS ? v[1] : 0.0
end

@inline function _porosity_sum(v, comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	φ = 0.0
	if comps[1].phase == GAS
		;
		φ += v[1];
	end
	if comps[2].phase == GAS
		;
		φ += v[2];
	end
	return φ
end

@inline function _porosity_sum(v, comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	φ = 0.0
	if comps[1].phase == GAS
		;
		φ += v[1];
	end
	if comps[2].phase == GAS
		;
		φ += v[2];
	end
	if comps[3].phase == GAS
		;
		φ += v[3];
	end
	return φ
end

# Explicit unrolled version for 4 components (common case)
@inline function _porosity_sum(v, comps::Tuple{C1, C2, C3, C4}, ::Val{4}) where {C1, C2, C3, C4}
	φ = 0.0
	if comps[1].phase == GAS
		φ += v[1]
	end
	if comps[2].phase == GAS
		φ += v[2]
	end
	if comps[3].phase == GAS
		φ += v[3]
	end
	if comps[4].phase == GAS
		φ += v[4]
	end
	return φ
end

# Generic fallback using generated function for NC > 4
@inline function _porosity_sum(v, comps::Tuple, ::Val{NC}) where NC
	return _porosity_sum_gen(v, comps)
end

@generated function _porosity_sum_gen(v::NTuple{NC, Float64}, comps::Tuple) where NC
	expr = quote
		φ = 0.0
	end
	for i in 1:NC
		push!(expr.args, quote
			if comps[$i].phase == GAS
				φ += v[$i]
			end
		end)
	end
	push!(expr.args, :(return φ))
	return expr
end

"""
	porosity(ξ::NTuple{NC,Float64}, T::Float64, material::Material) -> Float64

Compute porosity (gas volume fraction).
"""
@inline function porosity(ξ::NTuple{NC, Float64}, T::Float64, material) where NC
	v = volume_fractions(ξ, T, material)
	return _porosity_sum(v, material.components, Val(NC))
end

"""
	porosity_derivatives(ξ::NTuple{NC,Float64}, T::Float64, material::Material) -> NTuple{NC,Float64}

Compute derivatives of porosity with respect to each component concentration.

For porosity φ = 1 - Σ(ξⱼ/ρⱼ) for solid/liquid components:
- ∂φ/∂ξⱼ = -1/ρⱼ for solid/liquid components
- ∂φ/∂ξⱼ = 0 for gas components (gas fills the pore space but doesn't affect solid volume)

# Arguments
- `ξ`: Mass concentrations [kg/m³]
- `T`: Temperature [K]
- `material`: Material definition

# Returns
Tuple of ∂φ/∂ξⱼ for each component.
"""
@inline function porosity_derivatives(ξ::NTuple{NC, Float64}, T::Float64,
	material::Material{NC}) where NC
	return _porosity_deriv_gen(T, material.components, Val(NC))
end

# Generated function for allocation-free porosity derivatives
@generated function _porosity_deriv_gen(T::Float64, comps::Tuple, ::Val{NC}) where NC
	var_exprs = []
	for i in 1:NC
		sym = Symbol("dφ_dξ$i")
		push!(var_exprs, quote
			comp = comps[$i]
			if comp.phase == GAS
				$sym = 0.0
			else
				# ∂φ/∂ξⱼ = -∂(Σ ξⱼ/ρⱼ)/∂ξⱼ = -1/ρⱼ
				ρ = comp.density(T)
				$sym = ρ > 0.0 ? -1.0 / ρ : 0.0
			end
		end)
	end

	tuple_symbols = [Symbol("dφ_dξ$i") for i in 1:NC]
	tuple_expr = Expr(:tuple, tuple_symbols...)

	return quote
		$(var_exprs...)
		return $tuple_expr
	end
end

"""
	porosity_derivatives!(∂φ_∂ξ::AbstractVector, T::Float64, material::Material{NC}) where NC

Compute derivatives of porosity in-place (allocation-free).

Same as `porosity_derivatives` but writes to a pre-allocated vector instead of
returning a tuple. This eliminates heap allocations when called in a loop.

# Arguments
- `∂φ_∂ξ`: Output vector of length NC for concentration derivatives [m³/kg]
- `T`: Temperature [K]  
- `material`: Material definition

# Notes
The output vector `∂φ_∂ξ` is modified in-place.
"""
@inline function porosity_derivatives!(∂φ_∂ξ::AbstractVector, T::Float64,
	material::Material{NC}) where NC
	return _porosity_deriv_inplace!(∂φ_∂ξ, T, material.components, Val(NC))
end

# Generated function for in-place porosity derivatives
@generated function _porosity_deriv_inplace!(∂φ_∂ξ::AbstractVector, T::Float64,
	comps::Tuple, ::Val{NC}) where NC
	expr = Expr(:block)

	for i in 1:NC
		push!(expr.args, quote
			comp = comps[$i]
			if comp.phase == GAS
				@inbounds ∂φ_∂ξ[$i] = 0.0
			else
				# ∂φ/∂ξⱼ = -∂(Σ ξⱼ/ρⱼ)/∂ξⱼ = -1/ρⱼ
				ρ = comp.density(T)
				@inbounds ∂φ_∂ξ[$i] = ρ > 0.0 ? -1.0 / ρ : 0.0
			end
		end)
	end

	push!(expr.args, :(return nothing))
	return expr
end

# ============================================================================
# Generated functions for allocation-free computation with many components
# These generate fully unrolled code at compile time, avoiding runtime
# indexing into heterogeneous tuples which would allocate.
# ============================================================================

"""
Generated function to compute sum of specific volumes.
Generates fully unrolled code like:
	v = 0.0
	ρ1 = comps[1].density(T); if ξ[1] > 0.0 && ρ1 > 0.0; v += ξ[1] / ρ1; end
	ρ2 = comps[2].density(T); if ξ[2] > 0.0 && ρ2 > 0.0; v += ξ[2] / ρ2; end
	...
	return v
"""
@generated function _sum_specific_vol_gen(ξ::NTuple{NC, Float64}, T::Float64,
	comps::Tuple) where NC
	expr = quote
		v = 0.0
	end
	for i in 1:NC
		push!(expr.args, quote
			ρ = comps[$i].density(T)
			if ξ[$i] > 0.0 && ρ > 0.0
				v += ξ[$i] / ρ
			end
		end)
	end
	push!(expr.args, :(return v))
	return expr
end

"""
Generated function to compute volume fractions tuple.
Generates fully unrolled code returning a compile-time sized tuple.
"""
@generated function _vol_fracs_gen(ξ::NTuple{NC, Float64}, T::Float64,
	inv_v_total::Float64, comps::Tuple) where NC
	# Generate variable declarations for each volume fraction
	var_exprs = []
	for i in 1:NC
		sym = Symbol("v$i")
		push!(var_exprs, quote
			ρ = comps[$i].density(T)
			$sym = ξ[$i] > 0.0 && ρ > 0.0 ? ξ[$i] / ρ * inv_v_total : 0.0
		end)
	end

	# Generate tuple return
	tuple_symbols = [Symbol("v$i") for i in 1:NC]
	tuple_expr = Expr(:tuple, tuple_symbols...)

	return quote
		$(var_exprs...)
		return $tuple_expr
	end
end

"""Generated function to compute solid volume fractions sum (excluding non-swelling gases).
"""
@generated function _sum_solid_specific_vol_gen(ξ::NTuple{NC, Float64}, T::Float64,
	comps::Tuple) where NC
	expr = quote
		v = 0.0
	end
	for i in 1:NC
		push!(expr.args, quote
			c = comps[$i]
			if !(c.phase == GAS && c.swelling_factor < 1e-10)
				ρ = c.density(T)
				if ξ[$i] > 0.0 && ρ > 0.0
					v += ξ[$i] / ρ
				end
			end
		end)
	end
	push!(expr.args, :(return v))
	return expr
end

"""
Generated function to compute solid volume fractions tuple.
"""
@generated function _solid_vol_fracs_gen(ξ::NTuple{NC, Float64}, T::Float64,
	inv_v_total::Float64, comps::Tuple) where NC
	var_exprs = []
	for i in 1:NC
		sym = Symbol("v$i")
		push!(var_exprs, quote
			c = comps[$i]
			if !(c.phase == GAS && c.swelling_factor < 1e-10)
				ρ = c.density(T)
				$sym = ξ[$i] > 0.0 && ρ > 0.0 ? ξ[$i] / ρ * inv_v_total : 0.0
			else
				$sym = 0.0
			end
		end)
	end

	tuple_symbols = [Symbol("v$i") for i in 1:NC]
	tuple_expr = Expr(:tuple, tuple_symbols...)

	return quote
		$(var_exprs...)
		return $tuple_expr
	end
end
