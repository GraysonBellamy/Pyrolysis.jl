"""
Reaction kinetics.

Reference: Design Document Section 2.4, 4.1
"""

"""
Depletion limiter settings for reaction rate computation.

The depletion limiter prevents numerical stiffness when reactant concentrations
approach zero by smoothly reducing the reaction rate using a hyperbolic tangent function.

# Fields
- `threshold::Float64`: Concentration below which rate limiting begins [kg/m³]
- `minimum::Float64`: Concentration below which rate is effectively zero [kg/m³]
- `enabled::Bool`: If false, no limiting is applied (may cause numerical issues)

# Notes
Lower thresholds are more physically accurate but may cause numerical stiffness.
For validation studies, temporarily disable with `enabled=false`.
"""
Base.@kwdef struct DepletionLimiter
	threshold::Float64 = 1.0    # kg/m³ - transition starts here
	minimum::Float64 = 0.01     # kg/m³ - rate effectively zero below this
	enabled::Bool = true        # Can disable for testing/validation
end

# Global default (can be overridden)
const DEFAULT_DEPLETION_LIMITER = Ref(DepletionLimiter())

"""
	set_depletion_limiter!(; threshold=1.0, minimum=0.01, enabled=true)

Configure the global depletion limiter settings.

# Arguments
- `threshold`: Concentration below which rate limiting begins [kg/m³]
- `minimum`: Concentration below which rate is effectively zero [kg/m³]
- `enabled`: If false, no limiting is applied (may cause numerical issues)

# Example
```julia
# More aggressive limiting (safer numerics, less accurate)
set_depletion_limiter!(threshold=5.0)

# Less limiting (more accurate, may be stiff)
set_depletion_limiter!(threshold=0.1)

# Disable for validation (use with caution!)
set_depletion_limiter!(enabled=false)
```
"""
function set_depletion_limiter!(; threshold::Float64 = 1.0, minimum::Float64 = 0.01,
	enabled::Bool = true)
	DEFAULT_DEPLETION_LIMITER[] = DepletionLimiter(threshold, minimum, enabled)
	return nothing
end

"""
	reaction_rate(rxn::Reaction, ξ::NTuple{NC,Float64}, T::Float64) -> Float64

Compute the reaction rate [kg/(m³·s)] for a single reaction.

r = A * exp(-E/RT) * ξ₁^n₁ * ξ₂^n₂ * H(T - T_min) * H(T_max - T)

# Depletion Protection
When reactant concentration falls below ξ_threshold (default 1.0 kg/m³),
the rate is smoothly reduced using a soft limiter to prevent numerical
instability from near-zero concentrations trying to react at full speed.

Use `set_depletion_limiter!()` to configure the limiter parameters globally.
"""

# PERFORMANCE: Fast power function for common reaction orders (0, 1)
# Avoids expensive x^y when y is a common integer value
@inline _fast_pow(x::Float64, n::Float64) = n == 1.0 ? x : (n == 0.0 ? 1.0 : x^n)

function reaction_rate(rxn::Reaction, ξ::NTuple{NC, Float64}, T::Float64;
	limiter::DepletionLimiter = DEFAULT_DEPLETION_LIMITER[]) where NC
	# Check temperature limits
	if T < rxn.T_min || T > rxn.T_max
		return 0.0
	end

	# Arrhenius rate constant
	k = rxn.A * exp(-rxn.E / (R_GAS * T))

	# Reactant concentrations with soft limiter for near-depletion
	ξ1_raw = rxn.reactant1 > 0 ? max(ξ[rxn.reactant1], 0.0) : 1.0
	ξ2_raw = rxn.reactant2 > 0 ? max(ξ[rxn.reactant2], 0.0) : 1.0

	# Apply soft depletion limiter if enabled
	if limiter.enabled
		# Smooth transition to zero as ξ → 0
		# This prevents the numerical stiffness when concentration is nearly depleted
		# f(ξ) = ξ * tanh(ξ/ξ_threshold) gives smooth approach to zero
		ξ1 = ξ1_raw < limiter.threshold ? ξ1_raw * tanh(ξ1_raw / limiter.threshold) : ξ1_raw
		ξ2 = ξ2_raw < limiter.threshold ? ξ2_raw * tanh(ξ2_raw / limiter.threshold) : ξ2_raw
	else
		ξ1 = ξ1_raw
		ξ2 = ξ2_raw
	end

	# Reaction rate (PERFORMANCE: use fast_pow for common orders 0.0, 1.0)
	r = k * _fast_pow(ξ1, rxn.order1) * _fast_pow(ξ2, rxn.order2)

	return r
end

"""
	evaluate_reactions!(rates::Vector{Float64}, heats::Vector{Float64},
					   ξ::NTuple{NC,Float64}, T::Float64, material::Material)

Evaluate all reaction rates and heats of reaction.

# Arguments
- `rates`: Output vector of reaction rates [kg/(m³·s)]
- `heats`: Output vector of heats of reaction [J/kg]
- `ξ`: Mass concentrations [kg/m³]
- `T`: Temperature [K]
- `material`: Material definition

This function is designed to be allocation-free using @generated for compile-time unrolling.

# State-Dependent Heats
Heats of reaction can depend on temperature only (standard) or on both
temperature and concentration state (e.g., moisture-dependent heat of sorption).
The concentration tuple ξ is passed to allow StateDependentProperty evaluation.
"""
function evaluate_reactions!(rates::Vector{Float64}, heats::Vector{Float64},
	ξ::NTuple{NC, Float64}, T::Float64, material::Material{NC, NR}) where {NC, NR}
	_evaluate_reactions_gen!(rates, heats, ξ, T, material.reactions, Val(NR))
	return nothing
end

# @generated version for allocation-free iteration over reactions
# Note: heats are evaluated with (T, ξ) to support StateDependentProperty
@generated function _evaluate_reactions_gen!(rates::Vector{Float64}, heats::Vector{Float64},
	ξ::NTuple{NC, Float64}, T::Float64, reactions::Tuple{Vararg{Any, NR}}, ::Val{NR}) where {NC, NR}
	exprs = Expr[]
	for i in 1:NR
		push!(exprs, quote
			@inbounds begin
				rxn = reactions[$i]
				rates[$i] = reaction_rate(rxn, ξ, T)
				heats[$i] = rxn.heat(T, ξ)  # Pass ξ for state-dependent heats
			end
		end)
	end
	push!(exprs, :(return nothing))
	return Expr(:block, exprs...)
end

"""
	species_source_terms!(sources::Vector{Float64}, rates::Vector{Float64},
						 material::Material)

Compute species source terms from reaction rates.

dξⱼ/dt = Σᵢ νᵢⱼ rᵢ

where νᵢⱼ is positive for products, negative for reactants.

# Arguments
- `sources`: Output vector of species source terms [kg/(m³·s)]
- `rates`: Reaction rates [kg/(m³·s)]
- `material`: Material definition
"""
function species_source_terms!(sources::NTuple{NC, Float64}, rates::Vector{Float64},
	material) where NC
	# Initialize sources to zero
	src = zeros(MVector{NC, Float64})

	NR = n_reactions(material)
	@inbounds for i in 1:NR
		rxn = material.reactions[i]
		r = rates[i]

		# Reactant consumption (negative)
		if rxn.reactant1 > 0
			src[rxn.reactant1] -= rxn.stoich[1] * r
		end
		if rxn.reactant2 > 0
			src[rxn.reactant2] -= rxn.stoich[2] * r
		end

		# Product generation (positive)
		if rxn.product1 > 0
			src[rxn.product1] += rxn.stoich[3] * r
		end
		if rxn.product2 > 0
			src[rxn.product2] += rxn.stoich[4] * r
		end
	end

	return NTuple{NC, Float64}(src)
end

"""
	species_source_terms!(sources::Vector{Float64}, rates::Vector{Float64},
						 material::Material, NC::Int)

In-place version that modifies a pre-allocated vector.
Uses @generated for allocation-free iteration over reactions.
"""
# Primary version that uses @generated for allocation-free iteration
function species_source_terms!(sources::Vector{Float64}, rates::Vector{Float64},
	material::Material{NC, NR}, ::Val{NC}) where {NC, NR}
	fill!(sources, 0.0)
	_species_source_terms_gen!(sources, rates, material.reactions, Val(NR))
	return nothing
end

# Overload accepting Int instead of Val{NC} - calls the same @generated implementation
function species_source_terms!(sources::Vector{Float64}, rates::Vector{Float64},
	material::Material{NC, NR}, NC_int::Int) where {NC, NR}
	fill!(sources, 0.0)
	_species_source_terms_gen!(sources, rates, material.reactions, Val(NR))
	return nothing
end

# @generated version for allocation-free iteration over reactions
@generated function _species_source_terms_gen!(sources::Vector{Float64}, rates::Vector{Float64},
	reactions::Tuple{Vararg{Any, NR}}, ::Val{NR}) where NR
	exprs = Expr[]
	for i in 1:NR
		push!(exprs, quote
			@inbounds begin
				rxn = reactions[$i]
				r = rates[$i]

				# Reactant consumption (negative)
				if rxn.reactant1 > 0
					sources[rxn.reactant1] -= rxn.stoich[1] * r
				end
				if rxn.reactant2 > 0
					sources[rxn.reactant2] -= rxn.stoich[2] * r
				end

				# Product generation (positive)
				if rxn.product1 > 0
					sources[rxn.product1] += rxn.stoich[3] * r
				end
				if rxn.product2 > 0
					sources[rxn.product2] += rxn.stoich[4] * r
				end
			end
		end)
	end
	push!(exprs, :(return nothing))
	return Expr(:block, exprs...)
end

# Fallback loop version for runtime NR
function _species_source_terms_loop!(sources::Vector{Float64}, rates::Vector{Float64},
	reactions, NR::Int)
	@inbounds for i in 1:NR
		rxn = reactions[i]
		r = rates[i]

		if rxn.reactant1 > 0
			sources[rxn.reactant1] -= rxn.stoich[1] * r
		end
		if rxn.reactant2 > 0
			sources[rxn.reactant2] -= rxn.stoich[2] * r
		end
		if rxn.product1 > 0
			sources[rxn.product1] += rxn.stoich[3] * r
		end
		if rxn.product2 > 0
			sources[rxn.product2] += rxn.stoich[4] * r
		end
	end
	return nothing
end

"""
	reaction_heat_source(rates::Vector{Float64}, heats::Vector{Float64}) -> Float64

Compute total heat source from reactions [W/m³].

Q = -Σᵢ hᵢ rᵢ

The negative sign converts from enthalpy convention (h < 0 exothermic) to
heat source convention (Q > 0 heats the system).

For exothermic reactions: h < 0, r > 0, so Q = -h*r > 0 (heating)
For endothermic reactions: h > 0, r > 0, so Q = -h*r < 0 (cooling)
"""
@inline function reaction_heat_source(rates::Vector{Float64}, heats::Vector{Float64})
	Q = 0.0
	@inbounds for i in eachindex(rates, heats)
		Q -= heats[i] * rates[i]  # Note the negative sign!
	end
	return Q
end

# ============================================================================
# Mass Conservation Validation
# ============================================================================

"""
	validate_species_source_conservation(sources::Vector{Float64}; rtol::Float64=1e-12) -> Bool

Check that species source terms sum to zero (mass conservation).

For properly balanced reactions, the total mass source should be exactly zero:
    Σⱼ dξⱼ/dt = 0

This function verifies this constraint within numerical tolerance.

# Arguments
- `sources`: Species source terms [kg/(m³·s)]
- `rtol`: Relative tolerance for sum compared to max source magnitude

# Returns
`true` if mass is conserved within tolerance, `false` otherwise.

# Example
```julia
sources = zeros(NC)
species_source_terms!(sources, rates, material, Val(NC))
if !validate_species_source_conservation(sources)
    @warn "Mass conservation violated in reaction sources"
end
```
"""
function validate_species_source_conservation(sources::Vector{Float64};
	rtol::Float64 = 1e-12)
	total = sum(sources)
	max_magnitude = maximum(abs, sources; init = 1.0)

	return abs(total) <= rtol * max_magnitude
end

"""
	species_source_terms_checked!(sources::Vector{Float64}, rates::Vector{Float64},
	                              material::Material{NC,NR}, ::Val{NC};
	                              throw_on_error::Bool=false) where {NC,NR}

Compute species sources with conservation check.

This is a wrapper around `species_source_terms!` that adds runtime validation
of mass conservation. Use during development and validation to catch errors early.

# Arguments
- `sources`: Output vector of species source terms [kg/(m³·s)]
- `rates`: Reaction rates [kg/(m³·s)]
- `material`: Material definition
- `Val(NC)`: Number of components (compile-time constant)
- `throw_on_error`: If true, raises an error on mass imbalance. Otherwise logs a warning.

# Performance Note
The conservation check adds minimal overhead (one sum + one max operation).
For production runs where reactions are known to be correct, use the unchecked
`species_source_terms!` directly.

# Example
```julia
# Development/validation mode
species_source_terms_checked!(sources, rates, material, Val(NC), throw_on_error=true)

# Production mode (skip check)
species_source_terms!(sources, rates, material, Val(NC))
```
"""
function species_source_terms_checked!(sources::Vector{Float64},
	rates::Vector{Float64},
	material::Material{NC, NR},
	::Val{NC};
	throw_on_error::Bool = false) where {NC, NR}
	# Compute sources normally
	species_source_terms!(sources, rates, material, Val(NC))

	# Validate conservation
	if !validate_species_source_conservation(sources)
		total = sum(sources)
		msg = "Reaction mass conservation violated: Σ(sources) = $total\n  Sources: $sources"

		if throw_on_error
			error(msg)
		else
			@warn msg
		end
	end

	return nothing
end
