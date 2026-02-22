"""
Material type definitions.

A Material combines components and reactions into a complete
description of a pyrolyzing substance.

Reference: Design Document Section 6.2
"""

# ============================================================================
# Mixing Rule Enumeration
# ============================================================================

"""
	MixingRule

Enumeration of property mixing rules for multi-component materials.

- `PARALLEL`: Parallel (volume-weighted) mixing: k = Σᵢ vᵢkᵢ
- `SERIES`: Series (harmonic) mixing: 1/k = Σᵢ vᵢ/kᵢ
- `WEIGHTED`: Weighted combination: k = β·k_parallel + (1-β)·k_series
- `BRUGGEMAN`: Bruggeman effective medium: Σᵢ vᵢ(kᵢ - k_eff)/(kᵢ + 2k_eff) = 0
"""
@enum MixingRule begin
	PARALLEL = 1
	SERIES = 2
	WEIGHTED = 3
	BRUGGEMAN = 4
end

# ============================================================================
# Material Type
# ============================================================================

"""
	Material{NC,NR,C,RT} <: AbstractMaterial

A complete material definition with components, reactions, mixing rules, and gas transport.

# Type Parameters
- `NC`: Number of components
- `NR`: Number of reactions
- `C`: Component tuple type (preserves heterogeneous component types)
- `RT`: Reaction tuple type (preserves heterogeneous reaction types)

# Fields
- `name::Symbol`: Material identifier
- `components::C`: Material components tuple
- `reactions::RT`: Decomposition reactions tuple
- `conductivity_mixing::MixingRule`: How to mix thermal conductivity
- `β::Float64`: Weighting parameter for WEIGHTED mixing [0-1]
- `gas_transfer_mixing::MixingRule`: How to mix gas transfer coefficient
- `β_λ::Float64`: Weighting parameter for gas transfer mixing
- `permeability_mixing::MixingRule`: How to mix component permeabilities
- `β_κ::Float64`: Weighting parameter for permeability mixing

# Gas Transport
Darcy-Fick gas transport is enabled when any solid/liquid component has κ > 0.
The effective permeability is computed from component values via mixing rules,
similar to thermal conductivity and gas transfer coefficient.

Specify κ on each component:
```julia
charring = Material(
	name = :CharringWood,
	components = (
		SolidComponent(:virgin, ρ=530, c=2000, k=0.12, κ=1e-15),  # tight structure
		SolidComponent(:char, ρ=150, c=1500, k=0.20, κ=1e-12),    # porous char
		GaseousComponent(:gas, M=0.044, c=1800, k=0.03, λ=1e-4),
	),
	reactions = (...),
	permeability_mixing = WEIGHTED,
	β_κ = 0.5,
)
```

# Example
```julia
# Simple PMMA (diffusion-only gas transport, κ=0 on all components)
pmma = Material(
	name = :PMMA,
	components = (
		SolidComponent(:virgin, ρ=1190, c=1420, k=0.21, ε=0.86),
		GaseousComponent(:MMA, M=0.100, c=1500, k=0.02, λ=1e-5)
	),
	reactions = (
		Reaction(:decomposition, 1 => 2, A=8.5e12, E=188e3, h=-870e3),
	),
	conductivity_mixing = WEIGHTED,
	β = 0.5
)
```
"""
struct Material{NC, NR, C <: Tuple, RT <: Tuple} <: AbstractMaterial
	name::Symbol
	components::C
	reactions::RT  # Full tuple type to preserve heterogeneous reaction types

	# Mixing parameters
	conductivity_mixing::MixingRule
	β::Float64                    # weighting for WEIGHTED rule
	gas_transfer_mixing::MixingRule
	β_λ::Float64
	permeability_mixing::MixingRule  # How to mix permeability
	β_κ::Float64                      # weighting for permeability mixing

	# Inner constructor with validation
	function Material{NC, NR, C, RT}(
		name::Symbol,
		components::C,
		reactions::RT,
		conductivity_mixing::MixingRule,
		β::Float64,
		gas_transfer_mixing::MixingRule,
		β_λ::Float64,
		permeability_mixing::MixingRule,
		β_κ::Float64,
	) where {NC, NR, C <: Tuple, RT <: Tuple}
		# Validate mixing parameters
		0 ≤ β ≤ 1 || throw(ArgumentError("β must be in [0,1], got $β"))
		0 ≤ β_λ ≤ 1 || throw(ArgumentError("β_λ must be in [0,1], got $β_λ"))
		0 ≤ β_κ ≤ 1 || throw(ArgumentError("β_κ must be in [0,1], got $β_κ"))
		# Validate component count matches
		length(components) == NC || throw(ArgumentError("Component count mismatch"))
		length(reactions) == NR || throw(ArgumentError("Reaction count mismatch"))
		new{NC, NR, C, RT}(name, components, reactions, conductivity_mixing, β, gas_transfer_mixing, β_λ, permeability_mixing, β_κ)
	end
end

# ============================================================================
# Constructors
# ============================================================================

"""
	Material(; name, components, reactions=(), conductivity_mixing=WEIGHTED, β=0.5, 
			  gas_transfer_mixing=WEIGHTED, β_λ=0.5, permeability_mixing=WEIGHTED, β_κ=0.5)

Create a material with keyword arguments.

# Arguments
- `name::Symbol`: Material identifier
- `components`: Tuple of components
- `reactions`: Tuple of reactions (can be empty tuple `()`)
- `conductivity_mixing=WEIGHTED`: Mixing rule for conductivity
- `β=0.5`: Weighting parameter for conductivity mixing
- `gas_transfer_mixing=WEIGHTED`: Mixing rule for gas transfer
- `β_λ=0.5`: Weighting parameter for gas transfer mixing
- `permeability_mixing=WEIGHTED`: Mixing rule for permeability
- `β_κ=0.5`: Weighting parameter for permeability mixing
"""
function Material(;
	name::Symbol,
	components::C,
	reactions::RT = (),
	conductivity_mixing::MixingRule = WEIGHTED,
	β::Real = 0.5,
	gas_transfer_mixing::MixingRule = WEIGHTED,
	β_λ::Real = 0.5,
	permeability_mixing::MixingRule = WEIGHTED,
	β_κ::Real = 0.5,
) where {C <: Tuple, RT <: Tuple}
	NC = length(components)
	NR = length(reactions)
	Material{NC, NR, C, RT}(
		name,
		components,
		reactions,
		conductivity_mixing,
		Float64(β),
		gas_transfer_mixing,
		Float64(β_λ),
		permeability_mixing,
		Float64(β_κ),
	)
end

# ============================================================================
# Accessor Functions
# ============================================================================

"""
	n_components(m::Material{NC}) where NC

Get the number of components in a material.
"""
n_components(::Material{NC}) where {NC} = NC

"""
	n_reactions(m::Material{NC,NR}) where {NC,NR}

Get the number of reactions in a material.
"""
n_reactions(::Material{NC, NR}) where {NC, NR} = NR

"""
	n_solids(m::Material)

Get the number of solid-phase components.
"""
function n_solids(m::Material)
	count(c -> c.phase == SOLID, m.components)
end

"""
	n_liquids(m::Material)

Get the number of liquid-phase components.
"""
function n_liquids(m::Material)
	count(c -> c.phase == LIQUID, m.components)
end

"""
	n_gases(m::Material)

Get the number of gas-phase components.
"""
function n_gases(m::Material)
	count(c -> c.phase == GAS, m.components)
end

"""
	has_darcy_flow(m::Material) -> Bool

Check if material uses Darcy-Fick gas transport.

Returns true if any solid/liquid component has non-zero permeability (κ > 0).

The component-based approach: define κ on each solid component,
and the effective permeability will be computed via mixing rules.
"""
function has_darcy_flow(m::Material)
	# Check component-based permeability
	for comp in m.components
		if has_permeability(comp)
			# Check if κ > 0 at a reference temperature
			if get_permeability(comp, 500.0) > 0.0
				return true
			end
		end
	end
	return false
end

"""
	get_component(m::Material, i::Int)

Get the i-th component of a material.
"""
@inline get_component(m::Material, i::Int) = m.components[i]

"""
	get_component(m::Material, name::Symbol)

Get a component by name.
"""
function get_component(m::Material, name::Symbol)
	for c in m.components
		c.name == name && return c
	end
	error("Component :$name not found in material $(m.name)")
end

"""
	get_reaction(m::Material, i::Int)

Get the i-th reaction of a material.
"""
@inline get_reaction(m::Material, i::Int) = m.reactions[i]

"""
	get_reaction(m::Material, name::Symbol)

Get a reaction by name.
"""
function get_reaction(m::Material, name::Symbol)
	for r in m.reactions
		r.name == name && return r
	end
	error("Reaction :$name not found in material $(m.name)")
end

"""
	component_index(m::Material, name::Symbol)

Get the index of a component by name.
"""
function component_index(m::Material, name::Symbol)
	for (i, c) in enumerate(m.components)
		c.name == name && return i
	end
	error("Component :$name not found in material $(m.name)")
end

"""
	solid_indices(m::Material)

Get indices of all solid-phase components.
"""
function solid_indices(m::Material{NC}) where {NC}
	indices = Int[]
	for i in 1:NC
		m.components[i].phase == SOLID && push!(indices, i)
	end
	return indices
end

"""
	gas_indices(m::Material)

Get indices of all gas-phase components.
"""
function gas_indices(m::Material{NC}) where {NC}
	indices = Int[]
	for i in 1:NC
		m.components[i].phase == GAS && push!(indices, i)
	end
	return indices
end

# ============================================================================
# Validation
# ============================================================================

"""
	validate_material(m::Material)

Validate material definition for consistency.

Checks:
- All reaction indices refer to valid components
- Stoichiometry is mass-balanced
- At least one component exists
"""
function validate_material(m::Material{NC, NR}) where {NC, NR}
	# Must have at least one component
	NC >= 1 || error("Material must have at least one component")

	# Validate reactions
	for (i, rxn) in enumerate(m.reactions)
		# Check reactant indices
		if rxn.reactant1 != 0
			1 <= rxn.reactant1 <= NC || error("Reaction $i: reactant1 index $(rxn.reactant1) out of range")
		end
		if rxn.reactant2 != 0
			1 <= rxn.reactant2 <= NC || error("Reaction $i: reactant2 index $(rxn.reactant2) out of range")
		end

		# Check product indices
		if rxn.product1 != 0
			1 <= rxn.product1 <= NC || error("Reaction $i: product1 index $(rxn.product1) out of range")
		end
		if rxn.product2 != 0
			1 <= rxn.product2 <= NC || error("Reaction $i: product2 index $(rxn.product2) out of range")
		end

		# Check mass balance
		verify_stoichiometry(rxn) || @warn "Reaction $i ($(rxn.name)) fails mass balance"
	end

	return true
end

# ============================================================================
# Display
# ============================================================================

function Base.show(io::IO, m::Material{NC, NR}) where {NC, NR}
	print(io, "Material(:$(m.name), $NC components, $NR reactions)")
end

function Base.show(io::IO, ::MIME"text/plain", m::Material{NC, NR}) where {NC, NR}
	println(io, "Material: $(m.name)")
	println(io, "  Components ($NC):")
	for (i, c) in enumerate(m.components)
		println(io, "    $i. $(c.name) ($(c.phase))")
	end
	if NR > 0
		println(io, "  Reactions ($NR):")
		for (i, r) in enumerate(m.reactions)
			println(io, "    $i. $(r.name): A=$(r.A), E=$(r.E) J/mol")
		end
	end
	println(io, "  Conductivity mixing: $(m.conductivity_mixing), β=$(m.β)")
end

# ============================================================================
# Generated Functions for Type-Stable Iteration over Heterogeneous Components
# ============================================================================
#
# These functions use @generated to unroll loops at compile time, avoiding
# the type instability that occurs when indexing heterogeneous tuples at runtime.
# This is critical for performance in hot paths like RHS evaluation.

"""
	is_gas_phase(material::Material, j::Int) -> Bool

Check if component j is in gas phase. Uses dispatch for type stability.
"""
@generated function is_gas_phase(material::Material{NC}, ::Val{J}) where {NC, J}
	if J < 1 || J > NC
		return :(error("Component index $J out of range 1:$NC"))
	end
	return :(material.components[$J].phase == GAS)
end

# Convenience method that takes Int (not type-stable, but used in non-critical paths)
is_gas_phase(material::Material, j::Int) = material.components[j].phase == GAS

"""
	is_solid_phase(material::Material, j::Int) -> Bool

Check if component j is in solid phase. Uses dispatch for type stability.
"""
@generated function is_solid_phase(material::Material{NC}, ::Val{J}) where {NC, J}
	if J < 1 || J > NC
		return :(error("Component index $J out of range 1:$NC"))
	end
	return :(material.components[$J].phase == SOLID)
end

# Convenience method
is_solid_phase(material::Material, j::Int) = material.components[j].phase == SOLID

"""
	get_swelling_factor(material::Material, ::Val{J}) -> Float64

Get swelling factor for component J. Type-stable via dispatch.
"""
@generated function get_swelling_factor(material::Material{NC}, ::Val{J}) where {NC, J}
	if J < 1 || J > NC
		return :(error("Component index $J out of range 1:$NC"))
	end
	return :(material.components[$J].swelling_factor)
end

"""
	get_component_density(material::Material, ::Val{J}, T) -> Float64

Get density of component J at temperature T. Type-stable via dispatch.
"""
@generated function get_component_density(material::Material{NC}, ::Val{J}, T) where {NC, J}
	if J < 1 || J > NC
		return :(error("Component index $J out of range 1:$NC"))
	end
	return :(material.components[$J].density(T))
end

"""
	foreach_component(f, material::Material{NC}) where NC

Apply function f to each component index, using dispatch for type stability.
f should take (material, Val{j}) as arguments.

Example:
	foreach_component(material) do mat, j_val
		if is_gas_phase(mat, j_val)
			# do something
		end
	end
"""
@generated function foreach_component(f::F, material::Material{NC}) where {F, NC}
	exprs = [:(f(material, Val{$j}())) for j in 1:NC]
	return Expr(:block, exprs...)
end

"""
	foreach_gas_component(f, material::Material{NC}) where NC

Apply function f only to gas component indices. Type-stable.
f should take (material, component_index::Int) as arguments.
"""
@generated function foreach_gas_component(f::F, material::Material{NC}) where {F, NC}
	# Build list of gas components at compile time
	exprs = Expr[]
	for j in 1:NC
		push!(exprs, :(
			if material.components[$j].phase == GAS
				;f(material, $j);
			end
		))
	end
	return Expr(:block, exprs...)
end

"""
	gas_component_indices(::Material{NC,NR,C,RT}) where {NC,NR,C,RT}

Return a tuple of gas component indices at compile time. Use this instead of
building an array at runtime.
"""
@generated function gas_component_indices(material::Material{NC}) where NC
	# This generates the expression to build a tuple of gas indices
	# We need to evaluate at generated function time, which is tricky
	# Instead, return code that builds the right tuple
	return quote
		_collect_gas_indices(material.components)
	end
end

# Helper that recursively collects gas indices - unrolled at compile time
@inline _collect_gas_indices(::Tuple{}) = ()
@inline function _collect_gas_indices(comps::Tuple)
	first_comp = comps[1]
	rest = Base.tail(comps)
	rest_indices = _collect_gas_indices(rest)
	if first_comp.phase == GAS
		return (1, map(x -> x + 1, rest_indices)...)
	else
		return map(x -> x + 1, rest_indices)
	end
end

# ============================================================================
# Unrolled Component Loops for Hot Paths
# ============================================================================

"""
	_unroll_gas_flux!(J::Matrix, j::Val{J}, ξ_L, ξ_R, f, flux_fn, material) where J

Unrolled gas flux computation for a single component.
If component J is GAS, compute flux and store in J[J, f].
Uses the rigorous λρ_face formulation.
"""
@generated function _unroll_gas_flux!(J::Matrix, ::Val{IDX}, ξ_L::NTuple{NC}, ξ_R::NTuple{NC},
	f::Int, λρ_face, ρ_L, ρ_R, d_LR, material::Material{NC}) where {IDX, NC}
	if IDX < 1 || IDX > NC
		return :(error("Component index $IDX out of range 1:$NC"))
	end
	return quote
		@inbounds if material.components[$IDX].phase == GAS
			ξ_L_j = ξ_L[$IDX]
			ξ_R_j = ξ_R[$IDX]
			# Rigorous formulation: J = -λρ_face × ΔY / d
			J[$IDX, f] = -λρ_face * ((ξ_R_j / ρ_R) - (ξ_L_j / ρ_L)) / d_LR
		end
	end
end

"""
	compute_gas_fluxes_unrolled!(J, f, ξ_L, ξ_R, λρ_face, ρ_L, ρ_R, d_LR, material::Material{NC})

Compute gas fluxes for all components at face f, using compile-time unrolling.

Uses the rigorous λρ_face formulation where λρ_face is the harmonic mean of (λ·ρ):
	J = -λρ_face × (Y_R - Y_L) / d_LR

This ensures physically correct flux at interfaces with sharp property jumps.
"""
@generated function compute_gas_fluxes_unrolled!(J::Matrix, f::Int,
	ξ_L::NTuple{NC}, ξ_R::NTuple{NC}, λρ_face, ρ_L, ρ_R, d_LR,
	material::Material{NC}) where NC

	exprs = Expr[]
	for j in 1:NC
		push!(exprs, quote
			@inbounds if material.components[$j].phase == GAS
				ξ_L_j = ξ_L[$j]
				ξ_R_j = ξ_R[$j]
				# Rigorous formulation: J = -λρ_face × ΔY / d
				J[$j, f] = -λρ_face * ((ξ_R_j / ρ_R) - (ξ_L_j / ρ_L)) / d_LR
			end
		end)
	end
	return Expr(:block, exprs...)
end

"""
	compute_gas_fluxes_species_scaled_unrolled!(J::Matrix, f::Int,
		ξ_L::NTuple{NC}, ξ_R::NTuple{NC}, λρ_face, ρ_L, ρ_R, d_LR,
		material::Material{NC})

Compute gas fluxes with species-specific molecular weight scaling.

Same as `compute_gas_fluxes_unrolled!` but scales diffusive flux by √(M_ref/Mⱼ):
- Light gases (M < M_ref) diffuse FASTER
- Heavy gases (M > M_ref) diffuse SLOWER

The reference molar mass M_ref = 0.029 kg/mol (air).
"""
@generated function compute_gas_fluxes_species_scaled_unrolled!(J::Matrix, f::Int,
	ξ_L::NTuple{NC}, ξ_R::NTuple{NC}, λρ_face, ρ_L, ρ_R, d_LR,
	material::Material{NC}) where NC

	M_REF = 0.029  # kg/mol (air)
	exprs = Expr[]
	for j in 1:NC
		push!(exprs, quote
			@inbounds if material.components[$j].phase == GAS
				ξ_L_j = ξ_L[$j]
				ξ_R_j = ξ_R[$j]

				# Species-specific scaling: scale_j = √(M_ref / Mⱼ)
				M_j = material.components[$j].molar_mass
				scale_j = M_j > 1e-6 ? sqrt($M_REF / M_j) : 1.0

				# Rigorous formulation: J = -λρ_face × scale_j × ΔY / d
				λρ_scaled = λρ_face * scale_j
				J[$j, f] = -λρ_scaled * ((ξ_R_j / ρ_R) - (ξ_L_j / ρ_L)) / d_LR
			end
		end)
	end
	return Expr(:block, exprs...)
end

"""
	compute_gas_fluxes_darcy_fick_unrolled!(J::Matrix, f::Int, v_g::Float64,
		ξ_L::NTuple{NC}, ξ_R::NTuple{NC}, λρ_face, ρ_L, ρ_R, d_LR,
		material::Material{NC})

Compile-time unrolled Darcy-Fick flux computation for type stability.

Computes J = J_adv + J_diff for each gas component:
- J_adv = v_g × ξ_upwind (first-order upwind)
- J_diff = -λρ_face × ΔY / d_LR (Fickian diffusion)
"""
@generated function compute_gas_fluxes_darcy_fick_unrolled!(J::Matrix, f::Int, v_g,
	ξ_L::NTuple{NC}, ξ_R::NTuple{NC}, λρ_face, ρ_L, ρ_R, d_LR,
	material::Material{NC}) where NC

	exprs = Expr[]
	for j in 1:NC
		push!(exprs, quote
			@inbounds if material.components[$j].phase == GAS
				ξ_L_j = ξ_L[$j]
				ξ_R_j = ξ_R[$j]

				# Diffusive flux: J_diff = -λρ_face × ΔY / d
				J_diff = -λρ_face * ((ξ_R_j / ρ_R) - (ξ_L_j / ρ_L)) / d_LR

				# Advective flux: J_adv = v_g × ξ_upwind (first-order upwind)
				ξ_upwind = v_g >= 0.0 ? ξ_L_j : ξ_R_j
				J_adv = v_g * ξ_upwind

				# Combined Darcy-Fick flux
				J[$j, f] = J_adv + J_diff
			end
		end)
	end
	return Expr(:block, exprs...)
end

"""
	compute_gas_fluxes_darcy_fick_species_scaled_unrolled!(J::Matrix, f::Int, v_g,
		ξ_L::NTuple{NC}, ξ_R::NTuple{NC}, λρ_face, ρ_L, ρ_R, d_LR,
		material::Material{NC})

Darcy-Fick flux with species-specific molecular weight scaling.

Same as `compute_gas_fluxes_darcy_fick_unrolled!` but scales diffusive flux by √(M_ref/Mⱼ).
Note: Advective flux is NOT scaled - all species move with the bulk gas velocity.
"""
@generated function compute_gas_fluxes_darcy_fick_species_scaled_unrolled!(J::Matrix, f::Int, v_g,
	ξ_L::NTuple{NC}, ξ_R::NTuple{NC}, λρ_face, ρ_L, ρ_R, d_LR,
	material::Material{NC}) where NC

	M_REF = 0.029  # kg/mol (air)
	exprs = Expr[]
	for j in 1:NC
		push!(exprs, quote
			@inbounds if material.components[$j].phase == GAS
				ξ_L_j = ξ_L[$j]
				ξ_R_j = ξ_R[$j]

				# Species-specific scaling for diffusion only
				M_j = material.components[$j].molar_mass
				scale_j = M_j > 1e-6 ? sqrt($M_REF / M_j) : 1.0

				# Diffusive flux: J_diff = -λρ_face × scale_j × ΔY / d
				λρ_scaled = λρ_face * scale_j
				J_diff = -λρ_scaled * ((ξ_R_j / ρ_R) - (ξ_L_j / ρ_L)) / d_LR

				# Advective flux: J_adv = v_g × ξ_upwind (NOT scaled)
				ξ_upwind = v_g >= 0.0 ? ξ_L_j : ξ_R_j
				J_adv = v_g * ξ_upwind

				# Combined Darcy-Fick flux
				J[$j, f] = J_adv + J_diff
			end
		end)
	end
	return Expr(:block, exprs...)
end

"""
	apply_species_residual_unrolled!(dξ::Matrix, i::Int, div_J::Matrix, dξ_rxn::Matrix, material::Material{NC})

Apply species residual: for gases dξ = -div_J + dξ_rxn, for solids dξ = dξ_rxn.
Uses compile-time unrolling for type stability.
"""
@generated function apply_species_residual_unrolled!(dξ::Matrix, i::Int,
	div_J::Matrix, dξ_rxn::Matrix, material::Material{NC}) where NC

	exprs = Expr[]
	for j in 1:NC
		push!(exprs, quote
			@inbounds if material.components[$j].phase == GAS
				dξ[$j, i] = -div_J[$j, i] + dξ_rxn[$j, i]
			else
				dξ[$j, i] = dξ_rxn[$j, i]
			end
		end)
	end
	return Expr(:block, exprs...)
end
