"""
Reaction type definitions.

Reactions define the kinetics of thermal decomposition, including
stoichiometry, Arrhenius parameters, and heat of reaction.

Reference: Design Document Section 6.2
"""

# ============================================================================
# Reaction Type
# ============================================================================

"""
	Reaction{Ph<:AbstractPropertyFunction} <: AbstractReaction

A single decomposition reaction with Arrhenius kinetics.

The reaction rate is computed as:
```
r = A * exp(-E/RT) * ξ₁^n₁ * ξ₂^n₂ * H(T - T_min) * H(T_max - T)
```

# Type Parameters
- `Ph`: Type of heat of reaction property function

# Fields
- `name::Symbol`: Reaction identifier
- `reactant1::Int`: Index of first reactant (0 = none)
- `reactant2::Int`: Index of second reactant (0 = none)
- `order1::Float64`: Reaction order for reactant 1
- `order2::Float64`: Reaction order for reactant 2
- `product1::Int`: Index of first product (0 = none)
- `product2::Int`: Index of second product (0 = none)
- `stoich::NTuple{4,Float64}`: Mass stoichiometry (ν₁, ν₂, ν₃, ν₄)
- `A::Float64`: Pre-exponential factor [s⁻¹] or [m³/(kg·s)]
- `E::Float64`: Activation energy [J/mol]
- `heat::Ph`: Heat of reaction [J/kg of reactant1]
- `T_min::Float64`: Minimum temperature for reaction [K]
- `T_max::Float64`: Maximum temperature for reaction [K]

# Stoichiometry Convention
The stoichiometry tuple (ν₁, ν₂, ν₃, ν₄) represents:
- ν₁: Mass of reactant1 consumed per unit reaction extent
- ν₂: Mass of reactant2 consumed (0 if single reactant)
- ν₃: Mass of product1 produced
- ν₄: Mass of product2 produced

Mass balance: ν₁ + ν₂ = ν₃ + ν₄

# Heat of Reaction Convention
- h > 0: Endothermic (absorbs heat, cools material)
- h < 0: Exothermic (releases heat, heats material)

# Heat of Reaction Definition

The heat of reaction `h` in this code follows the convention:

    h = H_products - H_reactants

where H is total enthalpy (formation + sensible). This means:

- For endothermic decomposition (e.g., pyrolysis): h > 0
- For exothermic reactions: h < 0

## Important: Sensible Enthalpy Inclusion

If your heat of reaction data is at a reference temperature T_ref (e.g., 298 K),
you may need to include a sensible enthalpy correction:

    h_actual(T) = h_ref + ∫_{T_ref}^{T} [Σ(ν_prod × cp_prod) - Σ(ν_react × cp_react)] dT

For temperature-dependent heats, use `StateDependentProperty` or `PolynomialProperty`.

## Gas Enthalpy at Boundaries

With the advection form of energy equation (∂T/∂t + w·∇T = ...), the enthalpy
of gas leaving the domain is automatically handled. No explicit boundary
enthalpy flux is needed PROVIDED that:

1. Heat of reaction includes sensible enthalpy of gas products
2. Gas heat capacity is properly accounted in effective ρcₚ

If using heats of reaction at 298 K without correction, be aware that the
sensible enthalpy of escaping gas may need to be explicitly tracked in the
energy balance.

# Example
```julia
# Single-step decomposition: virgin -> gas (endothermic)
decomp = Reaction(
	name = :decomposition,
	reactant1 = 1,      # virgin (component index)
	reactant2 = 0,      # no second reactant
	order1 = 1.0,
	order2 = 0.0,
	product1 = 2,       # gas (component index)
	product2 = 0,       # no second product
	stoich = (1.0, 0.0, 1.0, 0.0),
	A = 8.5e12,
	E = 188e3,
	heat = ConstantProperty(870e3),  # positive = endothermic
	T_min = 0.0,
	T_max = Inf
)
```
"""
struct Reaction{Ph <: AbstractPropertyFunction} <: AbstractReaction
	name::Symbol

	# Reactants (indices into component array, 0 = none)
	reactant1::Int
	reactant2::Int
	order1::Float64      # reaction order for reactant 1
	order2::Float64      # reaction order for reactant 2

	# Products (indices, 0 = none)
	product1::Int
	product2::Int

	# Stoichiometry: mass of product per mass of reactant consumed
	# ν₁ kg reactant1 + ν₂ kg reactant2 → ν₃ kg product1 + ν₄ kg product2
	stoich::NTuple{4, Float64}

	# Kinetic parameters
	A::Float64           # pre-exponential [s⁻¹] or [m³/(kg·s)]
	E::Float64           # activation energy [J/mol]

	# Heat of reaction (T-dependent)
	heat::Ph             # [J/kg of reactant1 consumed]

	# Temperature limits
	T_min::Float64       # K
	T_max::Float64       # K
end

# ============================================================================
# Validation Helpers
# ============================================================================

"""
	_validate_stoichiometry(stoich::NTuple{4,Float64}, name::Symbol)

Validate that stoichiometric coefficients conserve mass.

Throws an error if mass balance is violated: ν₁ + ν₂ ≠ ν₃ + ν₄

# Arguments
- `stoich`: Stoichiometric coefficients (ν_reactant1, ν_reactant2, ν_product1, ν_product2)
- `name`: Reaction name for error reporting

# Mass Balance
For mass conservation, total reactant mass must equal total product mass:
    ν₁ + ν₂ = ν₃ + ν₄

Tolerance: 1e-10 (relative to sum of coefficients)
"""
function _validate_stoichiometry(stoich::NTuple{4,Float64}, name::Symbol)
	reactant_mass = stoich[1] + stoich[2]
	product_mass = stoich[3] + stoich[4]

	# Use relative tolerance based on magnitude
	total_mass = reactant_mass + product_mass
	tol = max(1e-10 * total_mass, 1e-14)

	if !isapprox(reactant_mass, product_mass; atol=tol)
		error("""
			Reaction mass balance violated for reaction :$name
			  Reactant mass: $(stoich[1]) + $(stoich[2]) = $reactant_mass
			  Product mass:  $(stoich[3]) + $(stoich[4]) = $product_mass
			  Difference: $(reactant_mass - product_mass)

			Set `validate_mass=false` to disable this check (not recommended).
			""")
	end

	return nothing
end

# ============================================================================
# Convenience Constructors
# ============================================================================

"""
	Reaction(name, reactant => product; A, E, h, n=1.0, T_min=0.0, T_max=Inf, validate_mass=true)

Create a simple single-reactant, single-product reaction.

# Arguments
- `name::Symbol`: Reaction identifier
- `reactant => product`: Pair of component indices
- `A`: Pre-exponential factor [s⁻¹]
- `E`: Activation energy [J/mol]
- `h`: Heat of reaction [J/kg] (number or property function)
- `n=1.0`: Reaction order
- `T_min=0.0`: Minimum temperature [K]
- `T_max=Inf`: Maximum temperature [K]
- `validate_mass=true`: If true, enforce mass balance in stoichiometry

# Example
```julia
rxn = Reaction(:decomp, 1 => 2, A=8.5e12, E=188e3, h=870e3)  # positive h = endothermic
```
"""
function Reaction(name::Symbol, pair::Pair{Int, Int};
	A::Real, E::Real, h,
	n::Real = 1.0,
	T_min::Real = 0.0,
	T_max::Real = Inf,
	validate_mass::Bool = true)
	reactant, product = pair
	stoich = (1.0, 0.0, 1.0, 0.0)

	# Validate mass balance if requested
	if validate_mass
		_validate_stoichiometry(stoich, name)
	end

	Reaction(
		name,
		reactant, 0,
		Float64(n), 0.0,
		product, 0,
		stoich,
		Float64(A),
		Float64(E),
		to_property(h),
		Float64(T_min),
		Float64(T_max),
	)
end

"""
	Reaction(name, reactant => (product1, product2); A, E, h, yields, n=1.0, T_min=0.0, T_max=Inf, validate_mass=true)

Create a reaction with two products.

# Arguments
- `name::Symbol`: Reaction identifier
- `reactant => (product1, product2)`: Reactant and product indices
- `A`: Pre-exponential factor [s⁻¹]
- `E`: Activation energy [J/mol]
- `h`: Heat of reaction [J/kg]
- `yields`: Tuple (y₁, y₂) of mass yields for each product
- `n=1.0`: Reaction order
- `T_min=0.0`: Minimum temperature [K]
- `T_max=Inf`: Maximum temperature [K]
- `validate_mass=true`: If true, enforce mass balance in stoichiometry

# Example
```julia
# Charring reaction: virgin -> 0.2 char + 0.8 gas (endothermic)
rxn = Reaction(:charring, 1 => (2, 3), A=1e15, E=200e3, h=500e3, yields=(0.2, 0.8))
```
"""
function Reaction(name::Symbol, pair::Pair{Int, Tuple{Int, Int}};
	A::Real, E::Real, h,
	yields::Tuple{Real, Real},
	n::Real = 1.0,
	T_min::Real = 0.0,
	T_max::Real = Inf,
	validate_mass::Bool = true)
	reactant = pair.first
	product1, product2 = pair.second
	y1, y2 = yields
	stoich = (1.0, 0.0, Float64(y1), Float64(y2))

	# Validate mass balance if requested
	if validate_mass
		_validate_stoichiometry(stoich, name)
	end

	Reaction(
		name,
		reactant, 0,
		Float64(n), 0.0,
		product1, product2,
		stoich,
		Float64(A),
		Float64(E),
		to_property(h),
		Float64(T_min),
		Float64(T_max),
	)
end

# ============================================================================
# ReactionSet Type
# ============================================================================

"""
	ReactionSet{NR,Ph}

A collection of reactions for a material.

# Type Parameters
- `NR`: Number of reactions
- `Ph`: Property function type for heat of reaction

# Fields
- `reactions::NTuple{NR,Reaction{Ph}}`: Tuple of reactions
"""
struct ReactionSet{NR, Ph <: AbstractPropertyFunction} <: AbstractReaction
	reactions::NTuple{NR, Reaction{Ph}}
end

# Convenience constructor from varargs
ReactionSet(rxns::Reaction...) = ReactionSet(rxns)

# ============================================================================
# Accessor Functions
# ============================================================================

"""
	n_reactions(rs::ReactionSet{NR}) where NR

Get the number of reactions in a reaction set.
"""
n_reactions(::ReactionSet{NR}) where {NR} = NR

"""
	get_reaction(rs::ReactionSet, i::Int)

Get the i-th reaction from a reaction set.
"""
@inline get_reaction(rs::ReactionSet, i::Int) = rs.reactions[i]

# ============================================================================
# Reaction Rate Computation
# ============================================================================

"""
	rate_constant(rxn::Reaction, T::Real)

Compute the Arrhenius rate constant k(T) = A * exp(-E/RT).
"""
@inline function rate_constant(rxn::Reaction, T::Real)
	rxn.A * exp(-rxn.E / (R_GAS * T))
end

"""
	is_active(rxn::Reaction, T::Real)

Check if reaction is active at temperature T.
"""
@inline function is_active(rxn::Reaction, T::Real)
	T >= rxn.T_min && T <= rxn.T_max
end

"""
	get_heat(rxn::Reaction, T::Real)

Get heat of reaction at temperature T.
"""
@inline get_heat(rxn::Reaction, T::Real) = rxn.heat(T)

# ============================================================================
# Stoichiometry Helpers
# ============================================================================

"""
	mass_balance_error(rxn::Reaction)

Compute the mass balance error for a reaction (should be ≈ 0).
"""
function mass_balance_error(rxn::Reaction)
	ν = rxn.stoich
	reactant_mass = ν[1] + ν[2]
	product_mass = ν[3] + ν[4]
	return abs(reactant_mass - product_mass)
end

"""
	verify_stoichiometry(rxn::Reaction; tol=1e-10)

Verify that reaction stoichiometry satisfies mass balance.
"""
function verify_stoichiometry(rxn::Reaction; tol::Real = 1e-10)
	err = mass_balance_error(rxn)
	if err > tol
		@warn "Reaction $(rxn.name) has mass balance error: $err"
		return false
	end
	return true
end

"""
	check_material_mass_balance(reactions::Tuple; tol=1e-10) -> Bool

Verify all reactions in a material conserve mass.

Returns `true` if all reactions are mass-balanced within tolerance.
Prints warnings for any imbalanced reactions.

# Arguments
- `reactions`: Tuple of reactions (from Material.reactions)
- `tol`: Absolute tolerance for mass balance error

# Example
```julia
material = Material(...)
if !check_material_mass_balance(material.reactions)
    @warn "Material has reactions with mass balance issues"
end
```
"""
function check_material_mass_balance(reactions::Tuple; tol::Real = 1e-10)
	all_balanced = true

	for (i, rxn) in enumerate(reactions)
		if !verify_stoichiometry(rxn; tol=tol)
			all_balanced = false
		end
	end

	return all_balanced
end
