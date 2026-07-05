"""
Radiation model types and enumerations.

This file defines the fundamental types for volumetric radiation transport
using the P1 approximation. These types are used throughout the radiation
implementation.

Reference: P1_Volumetric_Reradiation_Implementation.md Section 4.2
"""

# ============================================================================
# Radiation Model Selection
# ============================================================================

"""
	RadiationModel

Enumeration of available radiation transport models for pyrolysis simulations.

# Values

- `NO_RADIATION`: No volumetric radiation absorption. External radiative boundary
  conditions (RadiativeBC, RadiativeFluxBC) still apply at surfaces. Use for
  transparent materials or when radiation effects are negligible.

- `SURFACE_ABSORPTION`: All incident radiation is absorbed in the surface cell.
  Computationally efficient O(1) with tridiagonal Jacobian. Appropriate for opaque
  materials where optical thickness τ >> 1 in the first cell. Component absorption
  coefficients (α) are ignored - a warning is issued if they are set.

- `BEER_LAMBERT`: In-depth volumetric absorption following Beer-Lambert law:
  I(z) = I₀ exp(-∫ α_eff dz). Radiation penetrates into the material based on
  component absorption coefficients. Creates non-local coupling in Jacobian
  (upper triangular pattern). Use for semi-transparent materials (polymers, foams).

- `P1_QUASI_STEADY`: P1 approximation with quasi-steady radiation field (∂G/∂t = 0).
  Solves for incident radiation intensity G including volumetric emission (4σαT⁴).
  Recommended when internal re-radiation is significant (T > 800K, optically thick).
  Uses Newton iteration to solve radiation field each timestep.

# Physics

## Beer-Lambert Model
Exponential attenuation through absorbing medium:
  I(z) = I₀ × exp(-∫ α_eff(z') dz')
  q'''(z) = α_eff(z) × I(z)

## P1 Approximation
Solves for incident radiation intensity G [W/m²]:
  ∇·(D_rad ∇G) = α G - 4σαT⁴
where D_rad = 1/(3α) is the radiation diffusion coefficient.

# Validity (Optical Thickness τ = αL)
- τ > 10: Optically thick - P1 excellent, surface absorption acceptable
- 1 < τ < 10: Intermediate - P1 good, Beer-Lambert recommended
- 0.1 < τ < 1: Optically thin - Beer-Lambert recommended
- τ < 0.1: Transparent - Use NO_RADIATION

For typical pyrolysis materials (τ ≈ 5-50), P1_QUASI_STEADY is ideal for
high-temperature re-radiation, BEER_LAMBERT for semi-transparent polymers.

# Example
```julia
# Surface absorption for opaque material
solution = solve(problem; radiation_model = SURFACE_ABSORPTION)

# In-depth absorption for semi-transparent PMMA
solution = solve(problem; radiation_model = BEER_LAMBERT)

# P1 for high-temperature char with re-radiation
solution = solve(problem; radiation_model = P1_QUASI_STEADY)
```

# Performance Characteristics
| Model              | Loop Cost | Jacobian Pattern   | Memory    |
|--------------------|-----------|-------------------|-----------|
| NO_RADIATION       | O(1)      | Tridiagonal       | Baseline  |
| SURFACE_ABSORPTION | O(1)      | Tridiagonal       | Baseline  |
| BEER_LAMBERT       | O(n)      | Upper triangular  | +O(n²)    |
| P1_QUASI_STEADY    | O(10n)    | Tridiagonal       | +O(n)     |
"""
@enum RadiationModel begin
    NO_RADIATION = 0
    SURFACE_ABSORPTION = 1
    BEER_LAMBERT = 2
    P1_QUASI_STEADY = 3
end

# ============================================================================
# Radiation Model Query Functions
# ============================================================================

"""
	requires_absorption_coefficients(model::RadiationModel) -> Bool

Returns true if the radiation model uses component absorption coefficients (α).

SURFACE_ABSORPTION does not use α values - radiation is deposited directly in
the surface cell. BEER_LAMBERT and P1_QUASI_STEADY require α to compute
optical depth and radiation diffusion coefficient.

# Example
```julia
if requires_absorption_coefficients(radiation_model)
	# Validate that material has α values set
end
```
"""
requires_absorption_coefficients(model::RadiationModel) =
    model == BEER_LAMBERT || model == P1_QUASI_STEADY

"""
	has_nonlocal_jacobian_coupling(model::RadiationModel) -> Bool

Returns true if the radiation model creates non-local (non-tridiagonal) Jacobian coupling.

BEER_LAMBERT creates non-local coupling because the radiation intensity at cell i
depends on the optical depth through all cells from i to the surface. This results
in an upper triangular sparsity pattern in the Jacobian.

All other models maintain tridiagonal structure.

# Example
```julia
if has_nonlocal_jacobian_coupling(radiation_model)
	# Use upper triangular sparsity pattern
else
	# Use tridiagonal sparsity pattern
end
```
"""
has_nonlocal_jacobian_coupling(model::RadiationModel) = model == BEER_LAMBERT

"""
	requires_radiation_cache(model::RadiationModel) -> Bool

Returns true if the radiation model requires a RadiationCache for P1 solver state.

Only P1_QUASI_STEADY needs the radiation cache to store incident radiation field G,
diffusion coefficients, and Newton iteration work arrays.

# Example
```julia
if requires_radiation_cache(radiation_model)
	rad_cache = RadiationCache(n_cells)
end
```
"""
requires_radiation_cache(model::RadiationModel) = model == P1_QUASI_STEADY

# ============================================================================
# Abstract Type for Radiation Boundary Conditions
# ============================================================================

"""
	AbstractRadiationBC <: AbstractBoundaryCondition

Abstract supertype for radiation boundary conditions used with P1 models.

Radiation BCs are separate from thermal BCs because they apply to the 
radiation intensity field G rather than temperature T. The energy exchange
at boundaries is handled through the Marshak condition which couples G 
and T at the surface.

# Concrete Types
- `MarshakRadiationBC`: Diffuse surface with external radiation source
- `ReflectiveRadiationBC`: Zero gradient (adiabatic for radiation)  
- `DirichletRadiationBC`: Prescribed incident radiation intensity

# See Also
- `MarshakRadiationBC`, `ReflectiveRadiationBC`, `DirichletRadiationBC`
- `AbstractThermalBC` for temperature boundary conditions
"""
abstract type AbstractRadiationBC <: AbstractBoundaryCondition end
