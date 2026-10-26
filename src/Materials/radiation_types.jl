"""
The radiation model enumeration.
"""

# ============================================================================
# Radiation Model Selection
# ============================================================================

"""
	RadiationModel

Where incident radiation is absorbed.

# Values

- `NO_RADIATION`: No volumetric radiation absorption. External radiative boundary
  conditions (RadiativeBC, RadiativeFluxBC) still apply at surfaces. Use for
  transparent materials or when radiation effects are negligible.

- `SURFACE_ABSORPTION`: All incident radiation is absorbed at the exposed
  surface. Appropriate for opaque materials, where the optical thickness of the
  first cell is large. Component absorption coefficients (α) are not read.

- `BEER_LAMBERT`: In-depth volumetric absorption following the Beer–Lambert law,
  `I(z) = I₀ exp(−∫ α_eff dz)` with `q‴(z) = α_eff(z) I(z)`: radiation penetrates
  the material according to the component absorption coefficients. Use for
  semi-transparent materials (polymers, foams).

# Validity (optical thickness τ = αL)
- τ > 10: optically thick; surface absorption is adequate
- 0.1 < τ < 10: Beer–Lambert
- τ < 0.1: transparent; `NO_RADIATION`

# Example
```julia
# Surface absorption for opaque material
solution = solve(problem; radiation_model = SURFACE_ABSORPTION)

# In-depth absorption for semi-transparent PMMA
solution = solve(problem; radiation_model = BEER_LAMBERT)
```
"""
@enum RadiationModel begin
    NO_RADIATION = 0
    SURFACE_ABSORPTION = 1
    BEER_LAMBERT = 2
end
