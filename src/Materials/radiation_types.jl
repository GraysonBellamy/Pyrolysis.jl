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

- `SURFACE_ABSORPTION`: All incident radiation is absorbed at the exposed
  surface, through the surface energy balance. Appropriate for opaque materials,
  where the optical thickness of the first cell is large. Component absorption
  coefficients (α) are not read.

- `BEER_LAMBERT`: In-depth volumetric absorption following the Beer–Lambert law,
  `I(z) = I₀ exp(−∫ α_eff dz)` with `q‴(z) = α_eff(z) I(z)`: radiation penetrates
  the material according to the component absorption coefficients. Use for
  semi-transparent materials (polymers, foams).

# Validity (optical thickness τ = αL)
- τ > 10: optically thick; surface absorption is adequate
- τ < 10: Beer–Lambert

# Example
```julia
# Surface absorption for opaque material
solution = solve(problem; radiation_model = SURFACE_ABSORPTION)

# In-depth absorption for semi-transparent PMMA
solution = solve(problem; radiation_model = BEER_LAMBERT)
```
"""
@enum RadiationModel begin
    SURFACE_ABSORPTION = 1
    BEER_LAMBERT = 2
end
