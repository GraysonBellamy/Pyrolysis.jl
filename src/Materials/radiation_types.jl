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

- `RADIATIVE_TRANSFER`: Radiative transfer through the solid. Every cell
  absorbs and emits at its absorption coefficient and temperature; the
  radiation field is carried by two hemispherical fluxes per face (gray
  two-stream, diffusivity 3/2, integrated exactly across each cell). Incoming
  radiation is a property of the boundary: the exposed face is a gray
  interface at the emissivity of its cell, open to the radiative members of
  its condition, and the back is a mirror unless its condition has a
  radiative member. The surface energy balance keeps convection, prescribed
  fluxes and blowing; radiation leaves it. Planar lines only.

# Validity (optical thickness τ = αL)
- τ > 10: optically thick; surface absorption is adequate
- τ < 10: Beer–Lambert

# Example
```julia
# Surface absorption for opaque material
solution = solve(problem; radiation_model = SURFACE_ABSORPTION)

# In-depth absorption for semi-transparent PMMA
solution = solve(problem; radiation_model = BEER_LAMBERT)

# Absorption and emission through the solid
solution = solve(problem; radiation_model = RADIATIVE_TRANSFER)
```
"""
@enum RadiationModel begin
    SURFACE_ABSORPTION = 1
    BEER_LAMBERT = 2
    RADIATIVE_TRANSFER = 3
end
