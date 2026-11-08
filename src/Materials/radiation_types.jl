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

- `RADIATIVE_TRANSFER`: Radiative transfer through the solid, for
  semi-transparent materials (polymers, foams). Every cell
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
- τ < 10: radiative transfer

# Example
```julia
# Surface absorption for opaque material
solution = solve(problem; radiation_model = SURFACE_ABSORPTION)

# Absorption and emission through semi-transparent PMMA
solution = solve(problem; radiation_model = RADIATIVE_TRANSFER)
```
"""
@enum RadiationModel begin
    SURFACE_ABSORPTION = 1
    RADIATIVE_TRANSFER = 2
end
