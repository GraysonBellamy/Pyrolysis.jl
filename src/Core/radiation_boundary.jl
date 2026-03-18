"""
Radiation boundary condition types for P1 radiation transport.

These boundary conditions apply to the incident radiation intensity field G
solved by the P1 approximation. They work alongside thermal BCs to couple
radiation and temperature at domain boundaries.

Reference: P1_Volumetric_Reradiation_Implementation.md Section 2.3, 4.2
"""

# ============================================================================
# Marshak Boundary Condition (Primary BC for Exposed Surfaces)
# ============================================================================

"""
	MarshakRadiationBC{FI, FT, Fε} <: AbstractRadiationBC

Marshak boundary condition for P1 radiation at a diffuse surface.

Models radiation exchange between a surface and its surroundings, including
external radiation sources (e.g., cone heater, fire). This is the primary
BC for exposed surfaces in pyrolysis simulations.

# Fields
- `I_external::FI`: External incident radiation flux [W/m²]. Can be constant
  or time-dependent function `I_external(t)`. This represents radiation from
  external sources like cone heaters, radiant panels, or flames.
  
- `T_ambient::FT`: Ambient temperature [K] for radiation exchange with 
  surroundings. Can be constant or time-dependent `T_ambient(t)`.
  
- `emissivity::Fε`: Surface emissivity ε ∈ [0,1]. Can be constant or 
  temperature-dependent `emissivity(T)`. For pyrolysis materials:
  - Virgin polymers: 0.85-0.95
  - Char/carbon: 0.90-0.98

# Physical Model
The Marshak condition for a diffuse surface relates incident radiation G
to external sources and surface emission:

	G_surf - (2/3α) × (∂G/∂n) = 2(I_external + εσT_amb⁴ - εσT_surf⁴)

This accounts for:
- External radiation input (I_external)
- Ambient radiation exchange (εσT_amb⁴)  
- Surface re-radiation loss (εσT_surf⁴)

The net flux into the domain is:
	q_rad = I_external + εσT_amb⁴ - εσT_surf⁴ - G_surf/2

# Examples
```julia
# Cone calorimeter with 25 kW/m² heater (constant)
bc = MarshakRadiationBC(
	I_external = 25e3,      # W/m²
	T_ambient = 300.0,      # K
	emissivity = 0.9
)

# Time-varying external flux (ramp-up)
bc = MarshakRadiationBC(
	I_external = t -> 50e3 * (1 - exp(-t/10)),
	T_ambient = 300.0,
	emissivity = 0.9
)

# Temperature-dependent emissivity
bc = MarshakRadiationBC(
	I_external = 35e3,
	T_ambient = 300.0,
	emissivity = T -> 0.85 + 0.0001 * (T - 300)  # Increases with T
)
```

# See Also
- `ReflectiveRadiationBC` for back surfaces
- `RadiativeBC` (thermal BC) for surface-only radiation without P1
"""
struct MarshakRadiationBC{FI, FT, Fε} <: AbstractRadiationBC
	I_external::FI    # External radiation flux [W/m²]
	T_ambient::FT     # Ambient temperature [K]
	emissivity::Fε    # Surface emissivity [-]
end

# Keyword constructor with defaults
function MarshakRadiationBC(;
	I_external,
	T_ambient = 300.0,
	emissivity = 0.9,
	validate::Bool = true,
)
	if validate
		I_val = _get_scalar_value(I_external)
		T_val = _get_scalar_value(T_ambient)
		ε_val = _get_scalar_value(emissivity)

		if I_val !== nothing && I_val < 0
			throw(ArgumentError("External radiation I_external must be non-negative, got $I_val"))
		end
		if T_val !== nothing && T_val <= 0
			throw(ArgumentError("Ambient temperature must be positive (in Kelvin), got $T_val"))
		end
		if ε_val !== nothing && (ε_val < 0 || ε_val > 1)
			throw(ArgumentError("Emissivity must be in [0, 1], got $ε_val"))
		end
	end
	MarshakRadiationBC(I_external, T_ambient, emissivity)
end

# ============================================================================
# Reflective Boundary Condition (For Back/Insulated Surfaces)
# ============================================================================

"""
	ReflectiveRadiationBC <: AbstractRadiationBC

Reflective (adiabatic) boundary condition for radiation.

Imposes zero normal gradient on the incident radiation field:
	∂G/∂n = 0

This is appropriate for:
- Insulated back faces (no radiation escape)
- Symmetry planes
- Reflective substrates
- Impermeable boundaries where radiation cannot escape

# Physical Interpretation
A reflective BC means all radiation incident on the boundary is reflected
back into the domain. No net radiative flux crosses the boundary.

# Example
```julia
# Back face of sample (insulated)
radiation_bcs = RadiationBCSet(
	top = MarshakRadiationBC(I_external = 25e3, T_ambient = 300.0),
	bottom = ReflectiveRadiationBC()  # Adiabatic
)
```

# See Also
- `MarshakRadiationBC` for exposed surfaces
- `AdiabaticBC` (thermal BC) for thermal insulation
"""
struct ReflectiveRadiationBC <: AbstractRadiationBC end

# ============================================================================
# Dirichlet Boundary Condition (For Prescribed Radiation)
# ============================================================================

"""
	DirichletRadiationBC{FG} <: AbstractRadiationBC

Prescribed incident radiation intensity at a boundary.

Sets the value of G directly at the boundary:
	G_boundary = G_prescribed

This is rarely used in practice but can be useful for:
- Coupling to external radiation field calculations
- Validation against analytical solutions
- Specialized boundary treatments

# Fields  
- `G::FG`: Prescribed G value [W/m²]. Can be constant or time-dependent `G(t)`.

# Physical Meaning
For any surface at temperature T_b in volumetric equilibrium:
	G = 4σT_b⁴ (blackbody, standard P1 approximation)

Note: In P1, the volumetric equilibrium is ALWAYS blackbody (ε=1 locally).
Surface emissivity ε < 1 only appears in boundary conditions (Marshak BC)
where it represents partial reflection at the boundary, not in G itself.

# Example
```julia
# Isothermal black surface at 500 K
bc = DirichletRadiationBC(G = 4 * STEFAN_BOLTZMANN * 500.0^4)

# Time-varying prescribed G
bc = DirichletRadiationBC(G = t -> 1000.0 * (1 + 0.5*sin(2π*t/60)))
```

# Note
Using DirichletRadiationBC requires careful consideration of the physics.
In most pyrolysis applications, `MarshakRadiationBC` is preferred as it
properly couples external sources with surface emission.

# See Also
- `MarshakRadiationBC` for typical exposed surfaces
- `DirichletThermalBC` for prescribed temperature
"""
struct DirichletRadiationBC{FG} <: AbstractRadiationBC
	G::FG  # Prescribed incident radiation [W/m²]
end

# Keyword constructor
DirichletRadiationBC(; G) = DirichletRadiationBC(G)

# ============================================================================
# Default Radiation BCs Factory
# ============================================================================

"""
	default_radiation_bcs(I_external::Float64=0.0; T_ambient::Float64=300.0)

Create default radiation boundary conditions for a 1D domain.

Returns a RadiationBCSet with:
- Top (boundary 1): MarshakRadiationBC with specified external flux
- Bottom (boundary 2): ReflectiveRadiationBC (adiabatic)

# Arguments
- `I_external`: External radiation flux at top surface [W/m²]
- `T_ambient`: Ambient temperature for radiation exchange [K]
- `emissivity`: Surface emissivity (default 0.9)

# Example
```julia
# Default for cone calorimeter at 50 kW/m²
radiation_bcs = default_radiation_bcs(50e3, T_ambient=300.0)
```
"""
function default_radiation_bcs(I_external::Float64 = 0.0;
	T_ambient::Float64 = 300.0,
	emissivity::Float64 = 0.9)
	RadiationBCSet(
		top = MarshakRadiationBC(I_external = I_external, T_ambient = T_ambient, emissivity = emissivity),
		bottom = ReflectiveRadiationBC(),
	)
end

# =========================================================================
# Type-stable radiation BC container
# =========================================================================

"""
	RadiationBCSet{B1, B2}

Type-stable container for radiation boundary conditions (top, bottom).

This mirrors `BoundaryConditionSet` for thermal/mass BCs and avoids
Dict lookups inside hot radiation kernels.
"""
struct RadiationBCSet{B1 <: AbstractRadiationBC, B2 <: AbstractRadiationBC}
	top::B1
	bottom::B2
end

"""
	RadiationBCSet(; top=ReflectiveRadiationBC(), bottom=ReflectiveRadiationBC())

Keyword constructor with defaults.
"""
function RadiationBCSet(; top::AbstractRadiationBC = ReflectiveRadiationBC(),
	bottom::AbstractRadiationBC = ReflectiveRadiationBC())
	RadiationBCSet(top, bottom)
end

@inline get_radiation_bc(bc_set::RadiationBCSet, ::Val{1}) = bc_set.top
@inline get_radiation_bc(bc_set::RadiationBCSet, ::Val{2}) = bc_set.bottom

@inline function get_radiation_bc(bc_set::RadiationBCSet, boundary_id::Int)
	boundary_id == 1 && return bc_set.top
	boundary_id == 2 && return bc_set.bottom
	error("Invalid boundary_id: $boundary_id (expected 1 or 2)")
end
