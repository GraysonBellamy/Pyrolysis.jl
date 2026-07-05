"""
SimulationMode trait bundle.

The residual specializes on a single `SimulationMode{Mesh,Rad,Trans,Chi}`
type parameter rather than fanning out across orthogonal axes. Adding a
new orthogonal axis (e.g. a future combustion model) extends the bundled
type, not every signature.

**`energy_form` is a runtime field, not a trait**: it adds one additive
term to the surface energy balance and one to the cell residual; a
runtime `if` is cheaper than a separate compiled method.
"""

# ============================================================================
# Mesh motion
# ============================================================================

"""
	abstract type MeshMode

Trait axis for mesh motion. Subtypes: [`Eulerian`](@ref), [`ALE`](@ref).
"""
abstract type MeshMode end

"Eulerian mesh — fixed nodes, no mesh motion."
struct Eulerian <: MeshMode end

"ALE mesh — nodes move, geometry is part of the state vector."
struct ALE <: MeshMode end

# ============================================================================
# Radiation transport
# ============================================================================

"""
	abstract type RadiationMode

Trait axis for radiation transport. Subtypes: [`NoRadiation`](@ref),
[`SurfaceAbsorption`](@ref), [`BeerLambert`](@ref).

The P1 model lives in `Pyrolysis.Experimental` and adds its own
`RadiationMode` subtype there.
"""
abstract type RadiationMode end

"No radiation transport."
struct NoRadiation <: RadiationMode end

"Surface absorption: incident flux deposited entirely at the top face."
struct SurfaceAbsorption <: RadiationMode end

"Beer-Lambert: in-depth absorption following exp(-κz)."
struct BeerLambert <: RadiationMode end

# ============================================================================
# Mass transport
# ============================================================================

"""
	abstract type TransportMode

Trait axis for gas-phase mass transport. Subtypes: [`Fickian`](@ref),
[`DarcyFick`](@ref).
"""
abstract type TransportMode end

"Fickian-only diffusion (no pressure-driven flow)."
struct Fickian <: TransportMode end

"Darcy + Fick: pressure-driven flow plus Fickian diffusion."
struct DarcyFick <: TransportMode end

# ============================================================================
# Variable cross-section
# ============================================================================

"""
	abstract type ChiMode

Trait axis for the variable-cross-section / lateral-shrinkage extension.
Subtypes: [`NoChi`](@ref), [`WithChi`](@ref).

When `NoChi` (the default), the χ block is not allocated, packed,
unpacked, or differentiated. `WithChi` enables full machinery.
"""
abstract type ChiMode end

"Cross-section is constant; no χ state."
struct NoChi <: ChiMode end

"Variable cross-section / lateral shrinkage active; χ is a packed state field."
struct WithChi <: ChiMode end

# ============================================================================
# Bundled simulation mode
# ============================================================================

"""
	SimulationMode{Mesh, Rad, Trans, Chi}

Bundled trait carrying all four orthogonal axes as type parameters.
The residual takes `SimulationMode` as a *single* parameter so adding new
axes does not multiply method-table size.
"""
struct SimulationMode{Mesh<:MeshMode,Rad<:RadiationMode,Trans<:TransportMode,Chi<:ChiMode} end

"""
	SimulationMode(; mesh=Eulerian(), radiation=NoRadiation(),
	                 transport=Fickian(), chi=NoChi())

Construct a `SimulationMode` from trait instances. This is the **one**
way to build a mode; each axis takes a `MeshMode` / `RadiationMode` /
`TransportMode` / `ChiMode` instance respectively.

```julia
using Pyrolysis.Problem: ALE, BeerLambert, DarcyFick, NoChi
mode = SimulationMode(; mesh=ALE(), radiation=BeerLambert(),
                       transport=DarcyFick(), chi=NoChi())
```
"""
@inline SimulationMode(;
    mesh::MeshMode = Eulerian(),
    radiation::RadiationMode = NoRadiation(),
    transport::TransportMode = Fickian(),
    chi::ChiMode = NoChi(),
) = SimulationMode{typeof(mesh),typeof(radiation),typeof(transport),typeof(chi)}()

# Trait extractors --- compile-time, free. These are the canonical accessors;
# they pull instances out of the type parameters so callers can dispatch with
# `mesh_mode(m) isa ALE` etc.
@inline mesh_mode(::SimulationMode{M,R,T,C}) where {M,R,T,C} = M()
@inline radiation_mode(::SimulationMode{M,R,T,C}) where {M,R,T,C} = R()
@inline transport_mode(::SimulationMode{M,R,T,C}) where {M,R,T,C} = T()
@inline chi_mode(::SimulationMode{M,R,T,C}) where {M,R,T,C} = C()
