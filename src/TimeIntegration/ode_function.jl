"""
ODE RHS assembly for DifferentialEquations.jl integration.

Reference: Design Document Section 7.2

This module provides the ODE function wrapper that interfaces between
DifferentialEquations.jl solvers and the Pyrolysis discretization.
"""

using ..Discretization: RHSCache, compute_rhs!, resize_rhs_cache!, compute_rhs_darcy_fick!, update_darcy_properties!
using ..Discretization: update_radiation_properties!, compute_radiation_flux_divergence!, solve_quasi_steady_radiation!
using ..Core: has_darcy_flow

# ============================================================================
# PyrolysisODEFunction
# ============================================================================

"""
	PyrolysisODEFunction{M,Mat,BCS,PC,RC,SV}

ODE right-hand side function for pyrolysis problems.

This struct is callable and computes du/dt given the current state u.
It interfaces with DifferentialEquations.jl by implementing the ODE function
interface: `f(du, u, p, t)`.

# Type Parameters
- `M`: Mesh type
- `Mat`: Material type
- `BCS`: BoundaryConditionSet type (type-stable BC storage)
- `PC`: PropertyCache type
- `RC`: RHSCache type
- `SV`: StateVector type

# Fields
- `mesh::M`: The computational mesh
- `material::Mat`: Material properties and reactions
- `bc_set::BCS`: BoundaryConditionSet with all boundary conditions (type-stable)
- `prop_cache::PC`: Pre-allocated property cache
- `rhs_cache::RC`: Pre-allocated RHS workspace
- `sv::SV`: State vector descriptor
- `I_radiation::Float64`: Surface radiation intensity [W/m²]
- `dT_work::Vector{Float64}`: Work array for temperature RHS
- `dξ_work::Matrix{Float64}`: Work array for concentration RHS

# Example
```julia
bc_set = BoundaryConditionSet(
	thermal_top = HeatFluxBC(50e3) + RadiativeBC(ε=0.9, T_inf=300),
	thermal_bottom = AdiabaticBC()
)
ode_fn = PyrolysisODEFunction(mesh, material, bc_set)
prob = ODEProblem(ode_fn, u0, tspan)
sol = solve(prob, Rodas5P())
```
"""
mutable struct PyrolysisODEFunction{M, Mat, BCS, PC, RC, SV}
	mesh::M
	material::Mat
	bc_set::BCS
	prop_cache::PC
	rhs_cache::RC
	sv::SV
	I_radiation::Float64
	radiation_model::RadiationModel   # Radiation transport model
	dT_work::Vector{Float64}
	dξ_work::Matrix{Float64}

	# Pre-cached boundary info (computed once at construction)
	top_face_ids::Vector{Int}      # Faces at top boundary (boundary_id=1)
	bottom_face_ids::Vector{Int}   # Faces at bottom boundary (boundary_id=2)
	top_cell_ids::Vector{Int}      # Cells adjacent to top boundary
	bottom_cell_ids::Vector{Int}   # Cells adjacent to bottom boundary

	# Flame state for FlameBC (stateless FlameBC, thread-safe)
	flame_state::Dict{Int, Bool}   # boundary_id => flame_active
end

# Helper function to extract boundary face and cell IDs
function _get_boundary_info(mesh, boundary_id::Int)
	face_ids = Int[]
	cell_ids = Int[]
	if haskey(mesh.boundary_faces, boundary_id)
		for f in mesh.boundary_faces[boundary_id]
			push!(face_ids, f)
			face = mesh.faces[f]
			left_id, right_id = face.cell_ids
			cell_id = left_id == 0 ? right_id : left_id
			push!(cell_ids, cell_id)
		end
	end
	return face_ids, cell_ids
end

"""
	PyrolysisODEFunction(mesh, material, bc_set::BoundaryConditionSet;
						 I_radiation::Float64=0.0,
						 radiation_model::RadiationModel=NO_RADIATION)

Create a PyrolysisODEFunction using a type-stable BoundaryConditionSet.

This is the preferred constructor for maximum performance and type stability.
The BoundaryConditionSet is stored directly, enabling compile-time dispatch
in the ODE function call.

# Arguments
- `mesh`: Mesh with initialized cell_states
- `material`: Material definition
- `bc_set`: BoundaryConditionSet with all boundary conditions
- `I_radiation`: Optional surface radiation intensity [W/m²]
- `radiation_model`: Radiation transport model (NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT, P1_QUASI_STEADY)

# Example
```julia
bc_set = BoundaryConditionSet(
	thermal_top = HeatFluxBC(50e3) + RadiativeBC(ε=0.9, T_inf=300),
	thermal_bottom = AdiabaticBC(),
	mass_top = DiffusiveMassBC(Inf),
	mass_bottom = ImpermeableBC()
)
ode_fn = PyrolysisODEFunction(mesh, material, bc_set;
	radiation_model = SURFACE_ABSORPTION)
```
"""
function PyrolysisODEFunction(mesh, material, bc_set::BoundaryConditionSet;
	I_radiation::Float64 = 0.0,
	radiation_model::RadiationModel = NO_RADIATION)
	n_cells_val = length(mesh.cells)
	NC = n_components(material)
	NR = n_reactions(material)
	n_faces_val = length(mesh.faces)

	# Create caches
	prop_cache = Properties.PropertyCache(n_cells_val, n_faces_val)
	rhs_cache = RHSCache(n_cells_val, n_faces_val, NC, NR)
	sv = StateVector(n_cells_val, NC)

	# Create work arrays
	dT_work = zeros(Float64, n_cells_val)
	dξ_work = zeros(Float64, NC, n_cells_val)

	# Pre-cache boundary face and cell IDs (computed once)
	top_face_ids, top_cell_ids = _get_boundary_info(mesh, 1)
	bottom_face_ids, bottom_cell_ids = _get_boundary_info(mesh, 2)

	# Initialize flame state (for FlameBC boundaries)
	flame_state = Dict{Int, Bool}()

	PyrolysisODEFunction(
		mesh, material, bc_set,
		prop_cache, rhs_cache, sv, I_radiation, radiation_model,
		dT_work, dξ_work,
		top_face_ids, bottom_face_ids, top_cell_ids, bottom_cell_ids,
		flame_state,
	)
end

# ============================================================================
# ODE Function Interface (in-place form)
# ============================================================================

"""
	(f::PyrolysisODEFunction)(du, u, p, t)

Compute the ODE right-hand side du/dt.

This is the main entry point for DifferentialEquations.jl. It:
1. Unpacks the state vector `u` into mesh cell states
2. Updates material properties at current state
3. Computes the RHS (heat flux divergence, reactions, gas transport)
4. Packs the result into `du`

# Arguments
- `du`: Output derivative vector (same layout as u)
- `u`: Current state vector
- `p`: Parameters (unused, for interface compatibility)
- `t`: Current time
"""
function (f::PyrolysisODEFunction)(du, u, p, t)
	mesh = f.mesh
	material = f.material
	sv = f.sv
	NC = n_components(material)
	n = sv.n_cells

	# Step 1: Unpack state vector into mesh
	unpack_state!(mesh, u, sv)

	# Step 2: Update material properties
	Properties.update_properties!(f.prop_cache, mesh, material)

	# Step 3: Compute RHS (dispatch based on BC type)
	_compute_rhs_dispatch!(f, t)

	# Step 4: Pack results into du
	@inbounds for i in 1:n
		du[i] = f.dT_work[i]
		for j in 1:NC
			du[n+(j-1)*n+i] = f.dξ_work[j, i]
		end
	end

	return nothing
end

# Dispatch helper for compute_rhs! using BoundaryConditionSet (type-stable)
function _compute_rhs_dispatch!(f::PyrolysisODEFunction{M, Mat, BCS}, t) where {M, Mat, BCS <: BoundaryConditionSet}
	compute_rhs!(f.dT_work, f.dξ_work,
		f.mesh, f.prop_cache, f.rhs_cache,
		f.material, f.bc_set, t,
		f.top_face_ids, f.bottom_face_ids;
		I_radiation = f.I_radiation,
		radiation_model = f.radiation_model)
end

# ============================================================================
# Out-of-place form (for testing/compatibility)
# ============================================================================

"""
	(f::PyrolysisODEFunction)(u, p, t)

Out-of-place ODE function: returns du/dt as a new vector.

Less efficient than in-place form but useful for testing.
"""
function (f::PyrolysisODEFunction)(u, p, t)
	du = similar(u)
	f(du, u, p, t)
	return du
end

# ============================================================================
# Resize for Adaptive Mesh
# ============================================================================

"""
	resize_ode_function!(f::PyrolysisODEFunction)

Resize caches after mesh refinement/coarsening.

Call this after any mesh topology change.
"""
function resize_ode_function!(f::PyrolysisODEFunction)
	n_cells_val = length(f.mesh.cells)
	n_faces_val = length(f.mesh.faces)
	NC = n_components(f.material)
	NR = n_reactions(f.material)

	# Resize property cache
	Properties.resize_cache!(f.prop_cache, n_cells_val, n_faces_val)

	# Resize RHS cache
	resize_rhs_cache!(f.rhs_cache, n_cells_val, n_faces_val, NC, NR)

	# Update state vector descriptor
	f.sv = StateVector(n_cells_val, NC)

	# Resize work arrays
	resize!(f.dT_work, n_cells_val)
	f.dξ_work = zeros(Float64, NC, n_cells_val)

	# Refresh pre-cached boundary info after mesh topology change
	top_face_ids, top_cell_ids = _get_boundary_info(f.mesh, 1)
	bottom_face_ids, bottom_cell_ids = _get_boundary_info(f.mesh, 2)
	f.top_face_ids = top_face_ids
	f.bottom_face_ids = bottom_face_ids
	f.top_cell_ids = top_cell_ids
	f.bottom_cell_ids = bottom_cell_ids

	return f
end

# ============================================================================
# Accessors
# ============================================================================

"""
	get_mesh(f::PyrolysisODEFunction)

Get the mesh from the ODE function.
"""
get_mesh(f::PyrolysisODEFunction) = f.mesh

"""
	get_material(f::PyrolysisODEFunction)

Get the material from the ODE function.
"""
get_material(f::PyrolysisODEFunction) = f.material

"""
	get_state_vector_descriptor(f::PyrolysisODEFunction)

Get the StateVector descriptor.
"""
get_state_vector_descriptor(f::PyrolysisODEFunction) = f.sv

"""
	set_radiation!(f::PyrolysisODEFunction, I::Float64)

Set the surface radiation intensity.
"""
function set_radiation!(f::PyrolysisODEFunction, I::Float64)
	f.I_radiation = I
	return f
end

# ============================================================================
# Time-Varying Radiation Wrapper
# ============================================================================

"""
	PyrolysisODEFunctionWithRadiation{F,ODE}

Wrapper that updates radiation intensity from a time-dependent function.

# Fields
- `ode_fn::ODE`: The underlying PyrolysisODEFunction
- `I_func::F`: Function I(t) -> radiation intensity [W/m²]

# Example
```julia
ode_fn = PyrolysisODEFunction(mesh, material, bcs)
I_func = t -> 50e3 * (1 - exp(-t/10))
ode_fn_rad = PyrolysisODEFunctionWithRadiation(ode_fn, I_func)
```
"""
struct PyrolysisODEFunctionWithRadiation{F, ODE}
	ode_fn::ODE
	I_func::F
end

# Auto-generated default constructor will handle type inference

function (f::PyrolysisODEFunctionWithRadiation)(du, u, p, t)
	# Update radiation intensity
	f.ode_fn.I_radiation = f.I_func(t)
	# Call underlying function
	f.ode_fn(du, u, p, t)
end

function (f::PyrolysisODEFunctionWithRadiation)(u, p, t)
	f.ode_fn.I_radiation = f.I_func(t)
	f.ode_fn(u, p, t)
end

# ============================================================================
# PyrolysisODEFunctionALE - Full ALE Integration
# ============================================================================

"""
	PyrolysisODEFunctionALE{M,Mat,TBC,MBC,PC,RC,SV,GC}

ODE function with integrated ALE mesh motion.

The state vector includes node positions, and the RHS includes:
- dT/dt: Temperature evolution with ALE convection
- dξ/dt: Species evolution with ALE convection
- dz/dt: Node velocities (mesh motion)

This is the mathematically correct ALE formulation where mesh motion is
integrated directly into the ODE rather than applied via callbacks.

# Type Parameters
- `M`: Mesh type
- `Mat`: Material type
- `TBC`: Thermal BC type
- `MBC`: Mass BC type
- `PC`: PropertyCache type
- `RC`: RHSCache type
- `SV`: StateVectorALE type
- `GC`: GeometryCache type
- `MVM`: Mesh velocity method type

# State Vector Layout
```
u = [T₁...Tₙ, ξ₁₁...ξ_NC,n, z₁...z_{n+1}]
du = [dT₁/dt...dTₙ/dt, dξ₁₁/dt...dξ_NC,n/dt, w₁...w_{n+1}]
```

Where wᵢ = dzᵢ/dt is the mesh velocity at node i.

# Example
```julia
ode_fn = PyrolysisODEFunctionALE(mesh, material, thermal_bcs, mass_bcs)
u0 = create_state_vector_ale(mesh)
prob = ODEProblem(ode_fn, u0, tspan)
sol = solve(prob, Rodas5P())
```
"""
mutable struct PyrolysisODEFunctionALE{M, Mat, BCS, PC, RC, SV, GC}
	mesh::M
	material::Mat
	bc_set::BCS                       # BoundaryConditionSet (type-stable)
	prop_cache::PC
	rhs_cache::RC
	sv::SV
	geom_cache::GC
	I_radiation::Float64
	radiation_model::RadiationModel   # Radiation transport model

	# Work arrays for physics RHS
	dT_work::Vector{Float64}
	dξ_work::Matrix{Float64}

	# Work arrays for ALE
	w_work::Vector{Float64}           # Mesh velocities at nodes
	w_face::Vector{Float64}           # Mesh velocities at faces (for reference)
	v_gas_face::Vector{Float64}       # Gas velocity at faces (unused in upwind formulation)
	θ_work::Vector{Float64}           # Volumetric strain rate per cell [1/s] for ALE dilation

	# ALE advection arrays (CORRECTED FORMULATION)
	# Using upwind advection (w·∇φ) instead of conservative flux (∇·(φw))
	# Reference: ALE_Energy_Reformulation.md
	ale_advection::Vector{Float64}    # ALE energy advection term [K/s]
	ale_species_advection::Matrix{Float64}  # ALE species advection [NC × n_cells, kg/(m³·s)]

	# Pre-allocated work arrays for reaction evaluation (PERFORMANCE: avoid allocations in hot loops)
	rates_work::Vector{Float64}       # Reaction rates [NR]
	heats_work::Vector{Float64}       # Heats of reaction [NR]
	dξ_dt_work::Vector{Float64}       # Species source terms [NC]

	# Convective heat flux work array
	q_conv_work::Vector{Float64}      # Convective heat flux work array [n_faces]

	# Pre-cached boundary info (computed once at construction)
	top_face_ids::Vector{Int}         # Faces at top boundary (boundary_id=1)
	bottom_face_ids::Vector{Int}      # Faces at bottom boundary (boundary_id=2)
	top_cell_ids::Vector{Int}         # Cells adjacent to top boundary
	bottom_cell_ids::Vector{Int}      # Cells adjacent to bottom boundary

	# Settings
	min_thickness::Float64
	min_density::Float64              # Minimum density for regression calc
end

"""
	PyrolysisODEFunctionALE(mesh, material, bc_set::BoundaryConditionSet;
						   I_radiation=0.0,
						   radiation_model=NO_RADIATION,
						   min_thickness=1e-6,
						   min_density=1.0)

Create ALE ODE function with integrated mesh motion.

This is the preferred constructor using type-stable BoundaryConditionSet
for maximum performance and zero allocations in the hot path.

# Arguments
- `mesh`: Mesh with initialized cell_states
- `material`: Material definition (includes permeability_model for Darcy flow)
- `bc_set`: BoundaryConditionSet with all boundary conditions (type-stable)
- `I_radiation`: Optional surface radiation intensity [W/m²]
- `radiation_model`: Radiation transport model (NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT, P1_QUASI_STEADY)
- `min_thickness`: Minimum domain thickness before termination [m]
- `min_density`: Minimum surface density for regression calculation [kg/m³]

# Gas Transport
Gas transport model is determined by `material.permeability_model`:
- `nothing` (default): Pure Fickian diffusion
- `ConstantPermeability(κ)`: Darcy-Fick with constant permeability
- `KozenyCarmanPermeability(κ₀, φ₀)`: Darcy-Fick with porosity-dependent permeability

# Returns
A callable ALE ODE function for use with DifferentialEquations.jl.

# Example
```julia
bc_set = BoundaryConditionSet(
	thermal_top = HeatFluxBC(50e3) + RadiativeBC(ε=0.9, T_inf=300),
	thermal_bottom = AdiabaticBC(),
	mass_top = DiffusiveMassBC(Inf),
	mass_bottom = ImpermeableBC()
)
ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set;
	radiation_model = BEER_LAMBERT)
```
"""
function PyrolysisODEFunctionALE(mesh::Mesh1D{NC}, material::Material{NC},
	bc_set::BoundaryConditionSet;
	I_radiation::Float64 = 0.0,
	radiation_model::RadiationModel = NO_RADIATION,
	min_thickness::Float64 = 1e-6,
	min_density::Float64 = 1.0) where NC
	n_cells = length(mesh.cells)
	n_faces = length(mesh.faces)
	n_nodes = length(mesh.nodes)
	NR = n_reactions(material)

	# Create standard caches
	prop_cache = Properties.PropertyCache(n_cells, n_faces)
	rhs_cache = RHSCache(n_cells, n_faces, NC, NR)
	sv = StateVectorALE(n_cells, Val(NC))

	# Create geometry cache for allocation-free geometry computation
	geom_cache = Geometry.GeometryCache(n_cells, n_faces)

	# Work arrays for physics RHS
	dT_work = zeros(n_cells)
	dξ_work = zeros(NC, n_cells)

	# Work arrays for ALE
	w_work = zeros(n_nodes)
	w_face = zeros(n_faces)
	v_gas_face = zeros(n_faces)

	# Work array for volumetric strain rate (ALE dilation)
	θ_work = zeros(n_cells)

	# ALE advection arrays (CORRECTED: cell-centered, not face-based)
	ale_advection = zeros(n_cells)
	ale_species_advection = zeros(NC, n_cells)

	# Pre-allocated work arrays for reaction evaluation (PERFORMANCE: avoid allocations)
	rates_work = zeros(NR)
	heats_work = zeros(NR)
	dξ_dt_work = zeros(NC)

	# Work array for convective heat flux (Darcy)
	q_conv_work = zeros(n_faces)

	# Pre-cache boundary face and cell IDs (computed once at construction)
	top_face_ids, top_cell_ids = _get_boundary_info(mesh, 1)
	bottom_face_ids, bottom_cell_ids = _get_boundary_info(mesh, 2)

	PyrolysisODEFunctionALE(
		mesh, material, bc_set,
		prop_cache, rhs_cache, sv, geom_cache,
		I_radiation, radiation_model,
		dT_work, dξ_work,
		w_work, w_face, v_gas_face, θ_work,
		ale_advection, ale_species_advection,
		rates_work, heats_work, dξ_dt_work,
		q_conv_work,
		top_face_ids, bottom_face_ids, top_cell_ids, bottom_cell_ids,
		min_thickness, min_density,
	)
end

# ============================================================================
# Main ALE ODE Function (In-Place)
# ============================================================================

"""
	(f::PyrolysisODEFunctionALE)(du, u, p, t)

Compute full ALE RHS including mesh motion.

This is the main entry point for DifferentialEquations.jl. It:
1. Unpacks the state vector (T, ξ, z) and updates mesh geometry
2. Updates material properties
3. Computes mesh velocity field
4. Computes ALE convective fluxes
5. Computes standard physics RHS
6. Packs result into du = [dT/dt, dξ/dt, dz/dt]

# Arguments
- `du`: Output derivative vector
- `u`: Current state vector [T, ξ, z]
- `p`: Parameters (unused, for interface compatibility)
- `t`: Current time
"""
function (f::PyrolysisODEFunctionALE{M, Mat, BCS, PC, RC, SV, GC})(
	du, u, p, t) where {M, Mat, BCS <: BoundaryConditionSet, PC, RC, SV <: StateVectorALE{NC}, GC} where NC

	mesh = f.mesh
	material = f.material
	sv = f.sv
	n = sv.n_cells
	n_nodes = sv.n_nodes

	# =========================================================================
	# Step 1: Unpack state and update geometry
	# =========================================================================
	unpack_state_ale!(mesh, u, sv)

	# =========================================================================
	# Step 2: Update material properties at current state
	# =========================================================================
	Properties.update_properties!(f.prop_cache, mesh, material)

	# =========================================================================
	# Step 3: Compute standard physics RHS (conduction, reactions, transport)
	# =========================================================================
	# PERFORMANCE: Compute physics RHS FIRST so we can reuse the reaction-derived
	# species source terms (dξ_rxn) in the mesh velocity calculation, avoiding
	# redundant reaction evaluation (~15% speedup for complex kinetics)

	# Choose between standard RHS and Darcy-Fick RHS based on component permeability values
	if has_darcy_flow(material)
		# Use Darcy-Fick formulation with pressure-driven gas transport
		# Uses component-based permeability via mixing rules
		compute_rhs_darcy_fick!(f.dT_work, f.dξ_work, mesh, f.prop_cache, f.rhs_cache,
			material, f.bc_set, t, f.top_face_ids, f.bottom_face_ids;
			I_radiation = f.I_radiation,
			radiation_model = f.radiation_model)
	else
		# Use standard diffusion-only RHS
		compute_rhs!(f.dT_work, f.dξ_work, mesh, f.prop_cache, f.rhs_cache,
			material, f.bc_set, t, f.top_face_ids, f.bottom_face_ids;
			I_radiation = f.I_radiation,
			radiation_model = f.radiation_model)
	end

	# =========================================================================
	# Step 4: Compute mesh velocity field using PRECOMPUTED species sources
	#         Also compute θ (volumetric strain rate) for ALE dilation correction
	#         PERFORMANCE: Pass precomputed dξ_rxn to avoid re-evaluating reactions
	#         PERFORMANCE: Pass pre-allocated work arrays to avoid allocations
	# =========================================================================
	_compute_ale_mesh_velocity!(f.w_work, mesh, f.prop_cache, material,
		t, f.min_density;
		θ_out = f.θ_work,
		precomputed_dξ_rxn = f.rhs_cache.dξ_rxn,
		rates = f.rates_work,
		heats = f.heats_work,
		dξ_dt = f.dξ_dt_work)

	# Interpolate node velocities to faces
	_interpolate_node_to_face!(f.w_face, f.w_work, mesh)

	# =========================================================================
	# Step 5: Compute gas velocity at faces (simplified model)
	# =========================================================================
	# In diffusion-dominated regime, bulk gas velocity ≈ 0
	# (kept for potential future use, not used in upwind advection)
	fill!(f.v_gas_face, 0.0)

	# =========================================================================
	# Step 6: Compute ALE advection terms (CORRECTED FORMULATION)
	# =========================================================================
	# Using upwind advection (w·∇φ) instead of conservative flux divergence (∇·(φw))
	# This avoids spurious energy sources/sinks at fixed boundaries
	# Reference: ALE_Energy_Reformulation.md
	Discretization.compute_ale_advection_term!(f.ale_advection, mesh, f.w_work)
	Discretization.compute_ale_species_advection!(f.ale_species_advection, mesh, material, f.w_work)

	# =========================================================================
	# Step 7: Subtract ALE advection contributions
	# =========================================================================
	# Energy: dT/dt = (physics terms) - w·∇T
	# Species: dξ/dt = (physics terms) - w·∇ξ
	_subtract_ale_advection!(f.dT_work, f.dξ_work,
		f.ale_advection, f.ale_species_advection, NC)

	# =========================================================================
	# Step 7.5: Apply ALE dilation correction for solid species
	# =========================================================================
	# For solids in ALE with volume change, mass conservation requires:
	#   dξ/dt = S - ξ × θ   (where θ is volumetric strain rate)
	# This corrects the species equation to account for concentration changes
	# due to cell volume changes (shrinkage/expansion).
	_apply_ale_dilation_correction!(f.dξ_work, f.θ_work, mesh, material)

	# NOTE on energy dilation: Unlike species, the energy equation does NOT need
	# an explicit -T×θ correction because the temperature formulation implicitly
	# handles this through the time-varying ρcₚ denominator. The terms cancel:
	#   Full equation: dT/dt = sources/ρcₚ - T×θ - (T/ρcₚ)(dρcₚ/dt)
	#   When ρcₚ changes due to reactions, the last two terms largely cancel.
	# See Enthalpy_Formulation_Refactoring_Plan.md for the proper enthalpy approach.

	# =========================================================================
	# Step 8: Pack results into du
	# =========================================================================

	# Temperature derivatives
	@inbounds for i in 1:n
		du[i] = f.dT_work[i]
	end

	# Concentration derivatives
	@inbounds for i in 1:n
		for j in 1:NC
			du[sv.ξ_offset+(j-1)*n+i] = f.dξ_work[j, i]
		end
	end

	# Node velocities (dz/dt = w)
	@inbounds for i in 1:n_nodes
		du[sv.z_offset+i] = f.w_work[i]
	end

	return nothing
end

# ============================================================================
# Out-of-place form
# ============================================================================

"""
	(f::PyrolysisODEFunctionALE)(u, p, t)

Out-of-place ALE ODE function: returns du/dt as a new vector.
"""
function (f::PyrolysisODEFunctionALE)(u, p, t)
	du = similar(u)
	f(du, u, p, t)
	return du
end

# ============================================================================
# ALE Helper Functions
# ============================================================================

"""
	_compute_ale_mesh_velocity!(w, mesh, prop_cache, material, mass_bcs, t, method, min_density)

Compute mesh velocity field for ALE using CUMULATIVE VOLUME CHANGE.

# Physics-Based Mesh Velocity
For pyrolysis with volume change, the mesh velocity at position z represents
the CUMULATIVE volume change from the fixed boundary (z=0) to z:

	w(z) = ∫₀ᶻ θ(z') dz'

Where θ(z) is the LOCAL volume change rate:
	θ = Σⱼ (1/ρⱼ) × (dξⱼ/dt)   for solid components

For discrete cells (nodes numbered 1 to n_nodes, cells 1 to n_cells):
	w[1] = 0           (bottom node fixed)
	w[i+1] = w[i] + θᵢ × Δzᵢ   for each cell i

This ensures:
1. Each cell shrinks according to ITS OWN reaction rate
2. Surface velocity = integral of all volume changes
3. Cells with high reaction rates contribute more to shrinkage
4. Total thickness correctly tracks total solid mass

# Why NOT Laplacian Smoothing
The previous implementation used Laplacian smoothing (linear interpolation)
for interior nodes. This is WRONG for pyrolysis because:
1. Interior cells DO change volume (reactions happening everywhere)
2. Linear interpolation ignores WHERE reactions are happening
3. When surface cell is depleted but interior is actively reacting,
   Laplacian smoothing grossly underestimates surface regression

Example: If surface cell is depleted (ξ≈10) but cell below has ξ≈800:
- Laplacian: w_surface ≈ θ_surface × Δz ≈ 0 (surface cell nearly inert)
- Correct: w_surface = Σ(θᵢ × Δzᵢ) includes the active interior cell
"""
function _compute_ale_mesh_velocity!(w::Vector{Float64}, mesh::Mesh1D{NC},
	prop_cache::PropertyCache,
	material::Material{NC},
	t::Real,
	min_density::Float64;
	θ_out::Union{Nothing, Vector{Float64}} = nothing,
	precomputed_dξ_rxn::Union{Nothing, Matrix{Float64}} = nothing,
	rates::Union{Nothing, Vector{Float64}} = nothing,
	heats::Union{Nothing, Vector{Float64}} = nothing,
	dξ_dt::Union{Nothing, Vector{Float64}} = nothing) where NC
	n_cells = length(mesh.cells)
	n_nodes = length(mesh.nodes)
	NR = n_reactions(material)

	# PERFORMANCE OPTIMIZATION: Use precomputed species source terms if available
	# This avoids redundant reaction evaluation (~15% speedup for complex kinetics)
	use_precomputed = precomputed_dξ_rxn !== nothing

	# Only allocate reaction arrays if we need to compute reactions
	if !use_precomputed
		if rates === nothing
			rates = zeros(NR)
		end
		if heats === nothing
			heats = zeros(NR)
		end
	end
	if dξ_dt === nothing
		dξ_dt = zeros(NC)
	end

	# =========================================================================
	# Compute CUMULATIVE mesh velocity from volume change of ALL cells
	# =========================================================================
	# 
	# The mesh velocity at node i represents how much the material BELOW that
	# node has shrunk (or expanded). This is a cumulative sum:
	#   w[i] = Σⱼ₌₁ⁱ⁻¹ θⱼ × Δzⱼ
	#
	# Where θⱼ is the volume change rate of cell j.

	# Bottom node is fixed (node 1, at z=0)
	w[1] = 0.0

	# Cumulative velocity (starts at 0 at bottom)
	cumulative_velocity = 0.0

	# Loop through all cells from bottom to top
	@inbounds for cell_idx in 1:n_cells
		cell_state = mesh.cell_states[cell_idx]
		T_cell = cell_state.T

		# Cell thickness
		Δz_cell = mesh.cell_volumes[cell_idx] / mesh.cross_section_area

		# Get species source terms (dξ/dt from reactions)
		# PERFORMANCE: Use precomputed values from RHS if available
		if use_precomputed
			# Copy from precomputed matrix (already evaluated in compute_rhs!)
			for j in 1:NC
				dξ_dt[j] = precomputed_dξ_rxn[j, cell_idx]
			end
		else
			# Evaluate reaction rates for this cell (fallback path)
			ξ_cell = cell_state.ξ
			fill!(rates, 0.0)
			fill!(heats, 0.0)
			Physics.evaluate_reactions!(rates, heats, ξ_cell, T_cell, material)
			Physics.species_source_terms!(dξ_dt, rates, material, NC)
		end

		# Compute volume change rate θ for this cell using dispatch-based helper
		# PERFORMANCE: Uses Physics.volume_change_rate which avoids heterogeneous
		# tuple allocations via dispatch-based unrolling
		θ_cell = Physics.volume_change_rate(dξ_dt, T_cell, material)

		# Store θ if output array provided (for ALE dilation correction)
		if θ_out !== nothing
			θ_out[cell_idx] = θ_cell
		end

		# Add this cell's contribution to cumulative velocity
		# θ < 0 when decomposing → cumulative_velocity decreases (regression)
		cumulative_velocity += θ_cell * Δz_cell

		# The node to the RIGHT of this cell gets the cumulative velocity
		# In our convention: cell i has nodes (i, i+1)
		right_node = mesh.cells[cell_idx].node_ids[2]
		w[right_node] = cumulative_velocity
	end

	return nothing
end

"""
	_interpolate_node_to_face!(w_face, w_nodes, mesh)

Interpolate node velocities to face positions.

For 1D, each face is at a node location, so this is a direct mapping.
"""
function _interpolate_node_to_face!(w_face::Vector{Float64},
	w_nodes::Vector{Float64},
	mesh::Mesh1D)
	@inbounds for f in eachindex(mesh.faces)
		face = mesh.faces[f]
		w_face[f] = w_nodes[face.node_id]
	end
	return nothing
end

"""
	_subtract_ale_advection!(dT, dξ, ale_advection, ale_species_advection, NC)

Subtract ALE advection contributions from RHS.

This is the CORRECT formulation using upwind advection (w·∇φ) instead of
conservative flux divergence (∇·(φw)).

The ALE equations are:
	For temperature (advected with mesh): dT/dt = (physics RHS) - w·∇T
	For gas species (v≈0, mesh moves away): dξ/dt = (physics RHS) + w·∇ξ
	For solid species: dξ/dt = (physics RHS) - ξ·θ  (handled separately via dilation correction)

# Arguments
- `dT`: Temperature RHS vector [K/s], modified in place
- `dξ`: Species RHS matrix [NC × n_cells, kg/(m³·s)], modified in place
- `ale_advection`: ALE energy advection term [K/s] for each cell
- `ale_species_advection`: ALE species advection [NC × n_cells, kg/(m³·s)]
- `NC`: Number of components

# Reference
ALE_Energy_Reformulation.md
"""
function _subtract_ale_advection!(dT::Vector{Float64}, dξ::Matrix{Float64},
	ale_advection::Vector{Float64},
	ale_species_advection::Matrix{Float64},
	NC::Int)
	n_cells = length(dT)

	@inbounds for i in 1:n_cells
		# Energy: dT/dt -= w·∇T  (mesh moves, so temperature advects in opposite direction)
		dT[i] -= ale_advection[i]

		# Species advection for GAS components:
		# The ALE equation is: ∂ξ/∂t|_χ + (v - w)·∇ξ = S
		# For gas with v ≈ 0: ∂ξ/∂t|_χ = S - (-w)·∇ξ = S + w·∇ξ
		# Since ale_species_advection stores w·∇ξ, we ADD it.
		# Note: For SOLID components, ale_species_advection[j,i] = 0 (set in compute function)
		for j in 1:NC
			dξ[j, i] += ale_species_advection[j, i]
		end
	end

	return nothing
end

"""
	_apply_ale_dilation_correction!(dξ, θ, mesh, material, NC)

Apply the ALE dilation correction for solid species.

In the ALE formulation with volume-changing cells, mass conservation requires:
	d(ξV)/dt = S × V    (mass changes only due to reactions)

Expanding: V × dξ/dt + ξ × dV/dt = S × V
Since dV/dt = V × θ (where θ is volumetric strain rate):
	dξ/dt = S - ξ × θ

This function applies the correction term -ξ × θ to solid species.
Without this correction, mass appears to deplete faster than volume shrinks.

# Arguments
- `dξ`: Species RHS matrix [NC × n_cells], modified in place
- `θ`: Volumetric strain rate per cell [1/s]
- `mesh`: Mesh with current cell states (provides ξ values)
- `material`: Material (provides component phases)
- `NC`: Number of components

# Note
This correction is ONLY applied to SOLID components because:
- Solids move with the mesh (v = w), so no advection term
- But they still experience dilution/concentration from volume change
- Gases have advection terms that already account for flux divergence
"""
function _apply_ale_dilation_correction!(dξ::Matrix{Float64},
	θ::Vector{Float64},
	mesh::Mesh1D{NC},
	material::Material{NC}) where NC
	n_cells = length(mesh.cells)

	@inbounds for i in 1:n_cells
		θ_i = θ[i]
		ξ_i = mesh.cell_states[i].ξ

		# Use unrolled version for type stability
		_apply_dilation_correction_unrolled!(dξ, i, ξ_i, θ_i, material)
	end

	return nothing
end

# Dispatch-based unrolled dilation correction for type stability
@generated function _apply_dilation_correction_unrolled!(dξ::Matrix{Float64}, i::Int,
	ξ_i::NTuple{NC}, θ_i::Float64, material::Material{NC}) where NC
	exprs = Expr[]
	for j in 1:NC
		push!(exprs, quote
			@inbounds if material.components[$j].phase == SOLID
				# Correction: dξ/dt -= ξ × θ
				# This ensures mass conservation: d(ξV)/dt = S × V
				dξ[$j, i] -= ξ_i[$j] * θ_i
			end
		end)
	end
	return Expr(:block, exprs...)
end

# ============================================================================
# ALE ODE Function Accessors and Utilities
# ============================================================================

"""
	get_mesh(f::PyrolysisODEFunctionALE)

Get the mesh from the ALE ODE function.
"""
get_mesh(f::PyrolysisODEFunctionALE) = f.mesh

"""
	get_material(f::PyrolysisODEFunctionALE)

Get the material from the ALE ODE function.
"""
get_material(f::PyrolysisODEFunctionALE) = f.material

"""
	get_state_vector_descriptor(f::PyrolysisODEFunctionALE)

Get the StateVectorALE descriptor.
"""
get_state_vector_descriptor(f::PyrolysisODEFunctionALE) = f.sv

"""
	set_radiation!(f::PyrolysisODEFunctionALE, I::Float64)

Set the surface radiation intensity.
"""
function set_radiation!(f::PyrolysisODEFunctionALE, I::Float64)
	f.I_radiation = I
	return f
end

"""
	get_mesh_velocity(f::PyrolysisODEFunctionALE)

Get the current mesh velocity field (from last RHS evaluation).
"""
get_mesh_velocity(f::PyrolysisODEFunctionALE) = copy(f.w_work)

# ============================================================================
# Resize for Adaptive Mesh (ALE)
# ============================================================================

"""
	resize_ode_function_ale!(f::PyrolysisODEFunctionALE, n_cells_new::Int, n_faces_new::Int)

Resize all caches and work arrays in an ALE ODE function after mesh topology change.

This function is called during adaptive mesh operations (h-refinement, h-coarsening,
surface cell depletion) when the number of cells/faces changes. It resizes:
- Property cache
- RHS cache
- Geometry cache
- State vector descriptor
- Physics work arrays (dT_work, dξ_work)
- ALE work arrays (w_work, w_face, v_gas_face)
- ALE flux arrays (F_ale_energy, J_ale_species)

# Arguments
- `f::PyrolysisODEFunctionALE`: The ODE function to resize
- `n_cells_new::Int`: New number of cells
- `n_faces_new::Int`: New number of faces (typically n_cells + 1 for 1D)

# Notes
- The mesh itself should be modified BEFORE calling this function
- After calling this, the state vector must also be resized separately
- All array contents are invalidated; recompute geometry/properties as needed

# Example
```julia
# After merging surface cell:
merge_surface_cell!(mesh)
n_cells_new = length(mesh.cells)
n_faces_new = length(mesh.faces)
resize_ode_function_ale!(ode_fn, n_cells_new, n_faces_new)
```
"""
function resize_ode_function_ale!(f::PyrolysisODEFunctionALE, n_cells_new::Int, n_faces_new::Int)
	n_nodes_new = n_faces_new  # In 1D, n_faces = n_nodes
	NC = n_components(f.material)
	NR = n_reactions(f.material)

	# Resize property cache
	Properties.resize_cache!(f.prop_cache, n_cells_new, n_faces_new)

	# Resize RHS cache
	resize_rhs_cache!(f.rhs_cache, n_cells_new, n_faces_new, NC, NR)

	# Resize geometry cache
	Geometry.resize_geometry_cache!(f.geom_cache, n_cells_new, n_faces_new)

	# Update state vector descriptor with new dimensions
	f.sv = StateVectorALE(n_cells_new, Val(NC))

	# Resize physics work arrays
	resize!(f.dT_work, n_cells_new)
	f.dξ_work = zeros(Float64, NC, n_cells_new)

	# Resize ALE work arrays
	resize!(f.w_work, n_nodes_new)
	resize!(f.w_face, n_faces_new)
	resize!(f.v_gas_face, n_faces_new)
	resize!(f.θ_work, n_cells_new)

	# Resize ALE advection arrays (CORRECTED: cell-centered)
	resize!(f.ale_advection, n_cells_new)
	f.ale_species_advection = zeros(Float64, NC, n_cells_new)

	# Refresh pre-cached boundary info after mesh topology change
	top_face_ids, top_cell_ids = _get_boundary_info(f.mesh, 1)
	bottom_face_ids, bottom_cell_ids = _get_boundary_info(f.mesh, 2)
	f.top_face_ids = top_face_ids
	f.bottom_face_ids = bottom_face_ids
	f.top_cell_ids = top_cell_ids
	f.bottom_cell_ids = bottom_cell_ids

	return f
end

# ============================================================================
# Time-Varying Radiation Wrapper for ALE
# ============================================================================

"""
	PyrolysisODEFunctionALEWithRadiation{F,ODE}

Wrapper that updates radiation intensity from a time-dependent function.

# Example
```julia
ode_fn = PyrolysisODEFunctionALE(mesh, material, bcs)
I_func = t -> 50e3 * (1 - exp(-t/10))
ode_fn_rad = PyrolysisODEFunctionALEWithRadiation(ode_fn, I_func)
```
"""
struct PyrolysisODEFunctionALEWithRadiation{F, ODE}
	ode_fn::ODE
	I_func::F
end

function (f::PyrolysisODEFunctionALEWithRadiation)(du, u, p, t)
	f.ode_fn.I_radiation = f.I_func(t)
	f.ode_fn(du, u, p, t)
end

function (f::PyrolysisODEFunctionALEWithRadiation)(u, p, t)
	f.ode_fn.I_radiation = f.I_func(t)
	f.ode_fn(u, p, t)
end

"""
	resize_ode_function_ale!(f::PyrolysisODEFunctionALEWithRadiation, n_cells_new::Int, n_faces_new::Int)

Resize the wrapped ALE ODE function inside a radiation wrapper.
"""
function resize_ode_function_ale!(f::PyrolysisODEFunctionALEWithRadiation, n_cells_new::Int, n_faces_new::Int)
	resize_ode_function_ale!(f.ode_fn, n_cells_new, n_faces_new)
	return f
end

# ============================================================================
# P1 Volumetric Radiation ODE Function
# ============================================================================

"""
	PyrolysisODEFunctionALEWithP1Radiation{ODE, RC, RBC}

ALE ODE function with P1 volumetric radiation model.

This wrapper extends PyrolysisODEFunctionALE to include P1 radiative transport:
1. Maintains radiation cache with G field and work arrays
2. At each timestep, solves quasi-steady P1 equation for G
3. Adds radiative flux divergence to energy equation

# Fields
- `ode_fn::ODE`: Wrapped ALE ODE function
- `rad_cache::RC`: Radiation cache with G field and work arrays
- `radiation_bcs::RBC`: Radiation boundary conditions (Dict or RadiationBCSet)
- `use_quasi_steady::Bool`: If true, solve quasi-steady P1 at each timestep
- `G_work::Vector{Float64}`: Work array for radiation intensity RHS

# Example
```julia
# Create base ALE ODE function
ode_fn = PyrolysisODEFunctionALE(mesh, material, thermal_bcs, mass_bcs)

# Create radiation cache
rad_cache = RadiationCache(mesh.n_active_cells)

# Define Marshak boundary conditions
radiation_bcs = RadiationBCSet(
	top = MarshakRadiationBC(I_external = 50e3, emissivity = 0.9),
	bottom = ReflectiveRadiationBC(),
)

# Create P1 radiation ODE function
ode_fn_p1 = PyrolysisODEFunctionALEWithP1Radiation(ode_fn, rad_cache, radiation_bcs)

# Solve
prob = ODEProblem(ode_fn_p1, u0, tspan)
sol = solve(prob, TRBDF2())
```
"""
mutable struct PyrolysisODEFunctionALEWithP1Radiation{ODE, RC, RBC <: RadiationBCSet}
	ode_fn::ODE
	rad_cache::RC
	radiation_bcs::RBC
	use_quasi_steady::Bool
	G_work::Vector{Float64}  # Work array for dG/dt (transient mode only)
end

"""
	PyrolysisODEFunctionALEWithP1Radiation(ode_fn, rad_cache, radiation_bcs;
										  use_quasi_steady=true)

Create a P1 radiation wrapper for ALE ODE function.

# Arguments
- `ode_fn`: Base PyrolysisODEFunctionALE to wrap
- `rad_cache`: RadiationCache with appropriate size
- `radiation_bcs`: RadiationBCSet for radiation boundary conditions
- `use_quasi_steady`: If true (default), solve quasi-steady P1 at each step

# Notes
- For quasi-steady mode (default), G is solved at each timestep
- For transient mode, G should be part of the state vector (use StateVectorALERadiation)
"""
function PyrolysisODEFunctionALEWithP1Radiation(ode_fn, rad_cache, radiation_bcs::RadiationBCSet;
	use_quasi_steady::Bool = true)
	n_cells = ode_fn.mesh.n_active_cells
	G_work = zeros(n_cells)

	PyrolysisODEFunctionALEWithP1Radiation(
		ode_fn, rad_cache, radiation_bcs, use_quasi_steady, G_work,
	)
end

"""
	(f::PyrolysisODEFunctionALEWithP1Radiation)(du, u, p, t)

Compute ALE RHS with P1 volumetric radiation.

This function:
1. Unpacks state and updates geometry (via wrapped ODE function)
2. Updates material and radiation properties
3. Solves quasi-steady P1 equation for G (if enabled)
4. Computes standard physics RHS
5. Adds P1 radiative flux divergence to energy equation
6. Packs results into du
"""
function (f::PyrolysisODEFunctionALEWithP1Radiation)(du, u, p, t)
	ode_fn = f.ode_fn
	mesh = ode_fn.mesh
	material = ode_fn.material
	sv = ode_fn.sv
	n = sv.n_cells
	NC = n_components(material)

	# =========================================================================
	# Step 1: Unpack state and update geometry
	# =========================================================================
	unpack_state_ale!(mesh, u, sv)

	# =========================================================================
	# Step 2: Update material properties
	# =========================================================================
	Properties.update_properties!(ode_fn.prop_cache, mesh, material)

	# =========================================================================
	# Step 3: Update radiation properties and solve quasi-steady P1
	# =========================================================================
	Discretization.update_radiation_properties!(f.rad_cache, mesh, ode_fn.prop_cache)

	if f.use_quasi_steady
		# Solve ∇·(D_rad∇G) = αG - 4σεαT⁴
		Discretization.solve_quasi_steady_radiation!(f.rad_cache, mesh, ode_fn.prop_cache,
			f.radiation_bcs, t)
	end

	# Compute radiative flux divergence: ∇·q_rad = α(G - 4σT⁴)
	Discretization.compute_radiation_flux_divergence!(f.rad_cache, mesh, ode_fn.prop_cache)

	# =========================================================================
	# Step 4: Compute standard physics RHS (conduction, reactions, transport)
	#         with I_radiation = 0 since P1 handles radiation
	# =========================================================================
	# Note: Pass I_radiation=0.0 directly instead of mutating ode_fn.I_radiation
	# This ensures autodiff compatibility (no mutation during function evaluation)

	# Choose between standard RHS and Darcy-Fick RHS based on component permeability values
	if has_darcy_flow(material)
		# Use component-based Darcy-Fick (no perm_model argument)
		compute_rhs_darcy_fick!(ode_fn.dT_work, ode_fn.dξ_work, mesh, ode_fn.prop_cache,
			ode_fn.rhs_cache, material, ode_fn.bc_set, t,
			ode_fn.top_face_ids, ode_fn.bottom_face_ids;
			I_radiation = 0.0)
	else
		compute_rhs!(ode_fn.dT_work, ode_fn.dξ_work, mesh, ode_fn.prop_cache,
			ode_fn.rhs_cache, material, ode_fn.bc_set, t,
			ode_fn.top_face_ids, ode_fn.bottom_face_ids;
			I_radiation = 0.0)
	end

	# =========================================================================
	# Step 5: Compute mesh velocity and ALE advection
	# =========================================================================
	_compute_ale_mesh_velocity!(ode_fn.w_work, mesh, ode_fn.prop_cache, material,
		t, ode_fn.min_density;
		θ_out = ode_fn.θ_work,
		precomputed_dξ_rxn = ode_fn.rhs_cache.dξ_rxn,
		dξ_dt = ode_fn.dξ_dt_work)

	# Interpolate node velocities to faces
	_interpolate_node_to_face!(ode_fn.w_face, ode_fn.w_work, mesh)

	# Compute ALE advection terms
	Discretization.compute_ale_advection_term!(ode_fn.ale_advection, mesh, ode_fn.w_work)
	Discretization.compute_ale_species_advection!(ode_fn.ale_species_advection, mesh,
		material, ode_fn.w_work)

	# Subtract ALE advection contributions
	_subtract_ale_advection!(ode_fn.dT_work, ode_fn.dξ_work,
		ode_fn.ale_advection, ode_fn.ale_species_advection, NC)

	# Apply ALE dilation correction for solid species
	_apply_ale_dilation_correction!(ode_fn.dξ_work, ode_fn.θ_work, mesh, material)

	# =========================================================================
	# Step 6: Add P1 radiation contribution to energy equation
	# =========================================================================
	# The P1 radiation model adds a volumetric source term to the energy equation:
	#   S_rad = α(G - 4σT⁴) = Q_absorption - Q_emission
	#   Positive = net absorption (heating), Negative = net emission (cooling)
	#
	# Note: Despite the name "div_q_rad", this is the NET RADIATION SOURCE, not
	# the divergence of flux. The naming comes from the steady-state balance
	# where ∇·q_rad + S_rad = 0, so S_rad = -∇·q_rad.
	#
	# Energy equation: ρcₚ dT/dt = (conduction) + (reactions) + S_rad
	@inbounds for i in 1:n
		ρcp = ode_fn.prop_cache.ρcp[i]
		# Add net radiation source (positive = heating)
		ode_fn.dT_work[i] += f.rad_cache.div_q_rad[i] / ρcp
	end

	# =========================================================================
	# Step 7: Pack results into du
	# =========================================================================
	n_nodes = sv.n_nodes

	# Temperature derivatives
	@inbounds for i in 1:n
		du[i] = ode_fn.dT_work[i]
	end

	# Concentration derivatives
	@inbounds for i in 1:n
		for j in 1:NC
			du[sv.ξ_offset+(j-1)*n+i] = ode_fn.dξ_work[j, i]
		end
	end

	# Node velocities (dz/dt = w)
	@inbounds for i in 1:n_nodes
		du[sv.z_offset+i] = ode_fn.w_work[i]
	end

	return nothing
end

"""
	(f::PyrolysisODEFunctionALEWithP1Radiation)(u, p, t)

Out-of-place P1 radiation ODE function.
"""
function (f::PyrolysisODEFunctionALEWithP1Radiation)(u, p, t)
	du = similar(u)
	f(du, u, p, t)
	return du
end

"""
	resize_ode_function_ale!(f::PyrolysisODEFunctionALEWithP1Radiation, n_cells_new::Int, n_faces_new::Int)

Resize the P1 radiation ODE function for mesh changes.
"""
function resize_ode_function_ale!(f::PyrolysisODEFunctionALEWithP1Radiation, n_cells_new::Int, n_faces_new::Int)
	# Resize wrapped ODE function
	resize_ode_function_ale!(f.ode_fn, n_cells_new, n_faces_new)

	# Resize radiation cache
	Properties.resize_radiation_cache!(f.rad_cache, n_cells_new)

	# Resize work array
	resize!(f.G_work, n_cells_new)

	return f
end

"""
	get_radiation_intensity(f::PyrolysisODEFunctionALEWithP1Radiation, cell_id::Int)

Get the radiation intensity G at a specific cell.
"""
function get_radiation_intensity(f::PyrolysisODEFunctionALEWithP1Radiation, cell_id::Int)
	return f.rad_cache.G[cell_id]
end

"""
	get_radiation_intensity(f::PyrolysisODEFunctionALEWithP1Radiation)

Get the full radiation intensity field G.
"""
function get_radiation_intensity(f::PyrolysisODEFunctionALEWithP1Radiation)
	return f.rad_cache.G
end
