"""
	ALE

Arbitrary Lagrangian-Eulerian mesh motion.

This module provides:
- Mesh velocity computation (Laplacian smoothing, spring analogy, weighted ALE)
- Mesh motion algorithms with GCL compliance
- Geometric Conservation Law verification
- Conservative solution remapping
- Interface tracking (char/virgin boundary, reaction fronts)

Reference: Design Document Sections 3 (ALE Formulation), 6.5 (ALE Types), 7.3 (Mesh Motion)
"""
module ALE

using ..Core
using ..Geometry
using ..Properties
using LinearAlgebra
using SparseArrays

# Include subfiles
include("mesh_velocity.jl")
include("mesh_motion.jl")
include("gcl.jl")
include("remapping.jl")
include("interface_tracking.jl")

# Exports - Mesh Velocity Methods
export LaplacianSmoothing, SpringAnalogy, LinearElasticity, WeightedALE
export compute_mesh_velocity!, compute_mesh_velocity_from_volume_change!
export surface_regression_rate

# Exports - Mesh Motion
export move_mesh!, move_mesh_safe!
export check_mesh_quality, mesh_is_valid
export mesh_displacement, max_mesh_displacement, reset_mesh_to_reference!
export update_face_positions!, update_cell_centers!

# Exports - GCL
export check_gcl, gcl_error, gcl_error_post_motion
export gcl_volume_update, update_volumes_gcl!

# Exports - Remapping
export FirstOrderRemap, SecondOrderRemap
export remap_solution!, remap_cell_states!
export check_remap_conservation

# Exports - Interface Tracking
export InterfaceTracker
export find_interface_position, update_interface!, update_interface_from_state!
export interface_velocity, interface_velocity_from_flux
export interface_tracking_velocity!

end # module ALE
