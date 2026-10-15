"""
	Pyrolysis.Experimental

Namespace for experimental or unstable code paths. The names below remain
functional and are accessible as `Pyrolysis.Experimental.X`, but they are
deliberately removed from the top-level public API. These interfaces may
change between releases; prefer the stable top-level exports when possible.

Contents:

- **P1 quasi-steady radiation** (residual wrapper + G kernels). Lives in
  `RadiationP1.jl` — built on top of the unified `(def, ws)` pair and the
  unified `StateLayout`; the quasi-steady G field lives in the
  `RadiationCache`, not in the ODE state vector.

- **h-AMR** (h-refinement, h-coarsening, AMR controller, adjacent-cell merging,
  surface-depletion utilities). Will move to a dedicated `AdaptivityH`
  extension once the active-cell desync invariants are reworked.
"""
module Experimental

using ..Materials
using ..BoundaryConditions
using ..Geometry
using ..Properties
using ..Physics
using ..Discretization
using ..Adaptivity
using ..Problem
using ..Residual
import ..Materials: n_components
# Unexported BoundaryConditions helper used by the moved radiation_boundary.jl.
using ..BoundaryConditions: _get_scalar_value
using LinearAlgebra: norm
using SparseArrays: spzeros
using ForwardDiff

# ============================================================================
# P1 quasi-steady radiation
# ============================================================================
#
# The P1 implementation lives here, separate from the core residual: the
# radiation BC types, the P1 property/flux/BC kernels, and the quasi-steady G
# solver + net volumetric source. Include order matters: BC types first
# (referenced in kernel signatures), then the kernels, then the (def, ws)
# wrapper. Every kernel has exactly one method — G lives in the
# RadiationCache, never in the ODE state vector.
include("radiation_boundary.jl")
include("radiation_p1.jl")
include("radiation_fluxes.jl")
include("RadiationP1.jl")

# Public radiation-mode enum value (the type stays in Core for trait dispatch;
# only the public-API export moves here).
export P1_QUASI_STEADY

# Marshak / reflective / Dirichlet BCs (used only by P1)
export AbstractRadiationBC, MarshakRadiationBC, ReflectiveRadiationBC, DirichletRadiationBC
export RadiationBCSet, get_radiation_bc, default_radiation_bcs

# Radiation cache + diagnostics
export RadiationCache, resize_radiation_cache!, initialize_radiation_intensity!
export radiation_energy_density, net_radiation_source, total_radiation_flux

# P1 kernels / diagnostics
export update_radiation_properties!, compute_radiation_flux_divergence!
export solve_quasi_steady_radiation!, initialize_g_from_temperature!
export compute_radiation_diffusion_coefficient
export compute_radiation_source_terms!, compute_radiation_fluxes!
export radiation_bc_flux
export optical_thickness, optical_thickness_cell, check_p1_validity
export radiation_energy_balance
export has_absorption_coefficients, requires_radiation_cache

# ============================================================================
# h-AMR
# ============================================================================
#
# The h-refinement / h-coarsening / AMR-controller cluster lives here. The
# depleted-cell park (`park_cell!`) and the find_* neighbour helpers stay in
# core `Adaptivity` (used by the core solve callback). The AMR-coarsening
# path borrows `merge_temperatures_conservatively` from core via
# `using ..Adaptivity`.

include("h_refinement.jl")
include("h_coarsening.jl")
include("amr_controller.jl")
include("hr_adapt.jl")

export AMRController
export adapt_mesh!, should_adapt, maybe_adapt!
export amr_statistics, reset_statistics!, create_default_controller

export refine_cell!, perform_h_refinement!
export mark_cells_for_refinement, enforce_2to1_balance!
export verify_refinement_conservation

export coarsen_cells!, perform_h_coarsening!
export mark_cells_for_coarsening
export can_coarsen, get_sibling
export verify_coarsening_conservation

export hr_adapt!

end # module Experimental
