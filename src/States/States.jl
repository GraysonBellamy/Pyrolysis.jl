"""
	States

The state of the finite-volume core and the mapping from it to the quantities
the kernels read.

The integrator carries per-cell component masses, an energy state (the cell
temperature, or the enthalpy inventory behind the same two functions), and on
a moving mesh the node positions and velocities, with the face intensities and
cross-section scalars where the physics needs them. [`derive!`](@ref) is the
only reader of the raw state: it maps the masses and node positions to
volumes and concentrations, the energy state to temperatures, and evaluates
the cell properties, every trial-state quantity through a smooth guard.
"""
module States

using LinearAlgebra
using StaticArrays
using ForwardDiff
using ..Materials
using ..Materials:
    ConstantProperty, LinearProperty, PolynomialProperty, TableProperty, T_REF
using ..Properties: effective_conductivity
using ..Meshes
import ..Materials: n_components
import ..Meshes: n_cells, n_nodes

include("guards.jl")
include("layout.jl")
include("caloric.jl")
include("derived.jl")
include("tolerances.jl")
include("fields.jl")

export smooth_clamp, trial_temperature, positive_part, smooth_floor
export EnergyState, TemperatureState, EnthalpyState
export StateLayout, state_length, energy_state, is_moving, has_intensities, has_chi
export n_components, node_width, cell_width, cell_rows, node_rows
export mass_index, energy_index, history_index, position_index, velocity_index
export intensity_index, area_index, lateral_rate_index
export is_differential, mass_matrix
export sensible_enthalpy, enthalpy_slope, enthalpy_law, TabulatedEnthalpy, CaloricLaws
export component_enthalpies, component_enthalpy_slopes, component_heat_capacities
export CaloricExtension, regularization, extended_energy, extended_capacity
export invert_caloric
export CellReference, StateSetup, DerivedState, derive!, energy_state_value
export abstol_vector
export CarriedField, CONVERSION_HISTORY, material_axes, physical_tensor

end # module States
