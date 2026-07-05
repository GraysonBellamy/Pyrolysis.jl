```@meta
CurrentModule = Pyrolysis
```

# Residual layer

The pure-residual surface underneath `solve`: a [`ProblemDef`](@ref) is an
immutable problem description, a [`Workspace`](@ref) holds the mutable
buffers, and [`residual!`](@ref) evaluates the semi-discrete right-hand side.
This is the layer to target for custom solver loops; see the User Guide's
[solver-configuration chapter](../guide/user_guide/11_solver-jacobian-config.md).

## Problem definition and workspaces

```@docs
SimulationMode
StateLayout
ProblemDef
Workspace
to_problem_def
build_workspace
build_ad_workspace
residual!
```

## State-vector layout and accessors

```@docs
state_length
create_state_vector
pack_state!
update_mesh_from_state!
block_view
get_temperature
get_concentration
get_node_position
get_node_positions
get_thickness_from_state
column_chi_bar_from_state
cross_section_area_from_state
```
