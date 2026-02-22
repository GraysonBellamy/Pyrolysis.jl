# Hybrid Discretization Framework for Pyrolysis Modeling

## From-Scratch Design: Vertex-Centered Temperature, Cell-Centered Species

**Version**: 1.0  
**Purpose**: Comprehensive specification for a physically optimal discretization scheme

---

## Executive Summary

This document describes a **hybrid discretization approach** that combines:
- **Vertex-centered temperature** (T at nodes) — for exact boundary conditions
- **Cell-centered species concentrations** (ξⱼ in cells) — for mass conservation

This approach is physically optimal for pyrolysis modeling where:
1. Boundary conditions directly involve surface temperature
2. Mass conservation is critical for species tracking
3. ALE mesh motion naturally tracks node positions

---

## 1. Mathematical Framework

### 1.1 Mesh Topology

```
    Node 0      Cell 0      Node 1      Cell 1      Node 2      Cell 2      Node 3
      •    |____________|     •    |____________|     •    |____________|     •
     z=0        V₀           z₁         V₁           z₂         V₂          z=L
   (surface)                                                              (back)
   
   T₀, T₁, T₂, T₃ — temperatures at nodes (N+1 values)
   ξ₀, ξ₁, ξ₂     — concentrations in cells (N values)
```

**Key relationships:**
- N cells with volumes V₀, V₁, ..., V_{N-1}
- N+1 nodes with positions z₀, z₁, ..., z_N
- Cell i spans from node i to node i+1
- Cell volume: V_i = (z_{i+1} - z_i) × A (where A = cross-sectional area)

### 1.2 State Variables

**Primary variables (in ODE state vector):**
```julia
struct HybridState{NC}
    # Vertex-centered (N+1 values)
    T::Vector{Float64}           # Temperature at each node [K]
    
    # Cell-centered (N values, NC components each)
    ξ::Vector{NTuple{NC,Float64}} # Mass concentrations [kg/m³]
    
    # Node positions (for ALE)
    z::Vector{Float64}           # Node positions [m]
end
```

**Total DOFs:** (N+1) + N×NC + (N+1) = 2N + 2 + N×NC

### 1.3 Derived Quantities

For each **cell i**:
```julia
# Cell-averaged temperature (for properties and reactions)
T_cell_i = (T[i] + T[i+1]) / 2

# Or more accurate weighted average for non-uniform mesh:
T_cell_i = (T[i] * Δz_R + T[i+1] * Δz_L) / (Δz_L + Δz_R)
where Δz_L = z_center - z[i], Δz_R = z[i+1] - z_center

# Total density
ρ_i = sum(ξ_i)

# Mass fractions  
Y_j,i = ξ_j,i / ρ_i

# Volume fractions
v_j,i = Y_j,i × ρ_i / ρ_j(T_cell_i)
```

For each **node i** (interior, i = 1, ..., N-1):
```julia
# Control volume for energy equation (dual mesh)
# Extends from midpoint of left cell to midpoint of right cell
V_dual_i = A × ((z[i] - z[i-1])/2 + (z[i+1] - z[i])/2)
         = A × (z[i+1] - z[i-1]) / 2

# Effective thermal mass at node
(ρc)_node_i = ((ρc)_cell_{i-1} + (ρc)_cell_i) / 2
```

For **boundary nodes** (i = 0 and i = N):
```julia
# Half control volume at boundaries
V_dual_0 = A × (z[1] - z[0]) / 2
V_dual_N = A × (z[N] - z[N-1]) / 2
```

---

## 2. Governing Equations

### 2.1 Energy Equation (at nodes)

For **interior node i** (i = 1, ..., N-1):
```
(ρc)_i × V_dual_i × dT_i/dt = A × (q_{i-1/2} - q_{i+1/2}) + Q_source_i
```

Where:
- `q_{i-1/2}` = heat flux from cell i-1 to node i
- `q_{i+1/2}` = heat flux from node i to cell i
- `Q_source_i` = heat source at node (from reactions in adjacent cells)

**Heat flux between nodes** (through cell i):
```
q_i = -k_i × (T_{i+1} - T_i) / (z_{i+1} - z_i)

where k_i = k(T_cell_i, ξ_i)  # Conductivity evaluated at cell center
```

**Source term at node** (contribution from adjacent cells):
```
Q_source_i = (V_{i-1} × Σⱼ h_j × r_j,{i-1} + V_i × Σⱼ h_j × r_j,i) / 2
```

For **boundary node 0** (surface):
```
(ρc)_0 × V_dual_0 × dT_0/dt = A × (q_BC - q_{1/2}) + Q_source_0
```

Where `q_BC` is the boundary heat flux:
```
q_BC = q_ext - ε × σ × T_0⁴ - h_c × (T_0 - T_∞)
```

**Key advantage:** T_0 IS the surface temperature — no solve needed!

For **boundary node N** (back face):
```
(ρc)_N × V_dual_N × dT_N/dt = A × (q_{N-1/2} - q_BC,back) + Q_source_N
```

### 2.2 Species Equations (in cells)

For **cell i** (i = 0, ..., N-1):
```
V_i × dξ_j,i/dt = Σ_reactions ν_ij × r_i × V_i - A × (J_j,{i+1/2} - J_j,{i-1/2})
```

Where:
- `ν_ij` = stoichiometric coefficient for species j in reaction
- `r_i` = reaction rate in cell i (evaluated at T_cell_i)
- `J_j` = mass flux of species j

**Reaction rate** (Arrhenius, evaluated at cell-averaged T):
```
r_i = A_rxn × exp(-E / (R × T_cell_i)) × ξ_reactant,i^n
```

**Species flux at cell face** (between cells i-1 and i):
```
# For gaseous species (diffusive transport)
J_j,{i-1/2} = -λ_j × ρ × ∂(ξ_j/ρ)/∂z |_{i-1/2}
            ≈ -λ_{i-1/2} × ρ_{i-1/2} × (Y_j,i - Y_j,{i-1}) / (z_center_i - z_center_{i-1})

# For solid species (no flux)
J_j = 0
```

**Boundary flux for species** (at surface):
```
J_j,{-1/2} = h_m × (ξ_j,0 - ξ_j,∞)  # Convective BC
           = 0                       # Impermeable BC
```

### 2.3 Mesh Motion (ALE)

For **node i**:
```
dz_i/dt = w_i    # Mesh velocity
```

**Surface node** (i = 0):
```
# Surface regresses based on mass loss
dz_0/dt = -ṁ_surface / ρ_surface
        = -Σⱼ J_j,surface / ρ_0
```

**Interior nodes**: Mesh velocity from smoothing (Laplacian, spring analogy, etc.)

**Back node** (i = N): Fixed or prescribed motion

---

## 3. Property Evaluation Strategy

### 3.1 Where to Evaluate Properties

| Property | Location | Formula |
|----------|----------|---------|
| k (conductivity) | Cell center | k(T_cell_i, ξ_i) with T_cell = (T_i + T_{i+1})/2 |
| c_p (heat capacity) | Cell center | c_p(T_cell_i, ξ_i) |
| ρ (density) | Cell center | From ξ_i (it's stored there) |
| λ (gas transfer) | Cell face | Interpolate from cells |
| ε (emissivity) | Boundary node | ε(T_0, ξ_0) with ξ_0 from surface cell |
| r (reaction rate) | Cell center | r(T_cell_i, ξ_i) |

### 3.2 Handling Nonlinear Properties

For conductivity mixing at cell center:
```julia
function effective_conductivity(T_cell, ξ, material, β)
    # Volume fractions from cell concentrations
    v = compute_volume_fractions(ξ, T_cell, material)
    
    k_components = [comp.k(T_cell) for comp in material.components]
    
    k_parallel = sum(v .* k_components)
    k_series = 1 / sum(v ./ k_components)
    
    return β * k_parallel + (1 - β) * k_series
end
```

### 3.3 Face Values (for fluxes between cells)

```julia
# Harmonic mean for conductivity at face between cells i and i+1
k_face = 2 * k_i * k_{i+1} / (k_i + k_{i+1})

# Arithmetic mean for other properties
λ_face = (λ_i + λ_{i+1}) / 2
```

---

## 4. Discretization Details

### 4.1 Energy Equation Assembly

```julia
function compute_energy_rhs!(dT, T, ξ, z, material, bcs, t)
    N = length(T) - 1  # Number of cells
    A = cross_sectional_area
    
    # Compute cell properties
    for i in 0:N-1
        T_cell = (T[i] + T[i+1]) / 2
        k[i] = effective_conductivity(T_cell, ξ[i], material)
        ρc[i] = effective_heat_capacity(T_cell, ξ[i], material)
        Q_rxn[i] = reaction_heat_source(T_cell, ξ[i], material)
    end
    
    # Compute fluxes through each cell
    for i in 0:N-1
        Δz = z[i+1] - z[i]
        q[i] = -k[i] * (T[i+1] - T[i]) / Δz
    end
    
    # Surface node (i = 0)
    q_BC = compute_surface_bc_flux(T[0], bcs.top, t)
    V_dual_0 = A * (z[1] - z[0]) / 2
    ρc_0 = ρc[0]  # Use surface cell properties
    dT[0] = (A * (q_BC - q[0]) + Q_rxn[0]/2) / (ρc_0 * V_dual_0)
    
    # Interior nodes
    for i in 1:N-1
        V_dual_i = A * (z[i+1] - z[i-1]) / 2
        ρc_i = (ρc[i-1] + ρc[i]) / 2
        Q_i = (Q_rxn[i-1] + Q_rxn[i]) / 2
        dT[i] = (A * (q[i-1] - q[i]) + Q_i) / (ρc_i * V_dual_i)
    end
    
    # Back node (i = N)
    q_BC_back = compute_back_bc_flux(T[N], bcs.bottom, t)
    V_dual_N = A * (z[N] - z[N-1]) / 2
    ρc_N = ρc[N-1]
    dT[N] = (A * (q[N-1] - q_BC_back) + Q_rxn[N-1]/2) / (ρc_N * V_dual_N)
end
```

### 4.2 Species Equation Assembly

```julia
function compute_species_rhs!(dξ, T, ξ, z, material, bcs, t)
    N = length(ξ)  # Number of cells
    A = cross_sectional_area
    NC = n_components(material)
    
    for i in 0:N-1
        V_i = A * (z[i+1] - z[i])
        T_cell = (T[i] + T[i+1]) / 2
        
        # Reaction source terms
        r_source = compute_reaction_rates(T_cell, ξ[i], material)
        
        # Gas transport fluxes (at cell faces)
        if i == 0
            J_left = compute_surface_mass_flux(ξ[0], bcs.mass_top, t)
        else
            J_left = compute_interior_mass_flux(ξ[i-1], ξ[i], z, i-1, material)
        end
        
        if i == N-1
            J_right = compute_back_mass_flux(ξ[N-1], bcs.mass_bottom, t)
        else
            J_right = compute_interior_mass_flux(ξ[i], ξ[i+1], z, i, material)
        end
        
        for j in 1:NC
            dξ[i][j] = r_source[j] - (A / V_i) * (J_right[j] - J_left[j])
        end
    end
end
```

### 4.3 Boundary Condition Implementation

```julia
function compute_surface_bc_flux(T_s, bc::CombinedThermalBC, t)
    q_total = 0.0
    
    for sub_bc in bc.components
        q_total += compute_bc_component(T_s, sub_bc, t)
    end
    
    return q_total
end

function compute_bc_component(T_s, bc::HeatFluxBC, t)
    return bc.flux(t)  # Positive = into domain
end

function compute_bc_component(T_s, bc::ConvectiveBC, t)
    return bc.h * (bc.T_inf(t) - T_s)  # Positive = into domain
end

function compute_bc_component(T_s, bc::RadiativeBC, t)
    return bc.ε * σ * (bc.T_inf(t)^4 - T_s^4)  # Positive = into domain
end
```

---

## 5. Conservation Verification

### 5.1 Energy Conservation

Total energy in domain:
```julia
E_total = sum(ρc[i] * V[i] * T_cell[i] for i in 0:N-1)
        = sum(ρc[i] * V[i] * (T[i] + T[i+1])/2 for i in 0:N-1)
```

Rate of change:
```julia
dE/dt = q_surface * A - q_back * A + sum(Q_rxn[i] for i in 0:N-1)
```

This should balance exactly (to machine precision) with proper implementation.

### 5.2 Mass Conservation

Total mass of species j:
```julia
M_j = sum(ξ_j[i] * V[i] for i in 0:N-1)
```

Rate of change:
```julia
dM_j/dt = J_j,surface * A - J_j,back * A + sum(ν_j * r * V for all reactions)
```

---

## 6. Data Structures

### 6.1 Type Definitions

```julia
"""
Hybrid mesh with vertex-centered T and cell-centered ξ.
"""
mutable struct HybridMesh1D{NC}
    # Topology
    n_cells::Int
    n_nodes::Int  # = n_cells + 1
    
    # Node data
    z::Vector{Float64}              # Node positions [m]
    z_ref::Vector{Float64}          # Reference positions [m]
    T::Vector{Float64}              # Temperature at nodes [K]
    node_boundary_id::Vector{Int}   # 0=interior, 1=surface, 2=back
    
    # Cell data
    ξ::Vector{NTuple{NC,Float64}}   # Concentrations [kg/m³]
    V::Vector{Float64}              # Cell volumes [m³]
    V_ref::Vector{Float64}          # Reference volumes [m³]
    
    # Precomputed geometry
    Δz::Vector{Float64}             # Cell sizes [m]
    z_center::Vector{Float64}       # Cell center positions [m]
    V_dual::Vector{Float64}         # Dual control volumes at nodes [m³]
    
    # ALE state
    w::Vector{Float64}              # Node velocities [m/s]
end

"""
Property cache for hybrid mesh.
"""
mutable struct HybridPropertyCache{NC}
    # Cell-centered (N values)
    T_cell::Vector{Float64}         # Cell-averaged temperature
    k::Vector{Float64}              # Effective conductivity
    ρc::Vector{Float64}             # Volumetric heat capacity
    ρ::Vector{Float64}              # Total density
    v::Vector{NTuple{NC,Float64}}   # Volume fractions
    
    # For reactions
    r::Vector{Float64}              # Total reaction rate
    Q_rxn::Vector{Float64}          # Reaction heat source
    
    # Node-centered (N+1 values)
    ρc_node::Vector{Float64}        # Interpolated to nodes
    
    # Face values (N values, between cells)
    k_face::Vector{Float64}         # Harmonic mean conductivity
    q::Vector{Float64}              # Heat flux through cell
end
```

### 6.2 Mesh Generation

```julia
function create_hybrid_mesh(L::Real, n_cells::Int, NC::Int;
                            stretching::Symbol=:uniform,
                            stretch_ratio::Real=1.0)
    n_nodes = n_cells + 1
    
    # Generate node positions
    if stretching == :uniform
        z = collect(range(0.0, L, length=n_nodes))
    elseif stretching == :geometric
        z = geometric_spacing(0.0, L, n_nodes, stretch_ratio)
    end
    
    # Initialize temperatures (uniform)
    T = fill(300.0, n_nodes)
    
    # Initialize concentrations (will be set by initial conditions)
    ξ = [ntuple(_ -> 0.0, NC) for _ in 1:n_cells]
    
    # Compute geometry
    Δz = [z[i+1] - z[i] for i in 1:n_cells]
    z_center = [(z[i] + z[i+1]) / 2 for i in 1:n_cells]
    V = [Δz[i] * A for i in 1:n_cells]  # A = unit area for 1D
    
    # Dual control volumes at nodes
    V_dual = zeros(n_nodes)
    V_dual[1] = Δz[1] / 2 * A
    for i in 2:n_nodes-1
        V_dual[i] = (Δz[i-1] + Δz[i]) / 2 * A
    end
    V_dual[n_nodes] = Δz[n_cells] / 2 * A
    
    # Boundary IDs
    node_boundary_id = zeros(Int, n_nodes)
    node_boundary_id[1] = 1      # Surface
    node_boundary_id[n_nodes] = 2  # Back
    
    return HybridMesh1D{NC}(
        n_cells, n_nodes,
        z, copy(z), T, node_boundary_id,
        ξ, V, copy(V),
        Δz, z_center, V_dual,
        zeros(n_nodes)  # Mesh velocities
    )
end
```

---

## 7. Time Integration

### 7.1 ODE Function

```julia
function hybrid_rhs!(du, u, p, t)
    mesh, material, bcs, cache = p
    
    # Unpack state
    T = @view u[1:mesh.n_nodes]
    ξ_flat = @view u[mesh.n_nodes+1:end]
    
    # Update properties at current state
    update_property_cache!(cache, mesh, material)
    
    # Compute energy RHS
    dT = @view du[1:mesh.n_nodes]
    compute_energy_rhs!(dT, mesh, cache, bcs, t)
    
    # Compute species RHS
    dξ_flat = @view du[mesh.n_nodes+1:end]
    compute_species_rhs!(dξ_flat, mesh, cache, bcs, t)
end
```

### 7.2 State Vector Layout

```julia
# State vector u has length: (N+1) + N×NC
# Layout:
#   u[1:N+1]           — T at nodes
#   u[N+2:N+1+N×NC]    — ξ in cells (flattened)

function unpack_state(u, n_nodes, n_cells, NC)
    T = @view u[1:n_nodes]
    ξ_flat = @view u[n_nodes+1:n_nodes+n_cells*NC]
    ξ = reshape(ξ_flat, NC, n_cells)  # NC × N_cells matrix
    return T, ξ
end

function pack_state!(u, T, ξ)
    n_nodes = length(T)
    u[1:n_nodes] .= T
    u[n_nodes+1:end] .= vec(ξ)
end
```

---

## 8. Advantages of This Approach

### 8.1 Physics Benefits

1. **Exact surface temperature**: T[0] IS T_surface, no approximation
2. **Direct BC application**: q_BC uses T[0] directly
3. **Proper mass conservation**: ξ in cells with cell-integrated equations
4. **Natural ALE**: Node positions are primary variables

### 8.2 Numerical Benefits

1. **Higher-order boundary gradients**: Can use (-3T_0 + 4T_1 - T_2)/(2Δz)
2. **No boundary temperature solve**: Eliminates implicit coupling
3. **Clear output**: Surface T is state variable, not derived

### 8.3 Validation Benefits

1. **Direct ThermaKin comparison**: ThermaKin reports T at z=0
2. **Unambiguous surface temperature**: No confusion about cell-center vs boundary

---

## 9. Implementation Complexity

### 9.1 Additional Bookkeeping

- Dual control volumes for energy equation
- Cell-averaged T for property evaluation
- Interpolation of some quantities between nodes and cells

### 9.2 Stencil Considerations

- Energy equation at interior node involves: T_{i-1}, T_i, T_{i+1}
- Same as cell-centered approach
- No increase in stencil width

### 9.3 Jacobian Structure

- Temperature block: tridiagonal (same as cell-centered)
- Species block: tridiagonal (same as cell-centered)
- Coupling: T affects species through reactions, species affect T through properties
- Overall structure similar to cell-centered approach

---

## 10. Migration Path

If implementing in existing codebase:

1. **Phase 1**: Add node-based T storage alongside cell-based
2. **Phase 2**: Modify energy equation to use node-based T
3. **Phase 3**: Update boundary condition implementation
4. **Phase 4**: Update output routines
5. **Phase 5**: Remove cell-based T storage

Total effort: Significant but incremental migration possible.
