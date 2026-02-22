# Pyrolysis.jl Implementation Plan

## Staged Implementation Guide for AI-Assisted Development

**Reference Document**: `Pyrolysis_ALE_Framework_Design.md`
**Target**: Julia 1.9+
**Approach**: Bottom-up, test-driven, modular implementation

---

## Overview

This document provides a staged implementation plan for the Pyrolysis.jl package. Each stage builds upon the previous one, with clear deliverables, test requirements, and verification steps. Stages are designed to be implementable in focused sessions with well-defined interfaces.

### Implementation Principles

1. **Test-First**: Write tests before or alongside implementation
2. **Type Stability**: Verify with `@code_warntype` for all hot paths
3. **Incremental Verification**: Each stage includes verification against analytical solutions
4. **Clean Interfaces**: Define and document public APIs before implementation
5. **Performance Awareness**: Avoid allocations in inner loops from the start

---

## Stage 0: Project Scaffolding

**Goal**: Set up the package structure and development environment

### Deliverables

1. Create package directory structure per Section 5.1 of the design doc:
   ```
   Pyrolysis.jl/
   ├── Project.toml
   ├── src/
   │   ├── Pyrolysis.jl          # Main module
   │   ├── Core/
   │   ├── Geometry/
   │   ├── Properties/
   │   ├── Physics/
   │   ├── ALE/
   │   ├── Adaptivity/
   │   ├── Discretization/
   │   ├── TimeIntegration/
   │   ├── Problem/
   │   ├── Solvers/
   │   └── PostProcess/
   ├── test/
   │   └── runtests.jl
   └── examples/
   ```

2. Create `Project.toml` with dependencies (Section 13.1):
   - LinearAlgebra, SparseArrays (stdlib)
   - StaticArrays
   - DifferentialEquations (or OrdinaryDiffEq + Sundials)
   - ForwardDiff
   - LinearSolve
   - RecipesBase (optional)

3. Set up main module `src/Pyrolysis.jl` with submodule includes

4. Create `test/runtests.jl` with basic test infrastructure

### Acceptance Criteria

- [ ] Package instantiates without errors: `using Pyrolysis` works
- [ ] Tests run: `Pkg.test()` passes (empty test suite OK)
- [ ] All submodule files exist (can be empty stubs)

---

## Stage 1: Core Types and Constants

**Goal**: Implement fundamental types that all other modules depend on

**Reference**: Design Doc Sections 6.1-6.2 (Type System), Appendix A (Constants)

### 1.1 Physical Constants (`src/Core/constants.jl`)

Implement the `Constants` module:
```julia
module Constants
    const R_GAS = 8.31446261815324           # J/(mol·K)
    const STEFAN_BOLTZMANN = 5.670374419e-8  # W/(m²·K⁴)
    const T_REF = 298.15                     # K
    const P_REF = 101325.0                   # Pa
end
```

### 1.2 Property Functions (`src/Core/property_functions.jl`)

Implement all property function types from Section 6.2:

1. `ConstantProperty{T}`
2. `LinearProperty{T}`
3. `PolynomialProperty{N,T}`
4. `ArrheniusProperty{T}`
5. `TableProperty{N,T}` (with linear interpolation)
6. `CallableProperty{F}`

Each must implement the callable interface: `(p::PropertyType)(T) -> value`

**Key Implementation Notes**:
- Use `NTuple` for fixed-size data (type stable)
- Implement `linear_interpolate` helper for TableProperty
- Ensure all are allocation-free when called

### 1.3 Phase Enum and Component Types (`src/Core/components.jl`)

```julia
@enum Phase SOLID=1 LIQUID=2 GAS=3

struct Component{Pρ,Pc,Pk,Pλ} <: AbstractComponent
    name::Symbol
    phase::Phase
    molar_mass::Float64
    density::Pρ
    heat_capacity::Pc
    thermal_conductivity::Pk
    gas_transfer_coeff::Pλ
    emissivity::Float64
    absorption_coeff::Float64
    swelling_factor::Float64
end
```

Implement convenience constructors:
- `SolidComponent(name; ρ, c, k, ε=0.9, α=0.0, γ=0.0)`
- `GaseousComponent(name; M, c, k, λ, ε=0.0, α=0.0)`

### 1.4 Reaction Types (`src/Core/reactions.jl`)

```julia
struct Reaction{Ph} <: AbstractReaction
    name::Symbol
    reactant1::Int
    reactant2::Int
    order1::Float64
    order2::Float64
    product1::Int
    product2::Int
    stoich::NTuple{4,Float64}
    A::Float64
    E::Float64
    heat::Ph
    T_min::Float64
    T_max::Float64
end

struct ReactionSet{NR,Ph}
    reactions::NTuple{NR,Reaction{Ph}}
end
```

### 1.5 Material Type (`src/Core/materials.jl`)

```julia
@enum MixingRule PARALLEL=1 SERIES=2 WEIGHTED=3

struct Material{NC,NR,C,R} <: AbstractMaterial
    name::Symbol
    components::NTuple{NC,C}
    reactions::NTuple{NR,R}
    conductivity_mixing::MixingRule
    β::Float64
    gas_transfer_mixing::MixingRule
    β_λ::Float64
end
```

Implement accessor functions: `n_components`, `n_reactions`, `n_solids`, `n_gases`

### 1.6 Abstract Type Hierarchy (`src/Core/Core.jl`)

Define all abstract types from Section 6.1:
- `AbstractPropertyFunction`, `AbstractThermalProperty`, `AbstractTransportProperty`
- `AbstractComponent`, `AbstractReaction`, `AbstractMaterial`
- `AbstractNode`, `AbstractCell`, `AbstractFace`, `AbstractMesh{DIM}`
- `AbstractBoundaryCondition`, `AbstractThermalBC`, `AbstractMassBC`, `AbstractMeshBC`
- etc.

### Tests for Stage 1 (`test/test_core.jl`)

1. Test each property function type:
   - Correct evaluation at various temperatures
   - Type stability: `@inferred prop(300.0)`
   - Zero allocations: `@allocated prop(300.0) == 0`

2. Test component creation and property access

3. Test reaction stoichiometry (mass balance checks)

4. Test material creation with multiple components

### Acceptance Criteria

- [ ] All property functions are type-stable and allocation-free
- [ ] Components can be created with various property types
- [ ] Materials with N components and M reactions compile to concrete types
- [ ] All tests pass

---

## Stage 2: Geometry Module (1D Focus)

**Goal**: Implement 1D mesh infrastructure

**Reference**: Design Doc Section 6.3 (Geometry Types), Section 5.2 (Mesh Structure)

### 2.1 Node Types (`src/Geometry/nodes.jl`)

```julia
struct Node1D <: AbstractNode
    id::Int
    z::Float64
    z_ref::Float64
    is_boundary::Bool
    boundary_id::Int  # 0=interior, 1=top, 2=bottom
end
```

### 2.2 Cell Types (`src/Geometry/cells.jl`)

```julia
struct Cell1D <: AbstractCell
    id::Int
    node_ids::Tuple{Int,Int}
    face_ids::Tuple{Int,Int}
    level::Int
    parent::Int
    children::Tuple{Int,Int}
end

mutable struct CellState{NC}
    T::Float64
    ξ::NTuple{NC,Float64}
    volume::Float64
    volume_ref::Float64
end
```

### 2.3 Face Types (`src/Geometry/faces.jl`)

```julia
struct Face1D <: AbstractFace
    id::Int
    node_id::Int
    cell_ids::Tuple{Int,Int}  # (left, right), 0 = boundary
    is_boundary::Bool
    boundary_id::Int
end
```

### 2.4 Mesh Container (`src/Geometry/mesh.jl`)

```julia
mutable struct Mesh1D <: AbstractMesh{1}
    nodes::Vector{Node1D}
    cells::Vector{Cell1D}
    faces::Vector{Face1D}
    cell_states::Vector{CellState}

    # Geometry cache
    cell_centers::Vector{Float64}
    cell_volumes::Vector{Float64}
    face_positions::Vector{Float64}
    face_areas::Vector{Float64}

    # Connectivity
    boundary_faces::Dict{Int,Vector{Int}}

    # Refinement tracking
    max_level::Int
    n_active_cells::Int
    active_cell_mask::BitVector

    # ALE state
    node_velocities::Vector{Float64}
    mesh_time::Float64
end
```

### 2.5 Mesh Generation (`src/Geometry/mesh_generation.jl`)

Implement mesh creation functions:

```julia
function create_uniform_mesh(L::Float64, n_cells::Int, NC::Int) -> Mesh1D
function create_stretched_mesh(L::Float64, n_cells::Int, stretch::Float64, NC::Int) -> Mesh1D
```

Parameters:
- `L`: Domain length (thickness)
- `n_cells`: Number of cells
- `stretch`: Geometric stretching factor (1.0 = uniform)
- `NC`: Number of components (for CellState sizing)

### 2.6 Geometric Quantities (`src/Geometry/geometric_quantities.jl`)

```julia
function update_geometry!(mesh::Mesh1D)
    # Recompute all cached geometric quantities from node positions
end

function cell_volume(mesh::Mesh1D, i::Int) -> Float64
function face_area(mesh::Mesh1D, f::Int) -> Float64  # Cross-sectional area
function cell_center(mesh::Mesh1D, i::Int) -> Float64
```

### Tests for Stage 2 (`test/test_mesh.jl`)

1. Create uniform mesh and verify:
   - Correct number of nodes, cells, faces
   - Cell volumes sum to total domain volume
   - Boundary faces correctly identified

2. Create stretched mesh and verify:
   - Stretching ratio between adjacent cells
   - Total volume conservation

3. Test geometry update after node movement

4. Type stability checks for all mesh operations

### Acceptance Criteria

- [ ] Uniform and stretched meshes can be created
- [ ] Geometry quantities are correctly computed
- [ ] Boundary connectivity is correct
- [ ] All operations type-stable

---

## Stage 3: Properties Module

**Goal**: Compute mixture properties from component data

**Reference**: Design Doc Section 2.4 (Constitutive Relations), Section 5.2 (Properties)

### 3.1 Mixture Rules (`src/Properties/mixture_rules.jl`)

Implement mixture property computations:

```julia
function mixture_density(ξ::NTuple{NC,Float64}, T::Float64, material::Material) -> Float64
    # Eq: 1/ρ = Σ_solid(Yⱼ/ρⱼ) + Σ_liquid(Yⱼ/ρⱼ) + γ·Σ_gas(Yⱼ/ρⱼ)
end

function volume_fractions(ξ::NTuple{NC,Float64}, T::Float64, material::Material) -> NTuple{NC,Float64}
end

function mass_fractions(ξ::NTuple{NC,Float64}) -> NTuple{NC,Float64}
    # Yⱼ = ξⱼ / Σξⱼ
end
```

### 3.2 Effective Properties (`src/Properties/effective_properties.jl`)

```julia
function effective_heat_capacity(ξ, T, material) -> Float64
    # ρcₚ = Σⱼξⱼcₚ,ⱼ(T)
end

function effective_conductivity(ξ, T, material) -> Float64
    # k_eff = β·k_parallel + (1-β)·k_series
end

function effective_emissivity(ξ, material) -> Float64
    # ε_eff = Σⱼvⱼεⱼ
end

function effective_absorption(ξ, material) -> Float64
    # α_eff = Σⱼξⱼαⱼ
end

function effective_gas_transfer(ξ, T, material) -> Float64
end
```

### 3.3 Property Cache (`src/Properties/property_cache.jl`)

Pre-allocate storage for evaluated properties to avoid allocations:

```julia
mutable struct PropertyCache{NC}
    # Per-cell properties
    ρ::Vector{Float64}
    ρcp::Vector{Float64}
    k::Vector{Float64}
    λ::Vector{Float64}
    ε::Vector{Float64}
    α::Vector{Float64}

    # Per-face properties (for flux computation)
    k_face::Vector{Float64}
    λ_face::Vector{Float64}
end

function PropertyCache(mesh::Mesh1D, NC::Int) -> PropertyCache{NC}

function update_properties!(cache::PropertyCache, mesh::Mesh1D, material::Material)
    # Evaluate all properties for current state
end
```

### Tests for Stage 3 (`test/test_properties.jl`)

1. Test mixture density for known compositions
2. Test conductivity mixing (parallel, series, weighted) against hand calculations
3. Test property cache update is allocation-free
4. Verify mass/volume fraction consistency

### Acceptance Criteria

- [ ] Mixture properties match analytical formulas
- [ ] Property cache updates without allocations
- [ ] Type stability verified

---

## Stage 4: Physics Module (No ALE)

**Goal**: Implement physical models for pyrolysis on a fixed mesh

**Reference**: Design Doc Section 2 (Mathematical Framework), Section 7.2 (Physics Step)

### 4.1 Reaction Kinetics (`src/Physics/kinetics.jl`)

```julia
function reaction_rate(rxn::Reaction, ξ, T) -> Float64
    # rᵢ = Aᵢ·exp(-Eᵢ/RT)·ξᵣ₁^n₁·ξᵣ₂^n₂·H(T - Tᵢ,min)·H(Tᵢ,max - T)
end

function evaluate_reactions!(rates::Vector{Float64}, heats::Vector{Float64},
                            ξ, T, reactions::ReactionSet)
    # Evaluate all reaction rates and heats
end

function species_source_terms!(sources::Vector{Float64}, rates::Vector{Float64},
                               reactions::ReactionSet)
    # Compute dξⱼ/dt from reaction rates using stoichiometry
end
```

### 4.2 Heat Conduction (`src/Physics/heat_conduction.jl`)

```julia
function conductive_flux(k_L, k_R, T_L, T_R, d_LR) -> Float64
    # q = -k_f · (T_R - T_L) / d_LR
    # k_f = harmonic mean
end

function compute_conductive_fluxes!(fluxes::Vector{Float64},
                                    mesh::Mesh1D, cache::PropertyCache)
end
```

### 4.3 Radiation (`src/Physics/radiation.jl`)

Implement Beer-Lambert in-depth absorption:

```julia
function compute_radiation_source!(q_rad::Vector{Float64},
                                   mesh::Mesh1D, cache::PropertyCache,
                                   I_surface::Float64)
    # I(z) = I₀·exp(-∫₀ᶻ α_eff(z')dz')
    # q'''(z) = α_eff · I(z)
end
```

### 4.4 Gas Transport (`src/Physics/gas_transport.jl`)

```julia
function gas_flux(λ_L, λ_R, ξ_L, ξ_R, ρ_L, ρ_R, d_LR) -> Float64
    # Darcy-Fick: Jⱼ = -ρⱼλⱼ∇(ξⱼ/ρⱼ)
end

function compute_gas_fluxes!(J::Matrix{Float64}, mesh::Mesh1D,
                             cache::PropertyCache, material::Material)
    # J[j,f] = flux of component j at face f
end
```

### 4.5 Source Terms (`src/Physics/source_terms.jl`)

Combine all source terms:

```julia
function compute_energy_sources!(Q::Vector{Float64}, mesh::Mesh1D,
                                 rates::Vector{Float64}, heats::Vector{Float64})
    # Q_rxn = -Σᵢ hᵢ rᵢ
    # where h < 0 is exothermic (heating) and h > 0 is endothermic (cooling)
end
```

### 4.6 Volume Change (`src/Physics/volume_change.jl`)

```julia
function volume_change_rate(dξ_dt, material::Material) -> Float64
    # dV/dt = V·Σⱼ(1/ρⱼ)(dξⱼ/dt)
end
```

### Tests for Stage 4 (`test/test_physics.jl`)

1. **Kinetics**:
   - Isothermal first-order decay matches analytical: ξ(t) = ξ₀·exp(-kt)
   - Mass conservation in reactions with products

2. **Conduction**:
   - Steady-state linear temperature profile
   - Harmonic mean gives correct interface conductivity

3. **Radiation**:
   - Beer-Lambert exponential decay matches analytical
   - Energy conservation: ∫q'''dz = I₀(1 - exp(-αL))

4. **Gas Transport**:
   - Diffusion from point source gives Gaussian profile

### Acceptance Criteria

- [ ] All physics functions type-stable and allocation-free
- [ ] Analytical test cases pass
- [ ] Conservation properties verified

---

## Stage 5: Boundary Conditions

**Goal**: Implement all boundary condition types

**Reference**: Design Doc Section 2.5 (Boundary Conditions), Section 6.4 (BC Types)

### 5.1 Thermal BCs (`src/Core/boundary_conditions.jl`)

```julia
struct AdiabaticBC <: AbstractThermalBC end
struct HeatFluxBC{F} <: AbstractThermalBC
    flux::F  # q"(t) or constant
end
struct ConvectiveBC{Fh,FT} <: AbstractThermalBC
    h::Fh
    T_inf::FT
end
struct RadiativeBC{Fε,FT} <: AbstractThermalBC
    ε::Fε
    T_inf::FT
end
struct CombinedThermalBC{BCs<:Tuple} <: AbstractThermalBC
    bcs::BCs
end
struct FlameBC{Fq,Fc} <: AbstractThermalBC
    q_flame_max::Fq
    critical_mass_fluxes::Fc
    extinction_criterion::Float64
end
```

Implement `+` operator for combining BCs.

### 5.2 Mass BCs

```julia
struct ImpermeableBC <: AbstractMassBC end
struct MassFluxBC{F} <: AbstractMassBC
    flux::F
end
struct ConvectiveMassBC{Fh,Fξ} <: AbstractMassBC
    h_m::Fh
    ξ_inf::Fξ
end
```

### 5.3 Boundary Flux Evaluation (`src/Discretization/boundary_fluxes.jl`)

```julia
function thermal_bc_flux(bc::AdiabaticBC, T_surface, t) -> Float64
    return 0.0
end

function thermal_bc_flux(bc::HeatFluxBC, T_surface, t) -> Float64
    return bc.flux isa Number ? bc.flux : bc.flux(t)
end

function thermal_bc_flux(bc::ConvectiveBC, T_surface, t) -> Float64
    h = bc.h isa Number ? bc.h : bc.h(t)
    T_inf = bc.T_inf isa Number ? bc.T_inf : bc.T_inf(t)
    return h * (T_inf - T_surface)
end

function thermal_bc_flux(bc::RadiativeBC, T_surface, t) -> Float64
    ε = bc.ε isa Number ? bc.ε : bc.ε(T_surface)
    T_inf = bc.T_inf isa Number ? bc.T_inf : bc.T_inf(t)
    return ε * STEFAN_BOLTZMANN * (T_inf^4 - T_surface^4)
end

function thermal_bc_flux(bc::CombinedThermalBC, T_surface, t) -> Float64
    return sum(thermal_bc_flux(b, T_surface, t) for b in bc.bcs)
end
```

### Tests for Stage 5

1. Adiabatic BC gives zero flux
2. Heat flux BC returns correct time-varying flux
3. Convective BC with h=0 gives zero flux
4. Radiative BC matches Stefan-Boltzmann formula
5. Combined BC sums contributions correctly
6. Flame BC activates when critical mass flux exceeded

### Acceptance Criteria

- [ ] All BC types implemented
- [ ] BC combination with `+` works
- [ ] Flux evaluation is type-stable for each BC type

---

## Stage 6: Discretization Module (Fixed Mesh)

**Goal**: Assemble the spatial discretization for the governing equations

**Reference**: Design Doc Section 8.2 (Spatial Discretization), Section 7.2 (Physics RHS)

### 6.1 Finite Volume Operators (`src/Discretization/finite_volume.jl`)

```julia
function divergence_1d!(div::Vector{Float64}, fluxes::Vector{Float64},
                        mesh::Mesh1D)
    # div[i] = (1/V_i) * (F_{i+1/2} - F_{i-1/2}) * A
end
```

### 6.2 Gradient Reconstruction (`src/Discretization/gradient_reconstruction.jl`)

```julia
function compute_gradients!(∇T::Vector{Float64}, mesh::Mesh1D)
    # Cell-centered gradients via finite difference
end

function face_value_upwind(φ_L, φ_R, v) -> Float64
    return v > 0 ? φ_L : φ_R
end
```

### 6.3 Flux Limiters (`src/Discretization/flux_limiters.jl`)

```julia
@enum Limiter MINMOD=1 VANLEER=2 SUPERBEE=3

function apply_limiter(r::Float64, limiter::Limiter) -> Float64
end

function limited_face_value(φ_LL, φ_L, φ_R, limiter::Limiter) -> Float64
end
```

### 6.4 RHS Assembly (`src/Discretization/finite_volume.jl` continued)

```julia
function compute_rhs!(dT::Vector{Float64}, dξ::Matrix{Float64},
                      mesh::Mesh1D, cache::PropertyCache,
                      material::Material, bcs, t::Float64)
    # Full RHS evaluation for fixed mesh case
    # Returns dT/dt and dξⱼ/dt for all cells
end
```

### Tests for Stage 6

1. **MMS Verification** (Method of Manufactured Solutions):
   - Choose smooth T(z,t), compute required source
   - Verify error convergence rate (should be 2nd order)

2. **Steady conduction**: Linear profile with constant flux BCs

3. **Divergence operator**: Conservation check

### Acceptance Criteria

- [ ] RHS assembly is type-stable and low-allocation
- [ ] MMS shows 2nd-order convergence
- [ ] Conservation verified numerically

---

## Stage 7: Time Integration (Fixed Mesh)

**Goal**: Integrate in time using DifferentialEquations.jl

**Reference**: Design Doc Section 8.1 (Time Integration), Section 7.1 (Main Algorithm)

### 7.1 State Vector (`src/TimeIntegration/state_vector.jl`)

```julia
function state_to_vector!(u::Vector{Float64}, mesh::Mesh1D)
    # Pack T and ξ into flat vector for ODE solver
end

function vector_to_state!(mesh::Mesh1D, u::Vector{Float64})
    # Unpack flat vector into mesh cell states
end

function create_state_vector(mesh::Mesh1D, NC::Int) -> Vector{Float64}
```

### 7.2 ODE Function (`src/TimeIntegration/ode_function.jl`)

```julia
struct PyrolysisODEFunction{M,Mat,BC,PC}
    mesh::M
    material::Mat
    bcs::BC
    cache::PC
end

function (f::PyrolysisODEFunction)(du, u, p, t)
    vector_to_state!(f.mesh, u)
    update_properties!(f.cache, f.mesh, f.material)
    compute_rhs!(du, f.mesh, f.cache, f.material, f.bcs, t)
end
```

### 7.3 Timestep Control (`src/TimeIntegration/timestep_control.jl`)

```julia
function compute_timestep(mesh::Mesh1D, cache::PropertyCache,
                          cfl::Float64=0.5, safety::Float64=0.9) -> Float64
    # dt = min(dt_cfl, dt_diff) * safety
end
```

### 7.4 Jacobian (Optional) (`src/TimeIntegration/jacobian.jl`)

For implicit methods, either:
- Use ForwardDiff automatic differentiation
- Provide sparsity pattern for efficient finite differences

```julia
function jacobian_sparsity(n_cells::Int, NC::Int) -> SparseMatrixCSC
    # Return sparsity pattern (tridiagonal blocks)
end
```

### Tests for Stage 7

1. **Heat conduction verification**:
   - Semi-infinite solid with step BC
   - Compare to erf analytical solution
   - Verify 2nd-order convergence in space and time

2. **Single reaction verification**:
   - Isothermal decay: ξ(t) = ξ₀·exp(-kt)
   - Verify mass conservation

3. **Combined conduction + reaction**:
   - Qualitative behavior matches expectations
   - Conservation checks

### Acceptance Criteria

- [ ] ODE integration runs with DifferentialEquations.jl solvers
- [ ] Implicit solvers (Rodas5P, QNDF) work with Jacobian
- [ ] Analytical solutions matched to expected accuracy
- [ ] Mass and energy conservation verified

---

## Stage 8: Problem and Solver Interface

**Goal**: Create user-friendly problem specification and solve interface

**Reference**: Design Doc Section 6.7 (Problem Types), Section 12 (API Design)

### 8.1 Problem Type (`src/Problem/pyrolysis_problem.jl`)

```julia
struct PyrolysisProblem{M,Mat,TBC,MBC}
    mesh::M
    material::Mat
    thermal_bcs::Dict{Int,TBC}
    mass_bcs::Dict{Int,MBC}
    T_initial::Function
    ξ_initial::Vector{Function}
    external_flux::Function
    ambient_temperature::Function
    tspan::Tuple{Float64,Float64}
    dt_initial::Float64
    output_times::Vector{Float64}
end
```

### 8.2 Initial Conditions (`src/Problem/initial_conditions.jl`)

```julia
function apply_initial_conditions!(mesh::Mesh1D, problem::PyrolysisProblem)
    for cell in mesh.cells
        z = cell_center(mesh, cell.id)
        state = mesh.cell_states[cell.id]
        state.T = problem.T_initial(z)
        for j in 1:NC
            state.ξ = setindex(state.ξ, problem.ξ_initial[j](z), j)
        end
    end
end
```

### 8.3 Solve Interface (`src/Solvers/solve.jl`)

```julia
function solve(problem::PyrolysisProblem;
               solver=Rodas5P(),
               abstol=1e-8, reltol=1e-6,
               callbacks=nothing) -> PyrolysisSolution
end
```

### 8.4 Solution Type (`src/Solvers/solve.jl`)

```julia
mutable struct PyrolysisSolution
    problem::PyrolysisProblem
    times::Vector{Float64}
    temperatures::Vector{Vector{Float64}}
    concentrations::Vector{Vector{Vector{Float64}}}
    mass_loss_rate::Vector{Float64}
    surface_temperature::Vector{Float64}
    retcode::Symbol
    solve_time::Float64
end
```

### Tests for Stage 8

1. Create problem with simple PMMA-like material
2. Solve for short time, verify solution structure
3. Test different solvers (explicit, implicit)

### Acceptance Criteria

- [ ] High-level API matches design doc Section 12.1
- [ ] Problems can be created and solved in few lines
- [ ] Solution contains all expected outputs

---

## Stage 9: Verification Suite

**Goal**: Comprehensive verification against analytical solutions

**Reference**: Design Doc Section 10 (Verification & Validation)

### 9.1 Semi-Infinite Solid (`test/test_verification.jl`)

```julia
function test_semi_infinite_solid()
    # Setup: T(z,0) = T₀, T(0,t) = T_s
    # Solution: T = T_s + (T₀-T_s)·erf(z/(2√(αt)))
    # Verify 2nd-order convergence
end
```

### 9.2 Beer-Lambert Verification

```julia
function test_beer_lambert()
    # Verify I(z) = I₀·exp(-αz)
    # Verify energy conservation
end
```

### 9.3 Single Reaction Verification

```julia
function test_isothermal_decay()
    # ξ(t) = ξ₀·exp(-kt)
    # Verify exponential decay rate
end
```

### 9.4 MMS Verification

```julia
function test_mms_convergence()
    # Manufactured solution with computed source
    # Run at h, h/2, h/4, h/8
    # Verify convergence rate
end
```

### 9.5 Conservation Checks

```julia
function test_mass_conservation()
    # Total mass in domain + mass lost = initial mass
end

function test_energy_conservation()
    # With adiabatic boundaries, only Q_rxn changes total energy:
    # - exothermic (h < 0) increases energy
    # - endothermic (h > 0) decreases energy
end
```

### Acceptance Criteria

- [ ] All analytical solutions matched within tolerance
- [ ] 2nd-order convergence demonstrated
- [ ] Mass conservation to machine precision
- [ ] Energy conservation to expected tolerance

---

## Stage 10: ALE Module

**Goal**: Add mesh motion capability

**Reference**: Design Doc Section 3 (ALE Formulation), Section 6.5 (ALE Types), Section 7.3 (Mesh Motion)

### 10.1 Mesh Velocity Methods (`src/ALE/mesh_velocity.jl`)

```julia
abstract type AbstractMeshVelocityMethod end

struct LaplacianSmoothing <: AbstractMeshVelocityMethod
    iterations::Int
    relaxation::Float64
end

struct SpringAnalogy <: AbstractMeshVelocityMethod
    stiffness_power::Float64
end

struct WeightedALE <: AbstractMeshVelocityMethod
    α_lagrangian::Float64
    α_smoothing::Float64
end

function compute_mesh_velocity!(w::Vector{Float64}, mesh::Mesh1D,
                                state, method::AbstractMeshVelocityMethod)
end
```

### 10.2 Mesh Motion (`src/ALE/mesh_motion.jl`)

```julia
function move_mesh!(mesh::Mesh1D, w::Vector{Float64}, dt::Float64)
    for i in eachindex(mesh.nodes)
        mesh.nodes[i] = Node1D(
            mesh.nodes[i].id,
            mesh.nodes[i].z + dt * w[i],
            mesh.nodes[i].z_ref,
            mesh.nodes[i].is_boundary,
            mesh.nodes[i].boundary_id
        )
    end
    update_geometry!(mesh)
end
```

### 10.3 GCL Verification (`src/ALE/gcl.jl`)

```julia
function check_gcl(mesh::Mesh1D, w::Vector{Float64}, dt::Float64) -> Float64
    # Verify: V_new - V_old ≈ dt * Σ_faces(w·n·A)
    # Return maximum relative error
end
```

### 10.4 ALE Fluxes (`src/Discretization/ale_fluxes.jl`)

```julia
function compute_ale_convective_flux(φ_L, φ_R, v, w) -> Float64
    # F = (v - w) · φ, upwinded
end

function compute_rhs_ale!(dT, dξ, mesh, cache, material, bcs, w, t)
    # RHS with ALE convective terms
end
```

### 10.5 Remapping (`src/ALE/remapping.jl`)

```julia
struct FirstOrderRemap <: AbstractRemappingMethod end
struct SecondOrderRemap <: AbstractRemappingMethod
    limiter::Symbol
end

function remap_solution!(mesh_new::Mesh1D, mesh_old::Mesh1D,
                         u_old::Vector{Float64}, method::AbstractRemappingMethod)
    # Transfer solution from old mesh to new mesh
end
```

### Tests for Stage 10

1. **Uniform compression/expansion**:
   - Constant mesh velocity
   - Verify GCL holds exactly

2. **Surface regression**:
   - Moving top boundary
   - Verify conservation during motion

3. **Pure advection through moving mesh**:
   - Verify solution quality after remapping

### Acceptance Criteria

- [ ] Mesh can move without tangling
- [ ] GCL satisfied to machine precision
- [ ] Conservation maintained during mesh motion
- [ ] ALE terms correctly modify flux computation

---

## Stage 11: Adaptivity Module

**Goal**: Implement h-refinement and r-refinement

**Reference**: Design Doc Section 4 (Adaptive Mesh Refinement), Section 6.6 (Adaptivity Types)

### 11.1 Error Indicators (`src/Adaptivity/error_indicators.jl`)

```julia
struct GradientIndicator <: AbstractErrorIndicator
    variable::Symbol
    weight::Float64
end

struct CurvatureIndicator <: AbstractErrorIndicator
    variable::Symbol
    weight::Float64
end

struct ReactionRateIndicator <: AbstractErrorIndicator
    weight::Float64
end

function compute_indicator(ind::GradientIndicator, mesh, state) -> Vector{Float64}
    # η_i = h_i · |∇φ|_i
end
```

### 11.2 h-Refinement (`src/Adaptivity/h_refinement.jl`)

```julia
function refine_cell!(mesh::Mesh1D, cell_id::Int)
    # Split cell into two
    # Create new node at center
    # Update connectivity
    # Interpolate state to children
end

function perform_h_refinement!(mesh::Mesh1D, marked_cells::Vector{Int})
end
```

### 11.3 h-Coarsening (`src/Adaptivity/h_coarsening.jl`)

```julia
function coarsen_cells!(mesh::Mesh1D, cell_id_a::Int, cell_id_b::Int)
    # Merge two sibling cells
    # Remove intermediate node
    # Average state conservatively
end

function perform_h_coarsening!(mesh::Mesh1D, marked_pairs::Vector{Tuple{Int,Int}})
end
```

### 11.4 r-Refinement (`src/Adaptivity/r_refinement.jl`)

```julia
function r_refine!(mesh::Mesh1D, monitor::Function;
                   max_iter::Int=10, tol::Float64=1e-3)
    # Equidistribution algorithm
    # Move nodes to equalize monitor function integral
end
```

### 11.5 AMR Controller (`src/Adaptivity/amr_controller.jl`)

```julia
mutable struct AMRController
    error_indicator::AbstractErrorIndicator
    θ_refine::Float64
    θ_coarsen::Float64
    h_min::Float64
    h_max::Float64
    max_level::Int
    check_every::Int
end

function adapt!(mesh::Mesh1D, controller::AMRController, state)
    # Compute indicators
    # Mark cells
    # Perform h-refinement
    # Perform h-coarsening
    # Ensure 2:1 balance
end
```

### Tests for Stage 11

1. **Refinement conservation**: Mass/energy conserved when refining
2. **Coarsening conservation**: Mass/energy conserved when coarsening
3. **Indicator correctness**: Gradient indicator high where gradient high
4. **2:1 balance**: Mesh remains balanced after adaptivity

### Acceptance Criteria

- [ ] Cells can be split and merged
- [ ] Conservation maintained during adaptation
- [ ] Error indicators correctly identify regions needing refinement
- [ ] Mesh quality maintained

---

## Stage 12: Integration and Examples

**Goal**: Create working examples demonstrating full capability

**Reference**: Design Doc Section 12.1 (High-Level API), examples/ directory

### 12.1 Heat Conduction Example (`examples/01_heat_conduction.jl`)

Simple verification case:
- Constant properties
- Step heat flux BC
- Compare to analytical solution

### 12.2 Single Reaction Example (`examples/02_single_reaction.jl`)

PMMA-like gasification:
- Single component material
- First-order decomposition
- Plot MLR vs time

### 12.3 Multi-Step Kinetics (`examples/03_multi_step.jl`)

Wood-like material:
- Multiple components (cellulose, char, volatiles)
- Multiple reactions
- Demonstrate char layer formation

### 12.4 Charring with Ablation (`examples/04_charring_ablation.jl`)

Ablative material:
- Surface regression
- ALE mesh motion
- Track char/virgin interface

### 12.5 Intumescent Material (`examples/05_intumescent.jl`)

Swelling material:
- Material expansion during pyrolysis
- Large mesh deformation
- Demonstrate ALE capability

### 12.6 Adaptive Showcase (`examples/06_adaptive_showcase.jl`)

AMR demonstration:
- Sharp reaction front
- Dynamic refinement tracking front
- Compare to uniform mesh

### Acceptance Criteria

- [ ] All examples run without error
- [ ] Results match expected physical behavior
- [ ] Examples are documented and readable
- [ ] Performance is acceptable

---

## Stage 13: Post-Processing and Output

**Goal**: Add visualization and export capabilities

**Reference**: Design Doc Section 5.1 (PostProcess directory)

### 13.1 Standard Outputs (`src/PostProcess/outputs.jl`)

```julia
function mass_loss_rate(solution::PyrolysisSolution, t::Float64) -> Float64
function total_mass(solution::PyrolysisSolution, t::Float64) -> Float64
function surface_temperature(solution::PyrolysisSolution, t::Float64) -> Float64
function heat_release_rate(solution::PyrolysisSolution, t::Float64) -> Float64
```

### 13.2 Conservation Checks (`src/PostProcess/conservation_checks.jl`)

```julia
function mass_balance(solution::PyrolysisSolution) -> Vector{Float64}
function energy_balance(solution::PyrolysisSolution) -> Vector{Float64}
```

### 13.3 Visualization (`src/PostProcess/visualization.jl`)

Plotting recipes using RecipesBase:

```julia
@recipe function f(sol::PyrolysisSolution, ::Val{:mlr})
    # MLR vs time plot
end

@recipe function f(sol::PyrolysisSolution, ::Val{:profile}, t::Float64)
    # Temperature profile at time t
end
```

### 13.4 Export (`src/PostProcess/export.jl`)

```julia
function export_vtk(solution::PyrolysisSolution, filename::String)
function export_csv(solution::PyrolysisSolution, filename::String)
```

### Acceptance Criteria

- [ ] Standard outputs computable from solution
- [ ] Plotting recipes work with Plots.jl
- [ ] VTK export readable by ParaView
- [ ] CSV export contains all requested data

---

## Stage 14: Performance Optimization

**Goal**: Optimize critical paths for production use

**Reference**: Design Doc Section 11 (Performance Considerations), Section 9 (Implementation Guidelines)

### 14.1 Profiling

1. Profile full simulation with `@profile`
2. Identify hot spots
3. Target functions consuming >10% time

### 14.2 Type Stability Audit

For each hot function:
```julia
@code_warntype function_name(args...)
```

Fix any type instabilities.

### 14.3 Allocation Elimination

For each hot function:
```julia
@allocated function_name(args...)
```

Target: zero allocations in inner loops.

### 14.4 SIMD and Threading

1. Add `@simd` to appropriate loops
2. Add `Threads.@threads` for embarrassingly parallel operations
3. Benchmark to verify improvement

### 14.5 Benchmarking Suite

Create benchmark tests:
```julia
@btime solve($problem)
@btime compute_rhs!($du, $mesh, $cache, $material, $bcs, $t)
```

### Acceptance Criteria

- [ ] All hot paths type-stable
- [ ] Inner loops allocation-free
- [ ] Performance regression tests in place
- [ ] Documented baseline performance numbers

---

## Stage 15: Documentation and Testing Polish

**Goal**: Production-ready documentation and test coverage

### 15.1 API Documentation

Docstrings for all public functions per Section 9.7.

### 15.2 Tutorial Documentation (`docs/src/tutorial.md`)

Step-by-step guide to:
1. Creating materials
2. Setting up problems
3. Running simulations
4. Analyzing results

### 15.3 Theory Documentation (`docs/src/theory.md`)

Mathematical background:
1. Governing equations
2. ALE formulation
3. Numerical methods

### 15.4 Test Coverage

Target: >80% line coverage

Add tests for:
- Edge cases
- Error conditions
- All public API functions

### 15.5 Documenter.jl Setup

```julia
# docs/make.jl
using Documenter, Pyrolysis

makedocs(
    modules = [Pyrolysis],
    sitename = "Pyrolysis.jl",
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Theory" => "theory.md",
        "API" => "api.md"
    ]
)
```

### Acceptance Criteria

- [ ] All public functions documented
- [ ] Tutorial enables new users to get started
- [ ] Theory explains mathematical background
- [ ] Documentation builds without warnings
- [ ] Test coverage >80%

---

## Implementation Order Summary

| Stage | Module | Dependencies | Priority |
|-------|--------|--------------|----------|
| 0 | Scaffolding | None | Required |
| 1 | Core | None | Required |
| 2 | Geometry | Core | Required |
| 3 | Properties | Core | Required |
| 4 | Physics | Core, Properties | Required |
| 5 | BCs | Core | Required |
| 6 | Discretization | Geometry, Physics, BCs | Required |
| 7 | Time Integration | Discretization | Required |
| 8 | Problem/Solver | All above | Required |
| 9 | Verification | Problem/Solver | Required |
| 10 | ALE | Geometry, Discretization | High |
| 11 | Adaptivity | Geometry, ALE | Medium |
| 12 | Examples | All | High |
| 13 | PostProcess | Problem/Solver | Medium |
| 14 | Optimization | All | Medium |
| 15 | Documentation | All | High |

**Minimum Viable Product (MVP)**: Stages 0-9 provide a working fixed-mesh pyrolysis solver with verification.

**Full Product**: All 15 stages provide the complete ALE-adaptive solver as specified.

---

## Notes for Implementers

1. **Test constantly**: Run tests after each significant change
2. **Commit often**: Small, focused commits with clear messages
3. **Check types**: Use `@code_warntype` liberally during development
4. **Avoid premature optimization**: Get correctness first, then optimize
5. **Reference the design doc**: Section numbers are provided for each stage
6. **Ask questions**: Ambiguities should be resolved by consulting the design doc or domain experts

---

**Document Version**: 1.0
**Created**: Based on Pyrolysis_ALE_Framework_Design.md v2.0
