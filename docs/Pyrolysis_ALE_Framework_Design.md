# Pyrolysis.jl: ALE Adaptive Meshing Framework
## Comprehensive Design Document for Implementation

**Version**: 2.0 (ALE Formulation)
**Purpose**: This document provides complete specifications for implementing a Julia package that simulates axisymmetric 3D pyrolysis using an Arbitrary Lagrangian-Eulerian (ALE) adaptive meshing framework.

**Target**: AI-assisted implementation with human review

---

# Table of Contents

1. [Executive Summary](#1-executive-summary)
2. [Mathematical Framework](#2-mathematical-framework)
3. [ALE Formulation](#3-ale-formulation)
4. [Adaptive Mesh Refinement](#4-adaptive-mesh-refinement)
5. [Package Architecture](#5-package-architecture)
6. [Type System Design](#6-type-system-design)
7. [Algorithm Specifications](#7-algorithm-specifications)
8. [Numerical Methods](#8-numerical-methods)
9. [Implementation Guidelines](#9-implementation-guidelines)
10. [Verification & Validation](#10-verification--validation)
11. [Performance Considerations](#11-performance-considerations)
12. [API Design](#12-api-design)
13. [Dependencies & Compatibility](#13-dependencies--compatibility)

---

# 1. Executive Summary

## 1.1 Problem Statement

Simulate the pyrolysis (thermal decomposition) of materials under high heat flux conditions, accounting for:
- Multi-component materials with temperature-dependent properties
- Arrhenius-type decomposition kinetics (potentially multi-step)
- Gas generation and transport through porous char
- In-depth radiation absorption
- **Material volume changes (shrinkage, swelling, intumescence)**
- Surface regression and internal interface motion

## 1.2 Why ALE?

The Arbitrary Lagrangian-Eulerian framework offers critical advantages:

| Aspect | Fixed Grid (Eulerian) | ALE Adaptive |
|--------|----------------------|--------------|
| Shrinkage/Swelling | Approximate via porosity | Exact geometry tracking |
| Surface Regression | Requires immersed boundary | Natural boundary following |
| Interface Tracking | Diffuse over cells | Sharp interface preservation |
| Computational Cost | Uniform everywhere | Concentrated where needed |
| Conservation | Good | Excellent (with GCL) |
| Moving Boundaries | Complex | Natural |

## 1.3 Key Design Principles

1. **Geometric Conservation Law (GCL) Compliance**: Mesh motion must satisfy discrete GCL to preserve conservation
2. **Arbitrary Lagrangian-Eulerian**: Mesh velocity independent of material velocity
3. **h-Adaptivity**: Local refinement/coarsening based on error indicators
4. **r-Adaptivity**: Node relocation for smoothness and interface tracking
5. **Conservative Remapping**: Mass, energy, species conserved during mesh changes
6. **Type Stability**: All hot paths must be type-stable for JIT optimization
7. **Zero Allocations**: Inner loops must be allocation-free

---

# 2. Mathematical Framework

## 2.1 Physical Domain

**Geometry**: Axisymmetric domain Ω(t) with time-varying boundaries
- Radial coordinate: r ∈ [0, R_outer]
- Axial coordinate: z ∈ [0, L(t)] where L(t) is the current thickness
- For 1D (slab): degenerate to z only with unit area

**Coordinate Convention**:
- z = 0: Bottom boundary (substrate, fixed reference)
- z = L: Top boundary (exposed surface, may regress)
- Heat enters from top (z = L), material shrinks toward bottom (z = 0)
- Cell 1 is near z = 0 (bottom), Cell N is near z = L (top)
- boundary_id = 1: Top surface (exposed), boundary_id = 2: Bottom surface (substrate)

**Domain Boundaries**:
- Γ_top(t): Exposed surface at z = L(t) (may regress/advance)
- Γ_bottom: Substrate interface at z = 0 (fixed or prescribed)
- Γ_side: Adiabatic symmetry (for axisymmetric)

## 2.2 State Variables

Primary variables (per unit volume):
- T: Temperature [K]
- ξⱼ: Mass concentration of component j [kg/m³]
- ρ = Σξⱼ: Total density [kg/m³]

Derived quantities:
- Yⱼ = ξⱼ/ρ: Mass fraction of component j [-]
- vⱼ: Volume fraction of component j [-]
- γ: Porosity (gas volume fraction) [-]

## 2.3 Conservation Equations (Eulerian Form)

### Energy Conservation
```
∂(ρcₚT)/∂t + ∇·(ρcₚTv_g) = ∇·(k_eff∇T) - ∇·q_rad + Σᵢhᵢrᵢ
```

Where:
- ρcₚ = Σⱼξⱼcₚ,ⱼ(T): Volumetric heat capacity [J/(m³·K)]
- k_eff: Effective thermal conductivity [W/(m·K)]
- q_rad: Radiative heat flux [W/m²]
- hᵢ: Heat of reaction i [J/kg]
- rᵢ: Reaction rate i [kg/(m³·s)]
- v_g: Gas velocity [m/s]

### Species Conservation (Component j)
```
∂ξⱼ/∂t + ∇·(ξⱼvⱼ) = Σᵢνᵢⱼrᵢ
```

Where:
- vⱼ: Velocity of component j (v_s for solids, v_g for gases)
- νᵢⱼ: Stoichiometric coefficient (negative for reactants, positive for products)

### Mass Conservation (Overall)
```
∂ρ/∂t + ∇·(ρv) = 0
```

Note: In pyrolysis, solid mass converts to gas which escapes, so local mass is not conserved.

## 2.4 Constitutive Relations

### Arrhenius Kinetics
```
rᵢ = Aᵢ·exp(-Eᵢ/RT)·ξᵣ₁^n₁·ξᵣ₂^n₂·H(T - Tᵢ,min)·H(Tᵢ,max - T)
```

Where:
- Aᵢ: Pre-exponential factor [s⁻¹] or [m³/(kg·s)]
- Eᵢ: Activation energy [J/mol]
- R = 8.31446 J/(mol·K): Universal gas constant
- n₁, n₂: Reaction orders (typically 1)
- H: Heaviside function for temperature limits

### Gas Transport (Darcy-Fick)

The full Darcy-Fick formulation combines advective (pressure-driven) and diffusive transport:

```
Jⱼ = ξⱼ·v_g - Dⱼ,eff·∇ξⱼ     (combined advection + diffusion)
     ↑ advection  ↑ diffusion
```

**Fickian Diffusion** (always active):
```
J_diff = -λⱼ·ρⱼ·∇Yⱼ = -λⱼ·ρⱼ·∇(ξⱼ/ρⱼ)    for gaseous components
```

Where:
- λⱼ: Gas transfer coefficient [m²/s] (diffusion-like parameter)
- Yⱼ = ξⱼ/ρⱼ: Mass fraction of gas component j

**Darcy Flow** (optional, pressure-driven):
```
v_g = -(κ/μ)·∇P    (Darcy velocity)
```

Where:
- κ: Permeability [m²] (from permeability model)
- μ: Dynamic gas viscosity [Pa·s] (from Sutherland's law)
- P: Gas pressure [Pa] (from ideal gas law)

**Pressure Field** (ideal gas law):
```
P = Σⱼ (ξⱼ·R·T / Mⱼ)    for gaseous components
```

Where:
- Mⱼ: Molar mass of component j [kg/mol]
- R = 8.31446 J/(mol·K): Universal gas constant

**Permeability Models** (implemented):

1. *Constant Permeability*:
```
κ(φ,T) = κ₀
```

2. *Kozeny-Carman* (for packed beds/porous char):
```
κ(φ,T) = κ₀·(φ/φ₀)³·((1-φ₀)/(1-φ))²
```

Where φ is porosity (gas volume fraction) and φ₀ is reference porosity.

**Gas Viscosity** (Sutherland's law):
```
μ(T) = μ_ref·(T/T_ref)^(3/2)·(T_ref + S)/(T + S)
```

Where S ≈ 110 K for air-like pyrolysis gases, μ_ref = 1.81×10⁻⁵ Pa·s at T_ref = 300 K.

### Effective Thermal Conductivity
```
k_eff = β·k_parallel + (1-β)·k_series
```

Where:
- k_parallel = Σⱼvⱼkⱼ (parallel arrangement)
- k_series = 1/(Σⱼvⱼ/kⱼ) (series arrangement)
- β ∈ [0,1]: Mixing parameter (typically 0.5)

### Mixture Density with Swelling
```
1/ρ = Σ_solid(Yⱼ/ρⱼ) + Σ_liquid(Yⱼ/ρⱼ) + γ·Σ_gas(Yⱼ/ρⱼ)
```

Where γ ∈ [0,1] is the swelling factor:
- γ = 0: Rigid matrix (gas fills existing pores only)
- γ = 1: Full expansion (gas creates new volume)

### Volume Change Rate
```
∂V/∂t = V·Σⱼ(1/ρⱼ)(∂ξⱼ/∂t)
```

This drives the mesh motion in the Lagrangian limit.

## 2.5 Boundary Conditions

### Top Surface Γ_top(t)

**Energy**:
```
-k_eff(∂T/∂n) = q"_ext - ε_eff·σ(T⁴ - T∞⁴) - h_c(T - T∞) + q"_flame
```

Where:
- q"_ext: External imposed heat flux [W/m²]
- ε_eff = Σⱼvⱼεⱼ: Effective emissivity
- σ = 5.67×10⁻⁸ W/(m²·K⁴): Stefan-Boltzmann constant
- h_c: Convective heat transfer coefficient [W/(m²·K)]
- q"_flame: Heat flux from flame (if ignited) [W/m²]

**Flame Heat Flux Model**:
```
q"_flame = q"_flame,max · H(CI - 1)
CI = Σⱼ(Jⱼ/Jⱼ,crit)    (Critical mass flux ignition)
```

**Mass (Gaseous Components)**:
```
Jⱼ·n = (gas release rate based on interior gradient)
```

### Bottom Surface Γ_bottom

**Energy** (typically adiabatic or prescribed):
```
-k_eff(∂T/∂n) = 0    (adiabatic)
or
T = T_substrate(t)    (Dirichlet)
```

**Mass** (impermeable):
```
Jⱼ·n = 0    for all gaseous j
```

### Lateral Surface (Axisymmetric)

**At r = 0** (symmetry axis):
```
∂T/∂r = 0
∂ξⱼ/∂r = 0
```

## 2.6 Radiation Models

### Surface-Only Model
All radiation absorbed/emitted at surface. Simple but inaccurate for semi-transparent materials.

### Beer-Lambert In-Depth Absorption
```
I(z) = I₀·exp(-∫₀ᶻ α_eff(z')dz')
```

Where:
- α_eff = Σⱼξⱼαⱼ: Effective absorption coefficient [m⁻¹]
- αⱼ: Mass absorption coefficient of component j [m²/kg]

**Volumetric absorption rate**:
```
q"'_rad = -dI/dz = α_eff·I(z)
```

### Two-Flux Model (Optional, More Accurate)
For materials with significant scattering:
```
dI⁺/dz = -(α + σ_s)I⁺ + σ_s·I⁻
dI⁻/dz = +(α + σ_s)I⁻ - σ_s·I⁺
```

---

# 3. ALE Formulation

## 3.1 Reference Frames

Three reference frames are used:
1. **Eulerian (Spatial)**: Fixed in space, coordinates x
2. **Lagrangian (Material)**: Moves with material, coordinates X
3. **ALE (Referential)**: Moves with mesh, coordinates χ

**Mappings**:
- φ: χ → x (ALE to spatial)
- Ψ: X → x (material to spatial)
- Φ: X → χ (material to ALE)

## 3.2 ALE Time Derivative

The time derivative in the ALE frame relates to the spatial frame as:
```
∂φ/∂t|_χ = ∂φ/∂t|_x + w·∇φ
```

Where w = ∂x/∂t|_χ is the mesh velocity.

## 3.3 ALE Conservation Laws

### General Conservation Law
For a conserved quantity φ with flux F:
```
∂φ/∂t|_χ + ∇·((v - w)φ) + ∇·F = S
```

Or in integral form over a control volume V(t):
```
d/dt ∫_V φ dV + ∮_∂V φ(v - w)·n dA + ∮_∂V F·n dA = ∫_V S dV
```

### ALE Energy Equation
```
∂(ρcₚT)/∂t|_χ + ∇·((v_g - w)ρcₚT) = ∇·(k_eff∇T) - ∇·q_rad + Σᵢhᵢrᵢ
```

### ALE Species Equation
```
∂ξⱼ/∂t|_χ + ∇·((vⱼ - w)ξⱼ) = Σᵢνᵢⱼrᵢ
```

For solids: vⱼ = w (Lagrangian tracking)
For gases: vⱼ = v_g (transport relative to mesh)

## 3.4 Geometric Conservation Law (GCL)

The mesh motion must satisfy the GCL to maintain conservation:
```
d/dt ∫_V dV = ∮_∂V w·n dA
```

In discrete form for cell i:
```
dVᵢ/dt = Σ_faces(w·n·A)_f
```

**Critical**: The discrete GCL must be satisfied to machine precision, otherwise spurious sources/sinks appear.

## 3.5 Mesh Velocity Strategies

### Pure Lagrangian (w = v_material)
- Mesh follows material exactly
- Perfect for tracking interfaces
- Can lead to mesh tangling with large deformations
- Use for: char/virgin interface tracking

### Pure Eulerian (w = 0)
- Mesh fixed in space
- Material flows through mesh
- Use for: far-field regions with little interest

### Arbitrary (w chosen for quality)
- Mesh velocity chosen to optimize mesh quality
- Balance between Lagrangian accuracy and Eulerian stability
- Use for: most of the domain

### Surface-Following
At the top surface:
```
w_surface = dz_surface/dt · n
```

Where the surface position evolves based on mass loss rate.

## 3.6 Mesh Motion Algorithms

### Laplacian Smoothing
Solve for mesh velocity w that smooths node positions:
```
∇²w = 0    in Ω
w = w_bc   on ∂Ω
```

### Linear Elasticity Analogy
Treat mesh as elastic solid:
```
∇·σ_mesh = 0
σ_mesh = λ_mesh(∇·d)I + μ_mesh(∇d + ∇dᵀ)
```

Where d is node displacement. Stiffer near boundaries, softer in interior.

### Spring Analogy
Each edge acts as a spring:
```
wᵢ = Σⱼ kᵢⱼ(xⱼ - xᵢ) / Σⱼ kᵢⱼ
```

Spring stiffness inversely proportional to edge length: kᵢⱼ = 1/|xⱼ - xᵢ|

### Recommended: Weighted Combination
```
w = α_lag·w_material + α_smooth·w_laplacian + α_track·w_interface
```

With weights chosen to prioritize:
1. Interface tracking (reaction fronts)
2. Surface following (boundary regression)
3. Quality maintenance (prevent tangling)

---

# 4. Adaptive Mesh Refinement

## 4.1 Error Indicators

### Gradient-Based Indicator
```
ηᵢ = hᵢ·|∇φ|ᵢ
```

Where hᵢ is cell size and φ is the variable of interest (T or ξⱼ).

### Curvature-Based Indicator
```
ηᵢ = hᵢ²·|∇²φ|ᵢ
```

More sensitive to sharp fronts.

### Reaction Rate Indicator
```
ηᵢ = hᵢ·max_j(|rⱼ|ᵢ / max(ξⱼ,ᵢ, ε))
```

Refine where reactions are active.

### Combined Indicator
```
η = (w_T·η_T² + w_ξ·η_ξ² + w_r·η_r²)^(1/2)
```

With weights w_T, w_ξ, w_r chosen by user.

### Interface Proximity Indicator
```
ηᵢ = 1 / (1 + d_interface/h_target)
```

Force refinement near tracked interfaces.

## 4.2 Refinement Criteria

**Mark for Refinement if**:
```
ηᵢ > θ_refine · max(η)    AND    hᵢ > h_min
```

**Mark for Coarsening if**:
```
ηᵢ < θ_coarsen · max(η)    AND    hᵢ < h_max
```

Typical values: θ_refine = 0.7, θ_coarsen = 0.1

## 4.3 h-Refinement (Cell Splitting)

### 1D Refinement
Split cell i into two cells:
- New node at cell center
- Redistribute state variables conservatively

### 2D Refinement (Quadtree for structured, red-green for unstructured)
For axisymmetric structured grids, use anisotropic refinement:
- Refine in z-direction preferentially (thermal gradients)
- Refine in r-direction near outer boundary

### Conservation During Refinement
When cell i splits into cells a, b:
```
φ_a·V_a + φ_b·V_b = φᵢ·Vᵢ
```

With sub-cell interpolation for smooth fields.

## 4.4 h-Coarsening (Cell Merging)

Merge cells only if:
1. Both marked for coarsening
2. Were created by a previous refinement (maintain conformity)
3. Resulting cell meets quality criteria

### Conservation During Coarsening
When cells a, b merge into cell i:
```
φᵢ = (φ_a·V_a + φ_b·V_b) / (V_a + V_b)
```

## 4.5 r-Refinement (Node Relocation)

Move nodes to equidistribute error:
```
∫₀^xᵢ ω(x)dx = i/N · ∫₀^L ω(x)dx
```

Where ω is a monitor function:
```
ω = √(1 + α|∂φ/∂x|²)
```

### Iterative Relocation Algorithm
1. Compute monitor function ω on current mesh
2. Compute target positions x̂ᵢ from equidistribution
3. Move nodes: xᵢ^(n+1) = xᵢ^n + β(x̂ᵢ - xᵢ^n)
4. Remap solution to new mesh
5. Repeat until converged

## 4.6 Combined hr-Adaptivity

**Strategy**:
1. Perform h-adaptivity first (change topology)
2. Perform r-adaptivity second (smooth and optimize)
3. Apply mesh motion (ALE step)

**Frequency**:
- h-adaptivity: Every N_h timesteps (expensive)
- r-adaptivity: Every N_r timesteps (moderate)
- Mesh motion: Every timestep (cheap)

Typical: N_h = 50, N_r = 10

## 4.7 Mesh Quality Metrics

### 1D Quality
- Aspect ratio: Δzᵢ / Δzᵢ₊₁ (should be < 2)
- Smoothness: |Δzᵢ - Δzᵢ₊₁| / Δzᵢ (should be < 0.5)

### 2D Quality
- Jacobian: det(J) > 0 everywhere (non-negative volume)
- Aspect ratio: max(Δr,Δz) / min(Δr,Δz) < AR_max
- Skewness: angle deviation from 90° < θ_max

### Quality-Based Limiting
If mesh motion would violate quality:
1. Reduce mesh velocity magnitude
2. Apply smoothing iterations
3. If still violated, trigger h-adaptivity

---

# 5. Package Architecture

## 5.1 Directory Structure

```
Pyrolysis.jl/
├── Project.toml
├── README.md
├── docs/
│   ├── make.jl
│   └── src/
│       ├── index.md
│       ├── theory.md
│       ├── tutorial.md
│       └── api.md
├── src/
│   ├── Pyrolysis.jl                 # Main module
│   │
│   ├── Core/                        # Fundamental types
│   │   ├── Core.jl                  # Submodule
│   │   ├── constants.jl             # Physical constants
│   │   ├── property_functions.jl    # Temperature-dependent properties
│   │   ├── components.jl            # Material components
│   │   ├── reactions.jl             # Reaction definitions
│   │   ├── materials.jl             # Composite materials
│   │   └── boundary_conditions.jl   # BC types
│   │
│   ├── Geometry/                    # Mesh and geometry
│   │   ├── Geometry.jl              # Submodule
│   │   ├── nodes.jl                 # Node types
│   │   ├── cells.jl                 # Cell types  
│   │   ├── faces.jl                 # Face types
│   │   ├── mesh.jl                  # Mesh container
│   │   ├── mesh_generation.jl       # Initial mesh creation
│   │   ├── mesh_io.jl               # Read/write meshes
│   │   └── geometric_quantities.jl  # Areas, volumes, centroids
│   │
│   ├── ALE/                         # ALE-specific functionality
│   │   ├── ALE.jl                   # Submodule
│   │   ├── mesh_velocity.jl         # Mesh velocity computation
│   │   ├── mesh_motion.jl           # Node position updates
│   │   ├── gcl.jl                   # Geometric conservation law
│   │   ├── remapping.jl             # Solution transfer
│   │   └── interface_tracking.jl    # Sharp interface methods
│   │
│   ├── Adaptivity/                  # Adaptive mesh refinement
│   │   ├── Adaptivity.jl            # Submodule
│   │   ├── error_indicators.jl      # Refinement criteria
│   │   ├── h_refinement.jl          # Cell splitting
│   │   ├── h_coarsening.jl          # Cell merging
│   │   ├── r_refinement.jl          # Node relocation
│   │   ├── mesh_quality.jl          # Quality metrics
│   │   └── amr_controller.jl        # Orchestrates adaptivity
│   │
│   ├── Properties/                  # Material property evaluation
│   │   ├── Properties.jl            # Submodule
│   │   ├── mixture_rules.jl         # Density, volume fractions
│   │   ├── effective_properties.jl  # Thermal, transport, gas viscosity
│   │   ├── permeability.jl          # Permeability models (κ)
│   │   └── property_cache.jl        # Cached evaluations (incl. P, φ, κ, v_gas)
│   │
│   ├── Physics/                     # Physical models
│   │   ├── Physics.jl               # Submodule
│   │   ├── kinetics.jl              # Reaction rates
│   │   ├── heat_conduction.jl       # Fourier's law
│   │   ├── radiation.jl             # Beer-Lambert, surface
│   │   ├── gas_transport.jl         # Darcy-Fick gas fluxes
│   │   ├── pressure.jl              # Ideal gas pressure field
│   │   ├── convection.jl            # Gas convective heat
│   │   ├── source_terms.jl          # Combined sources
│   │   └── volume_change.jl         # Swelling/shrinkage rates
│   │
│   ├── Discretization/              # Spatial discretization
│   │   ├── Discretization.jl        # Submodule
│   │   ├── finite_volume.jl         # FVM operators, compute_rhs_darcy_fick!
│   │   ├── darcy_velocity.jl        # Darcy velocity computation
│   │   ├── flux_limiters.jl         # TVD limiters
│   │   ├── gradient_reconstruction.jl # Face gradients
│   │   ├── ale_fluxes.jl            # ALE convective fluxes
│   │   └── boundary_fluxes.jl       # BC implementation
│   │
│   ├── TimeIntegration/             # Temporal discretization
│   │   ├── TimeIntegration.jl       # Submodule
│   │   ├── state_vector.jl          # DOF management
│   │   ├── ode_function.jl          # RHS assembly
│   │   ├── ale_timestep.jl          # Combined ALE+physics step
│   │   ├── jacobian.jl              # Analytical Jacobian
│   │   └── timestep_control.jl      # Adaptive dt
│   │
│   ├── Problem/                     # Problem setup
│   │   ├── Problem.jl               # Submodule
│   │   ├── pyrolysis_problem.jl     # Main problem type
│   │   ├── initial_conditions.jl    # IC specification
│   │   └── problem_validation.jl    # Input checking
│   │
│   ├── Solvers/                     # Solution algorithms
│   │   ├── Solvers.jl               # Submodule
│   │   ├── solve.jl                 # Main interface
│   │   ├── operator_splitting.jl    # Physics/mesh splitting
│   │   ├── callbacks.jl             # Events, output
│   │   └── steady_state.jl          # Equilibrium solver
│   │
│   └── PostProcess/                 # Output and analysis
│       ├── PostProcess.jl           # Submodule
│       ├── outputs.jl               # Standard outputs (MLR, T)
│       ├── conservation_checks.jl   # Mass, energy balance
│       ├── visualization.jl         # Plotting recipes
│       └── export.jl                # VTK, CSV output
│
├── test/
│   ├── runtests.jl
│   ├── test_core.jl
│   ├── test_mesh.jl
│   ├── test_ale.jl
│   ├── test_adaptivity.jl
│   ├── test_physics.jl
│   ├── test_discretization.jl
│   ├── test_integration.jl
│   ├── test_darcy_infrastructure.jl # Darcy-Fick unit tests
│   └── test_verification.jl         # MMS, analytical solutions, Darcy verification
│
└── examples/
    ├── 01_heat_conduction.jl        # Simple verification
    ├── 02_single_reaction.jl        # Basic pyrolysis
    ├── 03_multi_step.jl             # Complex kinetics
    ├── 04_charring_ablation.jl      # With surface regression
    ├── 05_intumescent.jl            # Swelling material
    ├── 06_cone_calorimeter.jl       # Standard test
    ├── 07_adaptive_showcase.jl      # AMR demonstration
    └── 08_darcy_fick_demo.jl        # Darcy-Fick gas transport
```

## 5.2 Module Dependencies

```
Pyrolysis (main)
├── Core (no deps)
├── Geometry → Core
├── Properties → Core
├── Physics → Core, Properties
├── ALE → Geometry, Properties
├── Adaptivity → Geometry, ALE, Physics
├── Discretization → Geometry, Physics, ALE
├── TimeIntegration → Discretization, Properties, Physics
├── Problem → all above
├── Solvers → Problem, TimeIntegration, Adaptivity
└── PostProcess → Geometry, Properties
```

## 5.3 Submodule Pattern

Each subdirectory should have a submodule file that:
1. Includes all files in the directory
2. Exports public symbols
3. Documents the submodule purpose

Example `src/ALE/ALE.jl`:
```julia
module ALE

using ..Core
using ..Geometry
using ..Properties
using LinearAlgebra
using SparseArrays

include("mesh_velocity.jl")
include("mesh_motion.jl")
include("gcl.jl")
include("remapping.jl")
include("interface_tracking.jl")

export compute_mesh_velocity!, move_mesh!, check_gcl
export remap_solution!, InteraceTracker
export MeshVelocityMethod, Laplacian, SpringAnalogy, LinearElasticity

end # module
```

---

# 6. Type System Design

## 6.1 Abstract Type Hierarchy

```julia
# Property functions
abstract type AbstractPropertyFunction end
abstract type AbstractThermalProperty <: AbstractPropertyFunction end
abstract type AbstractTransportProperty <: AbstractPropertyFunction end

# Material system
abstract type AbstractComponent end
abstract type AbstractReaction end
abstract type AbstractMaterial end

# Geometry
abstract type AbstractNode end
abstract type AbstractCell end
abstract type AbstractFace end
abstract type AbstractMesh{DIM} end

# Boundary conditions
abstract type AbstractBoundaryCondition end
abstract type AbstractThermalBC <: AbstractBoundaryCondition end
abstract type AbstractMassBC <: AbstractBoundaryCondition end
abstract type AbstractMeshBC <: AbstractBoundaryCondition end

# ALE
abstract type AbstractMeshVelocityMethod end
abstract type AbstractRemappingMethod end
abstract type AbstractInterfaceTracker end

# Adaptivity
abstract type AbstractErrorIndicator end
abstract type AbstractRefinementCriterion end
abstract type AbstractQualityMetric end

# Problem and solution
abstract type AbstractPyrolysisProblem end
abstract type AbstractPyrolysisSolution end
```

## 6.2 Core Types

### Physical Constants
```julia
module Constants
    const R_GAS = 8.31446261815324           # J/(mol·K)
    const STEFAN_BOLTZMANN = 5.670374419e-8  # W/(m²·K⁴)
    const T_REF = 298.15                     # K (reference temperature)
    const P_REF = 101325.0                   # Pa (reference pressure)
end
```

### Property Functions
```julia
struct ConstantProperty{T<:Real} <: AbstractPropertyFunction
    value::T
end

struct LinearProperty{T<:Real} <: AbstractPropertyFunction
    a::T  # intercept
    b::T  # slope
end

struct PolynomialProperty{N,T<:Real} <: AbstractPropertyFunction
    coeffs::NTuple{N,T}  # (c₀, c₁, ..., cₙ₋₁) for c₀ + c₁T + ...
end

struct ArrheniusProperty{T<:Real} <: AbstractPropertyFunction
    A::T       # pre-exponential
    E::T       # activation energy [J/mol]
    T_ref::T   # reference temperature
end

struct TableProperty{N,T<:Real} <: AbstractPropertyFunction
    temperatures::NTuple{N,T}
    values::NTuple{N,T}
end

struct CallableProperty{F} <: AbstractPropertyFunction
    func::F  # Any callable f(T) -> value
end

# Evaluation interface
(p::ConstantProperty)(T) = p.value
(p::LinearProperty)(T) = p.a + p.b * T
(p::PolynomialProperty)(T) = evalpoly(T, p.coeffs)
(p::ArrheniusProperty)(T) = p.A * exp(-p.E / (R_GAS * T))
(p::TableProperty)(T) = linear_interpolate(p.temperatures, p.values, T)
(p::CallableProperty)(T) = p.func(T)
```

### Components
```julia
@enum Phase SOLID=1 LIQUID=2 GAS=3

struct Component{
    Pρ<:AbstractPropertyFunction,
    Pc<:AbstractPropertyFunction,
    Pk<:AbstractPropertyFunction,
    Pλ<:AbstractPropertyFunction,
    Pκ<:Union{Nothing,AbstractPropertyFunction}
} <: AbstractComponent
    name::Symbol
    phase::Phase
    molar_mass::Float64           # kg/mol (for gases)
    
    # Thermophysical properties (all T-dependent)
    density::Pρ                   # kg/m³
    heat_capacity::Pc             # J/(kg·K)
    thermal_conductivity::Pk      # W/(m·K)
    gas_transfer_coeff::Pλ        # m²/s
    permeability::Pκ              # m² (Nothing for gases)
    
    # Radiative properties
    emissivity::Float64           # [-]
    absorption_coeff::Float64     # m²/kg
    
    # Volume change
    swelling_factor::Float64      # 0-1
end

# Convenience constructors
function SolidComponent(name; ρ, c, k, ε=0.9, α=0.0, γ=1.0, κ=0.0)
    Component(name, SOLID, 0.0, ρ, c, k, ConstantProperty(0.0),
              ConstantProperty(κ), ε, α, γ)
end

function GaseousComponent(name; M, c, k, λ, ε=0.0, α=0.0)
    # For ideal gas: ρ = PM/(RT)
    ρ_func = CallableProperty(T -> P_REF * M / (R_GAS * T))
    Component(name, GAS, M, ρ_func, c, k, λ, nothing, ε, α, 1.0)
end
```

### Reactions
```julia
struct Reaction{Ph<:AbstractPropertyFunction} <: AbstractReaction
    name::Symbol
    
    # Reactants (indices into component array, 0 = none)
    reactant1::Int
    reactant2::Int
    order1::Float64      # reaction order for reactant 1
    order2::Float64      # reaction order for reactant 2
    
    # Products (indices, 0 = none)  
    product1::Int
    product2::Int
    
    # Stoichiometry: mass of product per mass of reactant consumed
    # ν₁ kg reactant1 + ν₂ kg reactant2 → ν₃ kg product1 + ν₄ kg product2
    stoich::NTuple{4,Float64}
    
    # Kinetic parameters
    A::Float64           # pre-exponential [s⁻¹] or [m³/(kg·s)]
    E::Float64           # activation energy [J/mol]
    
    # Heat of reaction (T-dependent)
    heat::Ph             # [J/kg of reactant1 consumed]
    
    # Temperature limits
    T_min::Float64       # K
    T_max::Float64       # K
end

struct ReactionSet{NR,Ph} <: AbstractReaction
    reactions::NTuple{NR,Reaction{Ph}}
end
```

### Materials
```julia
@enum MixingRule PARALLEL=1 SERIES=2 WEIGHTED=3

struct Material{NC,NR,C<:AbstractComponent,R<:AbstractReaction} <: AbstractMaterial
    name::Symbol
    components::NTuple{NC,C}
    reactions::NTuple{NR,R}
    
    # Mixing parameters
    conductivity_mixing::MixingRule
    β::Float64                    # weighting for WEIGHTED rule
    gas_transfer_mixing::MixingRule
    β_λ::Float64
    permeability_mixing::MixingRule  # How to mix component permeabilities
    β_κ::Float64                     # weighting for permeability mixing
end

# Accessor functions
n_components(m::Material{NC}) where NC = NC
n_reactions(m::Material{NC,NR}) where {NC,NR} = NR
n_solids(m::Material) = count(c -> c.phase == SOLID, m.components)
n_gases(m::Material) = count(c -> c.phase == GAS, m.components)
has_darcy_flow(m::Material) = any(c -> c.phase != GAS && c.permeability !== nothing && 
                                       c.permeability(300.0) > 0, m.components)
```

## 6.3 Geometry Types

### Nodes
```julia
struct Node1D <: AbstractNode
    id::Int
    z::Float64              # position (z=0 at bottom, z=L at top)
    z_ref::Float64          # reference (initial) position
    is_boundary::Bool
    boundary_id::Int        # 0 = interior, 1 = top (z=L), 2 = bottom (z=0)
end

struct Node2D <: AbstractNode
    id::Int
    r::Float64
    z::Float64
    r_ref::Float64
    z_ref::Float64
    is_boundary::Bool
    boundary_id::Int        # 0 = interior, 1 = top, 2 = bottom, 3 = axis, 4 = outer
end
```

### Cells
```julia
struct Cell1D <: AbstractCell
    id::Int
    node_ids::Tuple{Int,Int}     # (left, right) node indices
    face_ids::Tuple{Int,Int}     # (left, right) face indices
    level::Int                   # refinement level (0 = coarsest)
    parent::Int                  # parent cell (0 = original)
    children::Tuple{Int,Int}     # child cells (0,0 = leaf)
end

struct Cell2D <: AbstractCell
    id::Int
    node_ids::NTuple{4,Int}      # counter-clockwise from bottom-left
    face_ids::NTuple{4,Int}      # (bottom, right, top, left)
    level::Int
    parent::Int
    children::NTuple{4,Int}
end

# Cell state (separate for cache efficiency)
mutable struct CellState{NC}
    T::Float64                   # temperature
    ξ::NTuple{NC,Float64}        # mass concentrations
    volume::Float64              # current volume
    volume_ref::Float64          # reference volume
end
```

### Faces
```julia
struct Face1D <: AbstractFace
    id::Int
    node_id::Int                 # single node in 1D
    cell_ids::Tuple{Int,Int}     # (left, right), 0 = boundary
                                 # left_id=0 → bottom boundary (z=0)
                                 # right_id=0 → top boundary (z=L)
    is_boundary::Bool
    boundary_id::Int             # 0 = interior, 1 = top (z=L), 2 = bottom (z=0)
end

struct Face2D <: AbstractFace
    id::Int
    node_ids::Tuple{Int,Int}     # two nodes defining face
    cell_ids::Tuple{Int,Int}     # (owner, neighbor), 0 = boundary
    normal::Tuple{Float64,Float64}  # unit outward normal (from owner)
    is_boundary::Bool
    boundary_id::Int
end
```

### Mesh Container
```julia
mutable struct Mesh1D <: AbstractMesh{1}
    # Topology
    nodes::Vector{Node1D}
    cells::Vector{Cell1D}
    faces::Vector{Face1D}
    
    # State
    cell_states::Vector{CellState}
    
    # Geometry cache (recomputed after mesh motion)
    cell_centers::Vector{Float64}
    cell_volumes::Vector{Float64}
    face_positions::Vector{Float64}
    face_areas::Vector{Float64}      # = cross-sectional area
    
    # Connectivity helpers
    boundary_faces::Dict{Int,Vector{Int}}  # boundary_id => face indices
    
    # Refinement tracking
    max_level::Int
    n_active_cells::Int              # excludes refined-away cells
    active_cell_mask::BitVector
    
    # ALE state
    node_velocities::Vector{Float64}
    mesh_time::Float64               # time of current mesh configuration
end

mutable struct Mesh2D <: AbstractMesh{2}
    # Similar structure but with 2D types
    nodes::Vector{Node2D}
    cells::Vector{Cell2D}
    faces::Vector{Face2D}
    cell_states::Vector{CellState}
    
    # Geometry cache
    cell_centers::Vector{Tuple{Float64,Float64}}
    cell_volumes::Vector{Float64}
    face_centers::Vector{Tuple{Float64,Float64}}
    face_areas::Vector{Float64}
    face_normals::Vector{Tuple{Float64,Float64}}
    
    # Axisymmetric factors
    r_factors::Vector{Float64}       # 2πr for volume integration
    
    # Connectivity
    boundary_faces::Dict{Int,Vector{Int}}
    
    # Refinement
    max_level::Int
    n_active_cells::Int
    active_cell_mask::BitVector
    
    # ALE
    node_velocities::Vector{Tuple{Float64,Float64}}
    mesh_time::Float64
end
```

## 6.4 Boundary Condition Types

### Thermal BCs
```julia
struct AdiabaticBC <: AbstractThermalBC end

struct HeatFluxBC{F} <: AbstractThermalBC
    flux::F  # q"(t) or constant, positive = into domain
end

struct ConvectiveBC{Fh,FT} <: AbstractThermalBC
    h::Fh           # h(t) or constant
    T_inf::FT       # T∞(t) or constant
end

struct RadiativeBC{Fε,FT} <: AbstractThermalBC
    ε::Fε           # ε(T) or constant
    T_inf::FT       # T∞(t) or constant
end

struct CombinedThermalBC{BCs<:Tuple} <: AbstractThermalBC
    bcs::BCs
end

# Combine BCs with + operator
Base.:+(bc1::AbstractThermalBC, bc2::AbstractThermalBC) = 
    CombinedThermalBC((bc1, bc2))
Base.:+(cbc::CombinedThermalBC, bc::AbstractThermalBC) = 
    CombinedThermalBC((cbc.bcs..., bc))

struct FlameBC{Fq,Fc} <: AbstractThermalBC
    q_flame_max::Fq              # maximum flame heat flux
    critical_mass_fluxes::Fc     # Dict: component => J_crit
    extinction_criterion::Float64 # CI < this => extinction
end
```

### Mass BCs
```julia
struct ImpermeableBC <: AbstractMassBC end

struct DiffusiveMassBC <: AbstractMassBC
    h_m::Float64    # mass transfer coefficient [m/s]
end

struct MassFluxBC{F} <: AbstractMassBC
    flux::F  # J(t) for each gas component
end

struct ConvectiveMassBC{Fh,Fξ} <: AbstractMassBC
    h_m::Fh         # mass transfer coefficient
    ξ_inf::Fξ       # ambient concentration
end
```

**Boundary Condition Effects on Darcy Velocity**:

- `ImpermeableBC`: Zero flux at boundary (v_gas = 0)
- `DiffusiveMassBC`: Pressure-driven outflow, v_gas = -(κ/μ)·(P_atm - P_int)/d
- Surface boundary (id=1): Allows gas escape based on internal pressure gradient
- Back boundary (id=2): Typically impermeable (no gas escape through substrate)

**Note**: Permeability is specified per-component via the `κ` parameter in `SolidComponent`.
The effective permeability is computed from component values using mixing rules defined
in the `Material` struct (`permeability_mixing`, `β_κ`). Darcy flow is automatically
enabled when any solid/liquid component has κ > 0.

### Mesh BCs (for ALE)
```julia
struct FixedMeshBC <: AbstractMeshBC end

struct FreeSurfaceBC <: AbstractMeshBC
    # Surface follows material
end

struct PrescribedVelocityBC{F} <: AbstractMeshBC
    velocity::F  # w(t)
end

struct SlidingBC <: AbstractMeshBC
    # Node can move tangentially, not normally
    direction::Symbol  # :r or :z
end
```

## 6.5 ALE Types

### Mesh Velocity Methods
```julia
struct LaplacianSmoothing <: AbstractMeshVelocityMethod
    iterations::Int
    relaxation::Float64
end

struct SpringAnalogy <: AbstractMeshVelocityMethod
    stiffness_power::Float64  # k ∝ L^(-power)
end

struct LinearElasticity <: AbstractMeshVelocityMethod
    E_base::Float64           # base Young's modulus
    ν::Float64                # Poisson's ratio
    stiffening_near_boundary::Float64
end

struct WeightedALE <: AbstractMeshVelocityMethod
    α_lagrangian::Float64     # weight for material following
    α_smoothing::Float64      # weight for quality
    α_interface::Float64      # weight for interface tracking
    inner_method::AbstractMeshVelocityMethod
end
```

### Remapping Methods
```julia
struct FirstOrderRemap <: AbstractRemappingMethod end

struct SecondOrderRemap <: AbstractRemappingMethod
    limiter::Symbol  # :minmod, :vanleer, :superbee
end

struct ConservativeRemap <: AbstractRemappingMethod
    # Intersection-based, exactly conservative
end
```

### Interface Tracking
```julia
mutable struct InterfaceTracker <: AbstractInterfaceTracker
    # Track char/virgin interface
    interface_positions::Vector{Float64}
    interface_temperatures::Vector{Float64}
    tracking_variable::Symbol      # :T or component name
    threshold::Float64             # value defining interface
    
    # For level-set method
    use_level_set::Bool
    φ::Vector{Float64}             # level set field
end
```

## 6.6 Adaptivity Types

### Error Indicators
```julia
struct GradientIndicator <: AbstractErrorIndicator
    variable::Symbol       # :T or component name
    weight::Float64
end

struct CurvatureIndicator <: AbstractErrorIndicator
    variable::Symbol
    weight::Float64
end

struct ReactionRateIndicator <: AbstractErrorIndicator
    weight::Float64
end

struct CompositeIndicator <: AbstractErrorIndicator
    indicators::Vector{AbstractErrorIndicator}
    combination::Symbol   # :max, :l2, :sum
end
```

### Quality Metrics
```julia
struct AspectRatioMetric <: AbstractQualityMetric
    max_ratio::Float64
end

struct SmoothnessMetric <: AbstractQualityMetric
    max_size_jump::Float64
end

struct JacobianMetric <: AbstractQualityMetric
    min_jacobian::Float64
end

struct CompositeQualityMetric <: AbstractQualityMetric
    metrics::Vector{AbstractQualityMetric}
end
```

### AMR Controller
```julia
mutable struct AMRController{EI<:AbstractErrorIndicator, QM<:AbstractQualityMetric}
    error_indicator::EI
    quality_metric::QM
    
    # Refinement parameters
    θ_refine::Float64           # refine if η > θ_refine * max(η)
    θ_coarsen::Float64          # coarsen if η < θ_coarsen * max(η)
    h_min::Float64              # minimum cell size
    h_max::Float64              # maximum cell size
    max_level::Int              # maximum refinement levels
    
    # Frequency
    check_every::Int            # timesteps between AMR checks
    
    # Statistics
    n_refinements::Int
    n_coarsenings::Int
    last_check_time::Float64
end
```

## 6.7 Problem and Solution Types

### Pyrolysis Problem
```julia
struct PyrolysisProblem{
    M<:AbstractMesh,
    Mat<:AbstractMaterial,
    TBC,
    MBC,
    MeshBC
} <: AbstractPyrolysisProblem
    # Domain
    mesh::M
    material::Mat
    
    # Boundary conditions
    thermal_bcs::Dict{Int,TBC}      # boundary_id => BC
    mass_bcs::Dict{Int,MBC}
    mesh_bcs::Dict{Int,MeshBC}
    
    # Initial conditions
    T_initial::Function             # T(z) or T(r,z)
    ξ_initial::Vector{Function}     # ξⱼ(z) for each component
    
    # External forcing
    external_flux::Function         # q"(t)
    ambient_temperature::Function   # T∞(t)
    
    # Solver settings
    tspan::Tuple{Float64,Float64}
    dt_initial::Float64
    
    # ALE settings
    mesh_velocity_method::AbstractMeshVelocityMethod
    remapping_method::AbstractRemappingMethod
    interface_tracker::Union{Nothing,AbstractInterfaceTracker}
    
    # AMR settings
    amr_controller::Union{Nothing,AMRController}
    
    # Output settings
    output_times::Vector{Float64}
    save_mesh_history::Bool
end
```

### Solution Container
```julia
mutable struct PyrolysisSolution{M<:AbstractMesh} <: AbstractPyrolysisSolution
    # Reference to problem
    problem::PyrolysisProblem
    
    # Time history
    times::Vector{Float64}
    
    # Field history (only at output_times)
    temperatures::Vector{Vector{Float64}}
    concentrations::Vector{Vector{Vector{Float64}}}
    
    # Mesh history (if save_mesh_history)
    meshes::Vector{M}
    
    # Integral outputs (at all times)
    time_fine::Vector{Float64}
    mass_loss_rate::Vector{Float64}
    total_mass::Vector{Float64}
    surface_temperature::Vector{Float64}
    back_face_temperature::Vector{Float64}
    heat_release_rate::Vector{Float64}
    
    # Interface positions
    interface_history::Vector{Vector{Float64}}
    
    # Diagnostics
    dt_history::Vector{Float64}
    n_cells_history::Vector{Int}
    conservation_error::Vector{Float64}
    
    # Status
    retcode::Symbol               # :success, :dt_too_small, :max_iter, etc.
    solve_time::Float64           # wall-clock time
end
```

---

# 7. Algorithm Specifications

## 7.1 Main Solution Algorithm

```
ALGORITHM: Solve Pyrolysis Problem with ALE-AMR

INPUT: PyrolysisProblem P
OUTPUT: PyrolysisSolution

1. INITIALIZE
   a. Validate problem inputs
   b. Create initial mesh from P.mesh_spec
   c. Apply initial conditions: T(z), ξⱼ(z)
   d. Compute initial properties
   e. Initialize ALE state (w = 0)
   f. Initialize AMR controller
   g. Initialize solution storage

2. TIME LOOP: for t = t₀ to t_final
   
   2.1 COMPUTE TIMESTEP
       a. CFL condition for convection
       b. Diffusion stability limit
       c. Reaction timescale (stiffness indicator)
       d. Mesh motion stability
       e. dt = min(dt_cfl, dt_diff, dt_react, dt_mesh, dt_max)
   
   2.2 OPERATOR SPLITTING STEP
       
       Option A: Strang Splitting (second-order)
       a. Half-step mesh motion: mesh(t + dt/2) 
       b. Full-step physics: u(t + dt)
       c. Half-step mesh motion: mesh(t + dt)
       
       Option B: Lie Splitting (first-order, simpler)
       a. Full-step physics: u(t + dt) on fixed mesh
       b. Full-step mesh motion + remap
   
   2.3 PHYSICS STEP (on current mesh)
       a. Evaluate material properties at current state
       b. Compute reaction rates
       c. Compute conductive fluxes
       d. Compute radiative source (Beer-Lambert)
       e. Compute gas transport fluxes
       f. Compute ALE convective fluxes
       g. Assemble RHS: du/dt = f(u, t)
       h. Time integrate (implicit or IMEX)
       i. Update state: u^(n+1)
   
   2.4 MESH MOTION STEP
       a. Compute volume change rates from mass changes
       b. Compute surface regression rate
       c. Compute mesh velocity field w
       d. Verify GCL satisfaction
       e. Move nodes: x^(n+1) = x^n + dt·w
       f. Update geometric quantities
       g. Remap solution to new mesh (if needed)
   
   2.5 ADAPTIVITY STEP (every N_amr steps)
       a. Compute error indicators
       b. Mark cells for refinement/coarsening
       c. Check mesh quality metrics
       d. Perform h-refinement
       e. Perform h-coarsening
       f. Perform r-refinement (node smoothing)
       g. Conservative remap to new mesh
       h. Update connectivity structures
   
   2.6 OUTPUT AND DIAGNOSTICS
       a. Compute integral quantities (MLR, HRR)
       b. Check conservation
       c. If t ∈ output_times: save full fields
       d. Update interface tracking
       e. Check termination criteria

3. FINALIZE
   a. Compute final diagnostics
   b. Package solution
   c. Return PyrolysisSolution
```

## 7.2 Physics Step Detail

```
ALGORITHM: Physics RHS Evaluation

INPUT: Current state u = (T, ξ₁, ..., ξ_NC), mesh M, time t
OUTPUT: du/dt

1. PROPERTY EVALUATION (for each cell i)
   a. ρᵢ = mixture_density(ξᵢ, T_i)
   b. cₚ,ᵢ = effective_heat_capacity(ξᵢ, Tᵢ)
   c. kᵢ = effective_conductivity(ξᵢ, Tᵢ, β)
   d. λᵢ = effective_gas_transfer(ξᵢ, Tᵢ)
   e. εᵢ = effective_emissivity(ξᵢ)
   f. αᵢ = effective_absorption(ξᵢ)

2. REACTION RATES (for each cell i, reaction j)
   a. rⱼ,ᵢ = Aⱼ·exp(-Eⱼ/RTᵢ)·ξᵣ₁^n₁·ξᵣ₂^n₂
   b. Apply temperature limits
   c. Compute heat generation: Qᵢ = Σⱼ hⱼ(Tᵢ)·rⱼ,ᵢ

3. FACE FLUXES
   For each interior face f between cells L and R:
   
   3a. CONDUCTIVE FLUX
       k_f = harmonic_mean(k_L, k_R)
       q_cond,f = -k_f · (T_R - T_L) / δ_LR
   
   3b. RADIATIVE FLUX (Beer-Lambert)
       Compute I(z) by integrating from surface
       q_rad,f = I(z_f) - I(z_f+δ)
   
   3c. GAS CONVECTIVE FLUX
       For each gas component g:
       J_g,f = -ρ_g,f · λ_f · ∂(ξ_g/ρ_g)/∂z |_f
       Use upwind for stability
       q_conv,f = Σ_g c_p,g · J_g,f · (T_R - T_L)   (advection form)
       Note: Boundary faces have q_conv = 0 (advection only redistributes heat internally)

4. BOUNDARY FLUXES
   For each boundary face f:
   
   4a. THERMAL BC
       If adiabatic: q_f = 0
       If heat flux: q_f = q"(t)
       If convective: q_f = h(T_∞ - T_s)
       If radiative: q_f = εσ(T_∞⁴ - T_s⁴)
       If combined: sum all contributions
       If flame: q_f = q_base + q_flame·H(CI - 1)
   
   4b. MASS BC
       If impermeable: J_g,f = 0
       If convective: J_g,f = h_m(ξ_∞ - ξ_s)

5. ALE CONVECTIVE FLUXES
   For each face f with mesh velocity w_f:
   
   5a. ENERGY
       q_ale,f = (v_g - w) · ρcₚT |_f
       Use upwind based on sign(v_g - w)
   
   5b. SPECIES
       J_ale,j,f = (vⱼ - w) · ξⱼ |_f
       For solids: vⱼ = w → J_ale = 0
       For gases: use upwind

6. ASSEMBLE CELL RESIDUALS
   For each cell i:
   
   6a. ENERGY
       dTᵢ/dt = 1/(ρcₚ)ᵢ · [
           (1/Vᵢ)·Σ_faces(q_cond - q_ale) · A_f
           + q_rad,i"' 
           + Qᵢ
       ]
   
   6b. SPECIES (component j)
       dξⱼ,ᵢ/dt = [
           Σ_reactions νⱼ,ᵣ · rᵣ,ᵢ
           - (1/Vᵢ)·Σ_faces(J_j + J_ale,j) · A_f
       ]

7. RETURN du/dt vector
```

## 7.3 Mesh Motion Algorithm

```
ALGORITHM: Compute Mesh Velocity and Move Mesh

INPUT: Current mesh M, state u, dt
OUTPUT: Updated mesh M'

1. COMPUTE MATERIAL VELOCITY
   For each cell i:
   a. Volume change rate: dV/dt|_i = Vᵢ · Σⱼ (1/ρⱼ) · dξⱼ/dt
   b. Local expansion: θᵢ = (1/Vᵢ) · dV/dt
   
2. COMPUTE SURFACE VELOCITY
   At top boundary:
   a. Mass flux out: ṁ = Σ_g J_g · A_surface
   b. Surface regression: dz_s/dt = -ṁ / ρ_surface
   c. Or from volume change integral

3. COMPUTE INTERIOR MESH VELOCITY
   
   Method A: Laplacian Smoothing
   Solve: ∇²w = 0 in Ω
          w = w_surface on Γ_top
          w = 0 on Γ_bottom
   
   Method B: Spring Analogy
   For each interior node n:
   w_n = Σ_neighbors k_{n,m}(w_m - w_n) / Σ k_{n,m}
   Iterate to convergence
   
   Method C: Linear Elasticity
   Solve: ∇·σ = 0
   With σ = f(∇w), stiffer near boundaries
   
   Method D: Weighted ALE
   w = α_L · w_material + α_S · w_smooth + α_I · w_interface
   Where w_material tracks local volume change
         w_smooth from Laplacian
         w_interface tracks reaction front

4. VERIFY GEOMETRIC CONSERVATION LAW
   For each cell i:
   Check: dVᵢ/dt ≈ Σ_faces w·n·A |_f
   If violation > tolerance:
      a. Adjust mesh velocity
      b. Or use GCL-compliant discretization

5. MOVE NODES
   For each node n:
   x_n^(new) = x_n^(old) + dt · w_n
   
6. UPDATE GEOMETRY
   a. Recompute cell volumes
   b. Recompute face areas and normals
   c. Recompute cell centers
   d. Update geometric cache

7. CHECK MESH QUALITY
   For each cell i:
   a. Check aspect ratio < AR_max
   b. Check Jacobian > 0
   c. Check smoothness
   If failed:
      a. Reduce dt and retry
      b. Or trigger adaptivity

8. REMAP SOLUTION (if significant motion)
   If max(|x^new - x^old|) / h_min > remap_threshold:
   Use conservative remapping algorithm
```

## 7.4 h-Refinement Algorithm

```
ALGORITHM: Adaptive h-Refinement (1D)

INPUT: Mesh M, error indicators η, AMR parameters
OUTPUT: Refined mesh M'

1. COMPUTE ERROR INDICATORS
   For each active cell i:
   η_T,i = h_i · |∇T|_i
   η_ξ,i = h_i · max_j(|∇ξ_j|_i)
   η_r,i = h_i · max_j(|r_j|_i / max(ξ_j, ε))
   η_i = combine(η_T, η_ξ, η_r)

2. MARK CELLS
   η_max = max(η)
   For each cell i:
   If η_i > θ_refine · η_max AND h_i > h_min AND level_i < max_level:
      Mark i for REFINEMENT
   ElseIf η_i < θ_coarsen · η_max AND h_i < h_max AND level_i > 0:
      Mark i for COARSENING

3. ENSURE 2:1 BALANCE
   Repeat until no changes:
   For each cell i marked for refinement:
      For each neighbor n of i:
         If level_n < level_i:
            Mark n for refinement

4. REFINE MARKED CELLS
   For each cell i marked for REFINEMENT:
   a. Create new node at cell center
   b. Create two child cells
   c. Initialize child states by interpolation:
      - T_child = interpolate(T_parent, T_neighbors)
      - ξ_child = interpolate(ξ_parent, ξ_neighbors)
      - Ensure conservation: V_p·φ_p = V_a·φ_a + V_b·φ_b
   d. Update face connectivity
   e. Mark parent as inactive

5. COARSEN MARKED CELLS
   For each pair of siblings both marked for COARSENING:
   a. Verify both are leaves (no children)
   b. Verify parent exists
   c. Merge states conservatively:
      φ_parent = (V_a·φ_a + V_b·φ_b) / V_parent
   d. Delete child cells and intermediate node
   e. Mark parent as active

6. UPDATE MESH DATA STRUCTURES
   a. Renumber cells and nodes
   b. Update connectivity arrays
   c. Update boundary face lists
   d. Recompute geometric quantities

7. RETURN refined mesh M'
```

## 7.5 r-Refinement Algorithm

```
ALGORITHM: r-Refinement (Equidistribution)

INPUT: Mesh M, state u, monitor function parameters
OUTPUT: Mesh M' with relocated nodes

1. COMPUTE MONITOR FUNCTION
   For each cell i:
   ω_i = sqrt(1 + α_T·|∂T/∂z|² + α_ξ·Σ_j|∂ξ_j/∂z|²)
   
   Or interface-focused:
   ω_i = ω_base + ω_interface · exp(-d_i² / σ²)
   Where d_i = distance to nearest interface

2. COMPUTE TARGET POSITIONS
   Total integral: Ω = ∫₀^L ω(z) dz ≈ Σᵢ ω_i · Δz_i
   
   For each node n (interior only):
   Find z_n^target such that:
   ∫₀^{z_n^target} ω dz = (n/N) · Ω

3. ITERATIVE RELAXATION
   For iter = 1 to max_iter:
      For each interior node n:
         z_n^new = z_n^old + β · (z_n^target - z_n^old)
      
      Check convergence:
      If max|z^new - z^old| < tol: break
      
      Recompute monitor function on new mesh
      Recompute target positions

4. CHECK QUALITY CONSTRAINTS
   For each cell i:
   If Δz_i < Δz_min OR Δz_i / Δz_{i±1} > ratio_max:
      Apply limiting to prevent bad cells

5. CONSERVATIVE REMAP
   Use remapping algorithm to transfer solution

6. RETURN mesh M' with new node positions
```

## 7.6 Conservative Remapping Algorithm

```
ALGORITHM: Conservative Solution Remapping

INPUT: Old mesh M_old, new mesh M_new, old solution u_old
OUTPUT: New solution u_new

1. FOR EACH NEW CELL i_new:
   
   a. Find overlapping old cells:
      overlaps = {(i_old, V_overlap) : cell i_old overlaps i_new}
   
   b. For temperature (intensive):
      T_new = Σ (V_overlap / V_new) · T_old
      
      Or with gradient reconstruction (2nd order):
      T_new = Σ (V_overlap / V_new) · 
              [T_old + ∇T_old · (x_new - x_old)]
   
   c. For concentrations (extensive per volume):
      M_j,new = Σ V_overlap · ξ_j,old  (total mass)
      ξ_j,new = M_j,new / V_new
   
   d. Apply limiters to prevent new extrema:
      T_new = limit(T_new, min(T_neighbors), max(T_neighbors))
      ξ_new = max(ξ_new, 0)  # non-negative concentrations

2. VERIFY CONSERVATION
   Total mass before: M_old = Σ_i V_i · ρ_i |_old
   Total mass after:  M_new = Σ_i V_i · ρ_i |_new
   
   If |M_new - M_old| / M_old > tol:
      Apply mass correction:
      ξ_corrected = ξ_new · (M_old / M_new)

3. RETURN u_new
```

---

# 8. Numerical Methods

## 8.1 Time Integration

### Recommended: IMEX (Implicit-Explicit)

Treat stiff terms (reactions, diffusion) implicitly, non-stiff (convection, sources) explicitly:

```
u^(n+1) = u^n + dt·[f_exp(u^n) + f_imp(u^(n+1))]
```

**Implementation**: Use DifferentialEquations.jl with:
- `KenCarp4()` or `ARKODE()` for IMEX
- `Rodas5P()` or `QNDF()` for fully implicit
- `CVODE_BDF()` for high-accuracy stiff problems

### Splitting Approach

Alternative: operator splitting between physics and mesh motion

**Strang Splitting** (2nd order):
```
1. Mesh half-step: M^(n+1/2) = move_mesh(M^n, dt/2)
2. Physics full-step: u^(n+1) = physics_step(u^n, M^(n+1/2), dt)
3. Mesh half-step: M^(n+1) = move_mesh(M^(n+1/2), dt/2)
```

## 8.2 Spatial Discretization

### Finite Volume Basics

Cell-averaged values:
```
φ̄ᵢ = (1/Vᵢ) ∫_Vᵢ φ dV
```

Update equation:
```
dφ̄ᵢ/dt = -(1/Vᵢ) Σ_faces F_f · A_f + S̄ᵢ
```

### Flux Computation

**Diffusive fluxes** (e.g., heat conduction):
```
F_diff,f = -k_f · (φ_R - φ_L) / d_LR
```

Where k_f = harmonic mean for heterogeneous media:
```
k_f = 2 k_L k_R / (k_L + k_R)
```

**Convective fluxes** (ALE terms):
```
F_conv,f = (v - w)_f · φ_f
```

With upwinding:
```
φ_f = φ_L  if (v-w)·n > 0
    = φ_R  otherwise
```

### Gradient Reconstruction (2nd order)

For cell i with neighbors, solve least-squares:
```
∇φ|_i = argmin Σ_neighbors w_n · |φ_n - φ_i - ∇φ·(x_n - x_i)|²
```

### Flux Limiters (TVD)

Apply limiter to prevent oscillations:
```
φ_f = φ_L + ψ(r) · (φ_R - φ_L) / 2
```

Where r = (φ_L - φ_{LL}) / (φ_R - φ_L) and:

- Minmod: ψ(r) = max(0, min(1, r))
- Van Leer: ψ(r) = (r + |r|) / (1 + |r|)
- Superbee: ψ(r) = max(0, min(2r, 1), min(r, 2))

## 8.3 Stability Constraints

### CFL Condition (convection)
```
dt < CFL · h_min / |v - w|_max
```
Typical CFL = 0.5

### Diffusion Stability (if explicit)
```
dt < h_min² / (2 · α_max)
```
Where α = k / (ρ cₚ)

### Reaction Timescale
```
dt < 1 / (∂r/∂ξ)_max
```
Or use implicit integration

### Mesh Motion Stability
```
dt < h_min / |w|_max
```

### Combined Limit
```
dt = min(dt_cfl, dt_diff, dt_react, dt_mesh) · safety_factor
```
Typical safety = 0.9

## 8.4 Linear Algebra

### Sparse Matrix Storage

For implicit methods, Jacobian J = ∂f/∂u is sparse:
- Use CSR format (SparseMatrixCSC in Julia)
- Bandwidth ~ 3·N_components for 1D

### Linear Solvers

For Ax = b:
- Direct: LU factorization (small problems)
- Iterative: GMRES with ILU preconditioner (large problems)

Use LinearSolve.jl for automatic selection.

### Jacobian Computation

Options:
1. **Automatic differentiation**: ForwardDiff.jl (exact, type-stable)
2. **Finite differences**: Sparse coloring to minimize function evaluations
3. **Analytical**: Hand-derived, most efficient but error-prone

Recommend: Start with ForwardDiff, optimize later if needed.

---

# 9. Implementation Guidelines

## 9.1 Type Stability

**Critical for performance**. All functions in hot paths must be type-stable.

```julia
# BAD: type unstable
function compute_flux(cell)
    if cell isa Cell1D
        return flux_1d(cell)
    else
        return flux_2d(cell)  # Different return type
    end
end

# GOOD: dispatch-based
compute_flux(cell::Cell1D) = flux_1d(cell)
compute_flux(cell::Cell2D) = flux_2d(cell)
```

Use `@code_warntype` to verify type stability.

## 9.2 Allocation Avoidance

Pre-allocate all arrays used in loops:

```julia
# BAD: allocates in loop
function evaluate_properties!(mesh)
    for cell in mesh.cells
        props = compute_props(cell)  # Allocates new array
        cell.properties = props
    end
end

# GOOD: pre-allocated cache
struct PropertyCache
    k::Vector{Float64}
    cp::Vector{Float64}
    ρ::Vector{Float64}
end

function evaluate_properties!(cache::PropertyCache, mesh)
    @inbounds for i in eachindex(mesh.cells)
        cache.k[i] = compute_k(mesh.cells[i])
        cache.cp[i] = compute_cp(mesh.cells[i])
        cache.ρ[i] = compute_rho(mesh.cells[i])
    end
end
```

## 9.3 SIMD and Threading

Use `@simd` for simple loops:
```julia
@inbounds @simd for i in 1:n
    result[i] = a[i] * b[i]
end
```

Use `Threads.@threads` for embarrassingly parallel:
```julia
Threads.@threads for i in 1:n_cells
    compute_cell_residual!(residuals, i, mesh, state)
end
```

## 9.4 StaticArrays for Small Fixed-Size Data

```julia
using StaticArrays

# For up to ~10 components, use static arrays
struct CellState{NC}
    T::Float64
    ξ::SVector{NC,Float64}  # Stack-allocated, no heap
end

# Operations are fast and allocation-free
@inline function mix_density(state::CellState{NC}, material) where NC
    ρ_inv = zero(Float64)
    @inbounds for j in 1:NC
        ρ_inv += state.ξ[j] / material.components[j].density(state.T)
    end
    return 1 / ρ_inv
end
```

## 9.5 Function Barriers

Isolate type-unstable code with function barriers:

```julia
# Main function is type-stable
function solve!(problem)
    # Type-unstable setup
    mesh = setup_mesh(problem.mesh_spec)  # Could return Mesh1D or Mesh2D
    
    # Function barrier
    _solve_impl!(problem, mesh)
end

# Inner function is specialized for concrete mesh type
function _solve_impl!(problem, mesh::Mesh1D)
    # All code here is type-stable
end
```

## 9.6 Testing Strategy

### Unit Tests
- Each module should have corresponding test file
- Test property functions: verify values at known temperatures
- Test reactions: verify mass conservation in stoichiometry
- Test mesh: verify geometric quantities

### Integration Tests
- Full problem setup and single timestep
- Verify conservation properties
- Test boundary conditions

### Verification Tests (MMS)
- Method of Manufactured Solutions
- Verify expected convergence rates

### Validation Tests
- Compare against ThermaKin results
- Compare against analytical solutions

## 9.7 Documentation

Use docstrings for all public functions:

```julia
"""
    compute_mesh_velocity!(w, mesh, state, method)

Compute the mesh velocity field using the specified method.

# Arguments
- `w::Vector{Float64}`: Output mesh velocity at each node
- `mesh::AbstractMesh`: Current mesh configuration
- `state::StateVector`: Current solution state
- `method::AbstractMeshVelocityMethod`: Velocity computation method

# Returns
Nothing. Modifies `w` in-place.

# Notes
The mesh velocity satisfies the geometric conservation law:
``dV_i/dt = \\sum_f w \\cdot n \\cdot A_f``
"""
function compute_mesh_velocity!(w, mesh, state, method)
    # Implementation
end
```

---

# 10. Verification & Validation

## 10.1 Analytical Solutions for Verification

### Heat Conduction in Semi-Infinite Solid

**Setup**: Semi-infinite solid, initially at T₀, surface (at z = L) suddenly raised to T_s

**Solution** (with z measured from bottom, surface at z = L):
```
T(z,t) = T_s + (T₀ - T_s) · erf((L - z) / (2√(αt)))
```

Note: The argument `(L - z)` represents distance from the heated surface.

**Test**: 
- Compare numerical vs analytical at various times
- Verify 2nd-order convergence in space and time

### Beer-Lambert Radiation Absorption

**Setup**: Semi-transparent medium with absorption α, incident flux I₀ entering at surface (z = L)

**Solution** (radiation propagates from z = L toward z = 0):
```
I(z) = I₀ · exp(-α(L - z))
q'''(z) = α · I₀ · exp(-α(L - z))
```

Note: The argument `(L - z)` represents depth from the irradiated surface.

**Test**:
- Verify exponential decay of intensity from surface into depth
- Verify energy conservation: ∫q'''dz = I₀(1 - exp(-αL))

### Single First-Order Reaction

**Setup**: Isothermal decomposition, A → B

**Solution**:
```
ξ_A(t) = ξ_A,0 · exp(-k·t)
ξ_B(t) = ξ_A,0 · (1 - exp(-k·t))
```

**Test**:
- Verify mass conservation: ξ_A + ξ_B = ξ_A,0
- Verify exponential decay rate

### Diffusion from Thin Layer

**Setup**: Thin gas layer at z = L (surface), diffuses into domain

**Solution** (Gaussian, with distance from surface = L - z):
```
ξ(z,t) = (M₀ / √(4πDt)) · exp(-(L-z)² / (4Dt))
```

**Test**:
- Verify mass conservation
- Verify spreading rate from surface into material

## 10.2 Method of Manufactured Solutions (MMS)

Construct analytical solution and derive source terms:

1. Choose smooth manufactured solution: T̂(z,t), ξ̂(z,t)
2. Substitute into governing equations
3. Compute required source term S(z,t)
4. Solve numerical problem with source term
5. Compare numerical to manufactured solution

**Example manufactured solution**:
```
T̂(z,t) = T₀ + ΔT·sin(πz/L)·(1 - exp(-t/τ))
```

**Expected results**:
- 2nd-order schemes: error ∝ h²
- 1st-order schemes: error ∝ h

## 10.3 Conservation Checks

### Mass Conservation
```
dM_total/dt = -ṁ_surface + ṁ_bottom
```
Where ṁ = ∫ J·n dA

For impermeable bottom:
```
M(t) + ∫₀ᵗ ṁ_surface dt' = M(0)
```

### Energy Conservation
```
dE_total/dt = Q̇_in - Q̇_out + Q̇_reaction
```

Verify:
```
∫ρcₚT dV + ∫₀ᵗ q_surface dt' = E(0) + ∫₀ᵗ∫ Q_rxn dV dt'
```

### GCL Verification
```
V_i^(n+1) - V_i^n = dt · Σ_faces w·n·A
```
Should hold to machine precision.

## 10.4 Validation Against ThermaKin

### Test Case 1: PMMA Gasification
- Material: PMMA (well-characterized)
- Conditions: 50 kW/m² external flux
- Outputs: MLR, surface temperature
- Reference: ThermaKin publications

### Test Case 2: Charring Material
- Material: Wood analog
- Multiple reactions, char formation
- Track char layer thickness

### Test Case 3: Intumescent Coating
- Swelling material
- Test ALE mesh motion
- Compare thickness evolution

## 10.5 Convergence Studies

### Spatial Convergence
Run with h, h/2, h/4, h/8
Plot error vs h on log-log
Verify slope matches expected order

### Temporal Convergence
Run with dt, dt/2, dt/4, dt/8
Plot error vs dt on log-log

### Combined Space-Time
Refine both simultaneously
Verify combined error reduction

---

# 11. Performance Considerations

## 11.1 Expected Computational Cost

### Per-Timestep Operations
- Property evaluation: O(N · N_c) where N = cells, N_c = components
- Reaction rates: O(N · N_r) where N_r = reactions
- Flux computation: O(N)
- Matrix assembly: O(N · bw) where bw = bandwidth
- Linear solve: O(N · bw²) direct, O(N · k · bw) iterative (k = iterations)

### AMR Overhead
- Error indicator: O(N)
- Refinement/coarsening: O(N_modified)
- Remapping: O(N)

### ALE Overhead
- Mesh velocity: O(N) for Laplacian, O(N · k) for iterative methods
- Geometry update: O(N)
- GCL check: O(N)

## 11.2 Memory Requirements

### Per Cell
- State: (1 + N_c) × 8 bytes (T + concentrations)
- Properties: ~10 × 8 bytes (cached evaluations)
- Geometry: ~5 × 8 bytes (volume, centroid, etc.)
- Connectivity: ~8 × 4 bytes (indices)

### Total Estimate
```
Memory ≈ N × (N_c + 30) × 8 bytes + sparse matrix storage
```

For N = 1000, N_c = 10:
```
Memory ≈ 1000 × 40 × 8 = 320 KB + O(100 KB) for matrices
```

## 11.3 Parallelization Opportunities

### Thread-Level
- Property evaluation: embarrassingly parallel over cells
- Flux computation: parallel over faces (with care for race conditions)
- Reaction rates: parallel over cells

### Vectorization
- Property function evaluation (same function, different T)
- Flux assembly (regular memory access patterns)

### GPU Potential (future)
- Large problems (N > 10⁴) could benefit from GPU
- Would require rewriting core loops with CUDA.jl or KernelAbstractions.jl

## 11.4 Optimization Priorities

1. **Type stability first**: Biggest impact, check with @code_warntype
2. **Eliminate allocations**: Second biggest, check with @allocated
3. **SIMD**: Use @simd where applicable
4. **Cache efficiency**: Access arrays in column-major order
5. **Threading**: Only for large problems, overhead can dominate

## 11.5 Profiling Workflow

```julia
using Profile, ProfileView

# Warm up
solve(problem_small)

# Profile
@profile solve(problem)
ProfileView.view()

# Check allocations
@time solve(problem)
@allocated solve(problem)

# Benchmark critical functions
using BenchmarkTools
@btime compute_fluxes!($flux_cache, $mesh, $state)
```

---

# 12. API Design

## 12.1 High-Level API

```julia
using Pyrolysis

# Define material
# Note: κ=0.0 (default) means diffusion-only gas transport (no Darcy flow)
pmma = Material(
    name = :PMMA,
    components = (
        SolidComponent(:virgin, ρ=1190, c=1420, k=0.21, ε=0.86),  # κ=0 default
        GaseousComponent(:MMA, M=0.100, c=1500, k=0.02, λ=1e-5)
    ),
    reactions = (
        Reaction(:decomposition,
            :virgin => :MMA,
            A = 8.5e12, E = 188e3, h = -870e3),
    )
)

# Define geometry
mesh = create_mesh(
    Slab(thickness=0.01),  # 1 cm thick
    n_cells = 100,
    stretching = :geometric,
    stretch_factor = 1.05
)

# Define boundary conditions
top_bc = HeatFluxBC(t -> 50e3) + RadiativeBC(ε=0.9, T_inf=300)
bottom_bc = AdiabaticBC()

# Create problem
problem = PyrolysisProblem(
    mesh = mesh,
    material = pmma,
    thermal_bcs = Dict(1 => top_bc, 2 => bottom_bc),
    mass_bcs = Dict(1 => ConvectiveMassBC(h_m=0.01), 2 => ImpermeableBC()),
    T_initial = z -> 300.0,
    tspan = (0.0, 600.0)
)

# Configure ALE
set_ale!(problem,
    method = WeightedALE(α_lagrangian=0.3, α_smoothing=0.7),
    surface_tracking = true
)

# Configure AMR
set_amr!(problem,
    indicators = [GradientIndicator(:T), ReactionRateIndicator()],
    h_min = 1e-5,
    h_max = 1e-3,
    max_levels = 4
)

# Solve
solution = solve(problem, 
    solver = KenCarp4(),
    dt_initial = 0.01,
    callbacks = [MassLossCallback(), SaveEvery(10.0)]
)

# Post-process
plot_mlr(solution)
plot_temperature_profile(solution, times=[60, 120, 300, 600])
export_vtk(solution, "pmma_pyrolysis")
```

## 12.2 Darcy-Fick API (Pressure-Driven Gas Transport)

For simulations with pressure-driven gas transport, specify permeability `κ` on each
solid component. Darcy flow is automatically enabled when any component has κ > 0.

```julia
using Pyrolysis

# Define charring material with per-component permeability
# κ controls resistance to pressure-driven gas flow through the matrix
material = Material(
    name = :CharringMaterial,
    components = (
        # Virgin: tight cellular structure, low permeability
        SolidComponent(:virgin, ρ=500, c=2000, k=0.15, ε=0.9, κ=1e-15),
        # Char: porous with open channels, high permeability  
        SolidComponent(:char, ρ=150, c=1500, k=0.10, ε=0.95, κ=1e-12),
        # Gas component (no κ - gases don't have permeability)
        GaseousComponent(:volatiles, M=0.044, c=1200, k=0.03, λ=2e-5)
    ),
    reactions = (
        Reaction(:pyrolysis, 1 => (2, 3),
            stoich = (1.0, 0.0, 0.20, 0.80),  # 80% mass to gas
            A = 1e13, E = 150e3, h = 300e3),
    ),
    # Permeability mixing: how to combine component κ values
    permeability_mixing = WEIGHTED,
    β_κ = 0.5,
)

# Create mesh and problem
mesh = create_uniform_mesh(0.01, 50, Val(3))
problem = PyrolysisProblem(
    mesh = mesh,
    material = material,
    thermal_bcs = Dict(1 => HeatFluxBC(75e3) + RadiativeBC(ε=0.9, T_inf=300),
                       2 => AdiabaticBC()),
    mass_bcs = Dict(1 => DiffusiveMassBC(0.01), 2 => ImpermeableBC()),
    T_initial = z -> 300.0,
    tspan = (0.0, 120.0)
)

# Solve - Darcy flow automatically enabled since components have κ > 0
solution = solve(problem,
    use_ale = true,
    min_thickness = 0.5e-3
)

# Effective permeability κ_eff is computed from component values via mixing rules
# Gas velocity: v_g = -(κ_eff/μ)·∇P where P is from ideal gas law
```

## 12.3 Lower-Level API

For users needing more control:

```julia
# Manual mesh creation
mesh = Mesh1D()
add_node!(mesh, 0.0, boundary=:top)
for i in 1:99
    add_node!(mesh, i * 0.01 / 100)
end
add_node!(mesh, 0.01, boundary=:bottom)
generate_cells!(mesh)
generate_faces!(mesh)

# Manual property evaluation
cache = PropertyCache(mesh)
evaluate_properties!(cache, mesh, material)

# Manual timestep
dt = compute_timestep(mesh, cache)
rhs = compute_rhs(mesh, cache, material, bcs)
state_new = time_step(state, rhs, dt, solver)

# Manual mesh motion
w = compute_mesh_velocity(mesh, state, method)
move_mesh!(mesh, w, dt)
state_new = remap(state, mesh_old, mesh_new)
```

## 12.4 Extension Points

### Custom Property Functions
```julia
struct MyCustomProperty <: AbstractPropertyFunction
    params::Vector{Float64}
end

(p::MyCustomProperty)(T) = custom_formula(T, p.params)
```

### Custom Mesh Velocity Methods
```julia
struct MyMeshMethod <: AbstractMeshVelocityMethod
    params::Float64
end

function compute_mesh_velocity!(w, mesh, state, method::MyMeshMethod)
    # Custom implementation
end
```

### Custom Error Indicators
```julia
struct MyIndicator <: AbstractErrorIndicator end

function compute_indicator(ind::MyIndicator, mesh, state)
    # Return vector of error indicators per cell
end
```

## 12.5 Callbacks and Events

```julia
# Mass loss rate callback (called every timestep)
struct MassLossCallback <: AbstractCallback end

function (cb::MassLossCallback)(integrator)
    mlr = compute_mass_loss_rate(integrator)
    push!(integrator.sol.mlr, mlr)
end

# Event: surface temperature threshold
struct IgnitionEvent <: AbstractEvent
    T_ignition::Float64
end

function condition(event::IgnitionEvent, u, t, integrator)
    T_surface = u[1]  # First DOF is surface temperature
    return T_surface - event.T_ignition
end

function affect!(event::IgnitionEvent, integrator)
    println("Ignition at t = $(integrator.t)")
    enable_flame_bc!(integrator)
end
```

---

# 13. Dependencies & Compatibility

## 13.1 Required Dependencies

```toml
[deps]
# Core
LinearAlgebra = "37e2e46d-..."
SparseArrays = "2f01184e-..."

# Performance
StaticArrays = "90137ffa-..."  # Fixed-size arrays
Parameters = "d96e819e-..."    # @with_kw macro

# Time integration
DifferentialEquations = "0c46a032-..."  # ODE solvers
# OR individual packages:
# OrdinaryDiffEq = "1dea7af3-..."
# Sundials = "c3572dad-..."

# Automatic differentiation (for Jacobians)
ForwardDiff = "f6369f11-..."

# Linear solvers
LinearSolve = "7ed4a6bd-..."

# Optional: Visualization
RecipesBase = "3cdcf5f2-..."  # Plotting recipes

# Optional: Parallel
# Polyester = "f517fe37-..."  # Fast threading
```

## 13.2 Compatibility

- Julia ≥ 1.9 (for improved compilation and threading)
- 64-bit systems (Float64 throughout)
- Cross-platform (Linux, macOS, Windows)

## 13.3 Optional Dependencies

```toml
[extras]
# Testing
Test = "8dfed614-..."
SafeTestsets = "1bc83da4-..."

# Visualization (optional)
Plots = "91a5bcdd-..."
Makie = "ee78f7c6-..."

# I/O
WriteVTK = "64499a7a-..."  # VTK export
JLD2 = "033835bb-..."      # Save/load solutions

# Profiling
BenchmarkTools = "6e4b80f9-..."
ProfileView = "c46f51b8-..."
```

## 13.4 Build & Test

```bash
# Clone and setup
git clone https://github.com/user/Pyrolysis.jl
cd Pyrolysis.jl
julia --project -e 'using Pkg; Pkg.instantiate()'

# Run tests
julia --project -e 'using Pkg; Pkg.test()'

# Build documentation
julia --project=docs docs/make.jl
```

---

# Appendix A: Physical Constants

```julia
module Constants

# Fundamental constants
const R_GAS = 8.31446261815324           # J/(mol·K) - gas constant
const STEFAN_BOLTZMANN = 5.670374419e-8  # W/(m²·K⁴) - Stefan-Boltzmann
const AVOGADRO = 6.02214076e23           # mol⁻¹
const BOLTZMANN = 1.380649e-23           # J/K

# Reference conditions
const T_REF = 298.15                     # K - standard temperature
const P_REF = 101325.0                   # Pa - standard pressure

# Typical material property ranges (for validation)
const RHO_POLYMER_RANGE = (800.0, 2000.0)      # kg/m³
const CP_POLYMER_RANGE = (1000.0, 3000.0)      # J/(kg·K)
const K_POLYMER_RANGE = (0.1, 0.5)             # W/(m·K)
const EMISSIVITY_RANGE = (0.7, 1.0)            # -

# Typical kinetic parameter ranges
const A_RANGE = (1e8, 1e20)                    # s⁻¹
const E_RANGE = (1e5, 3e5)                     # J/mol

end
```

---

# Appendix B: Typical Material Properties

## B.1 Averaged Polymer Properties (from ThermaKin database)

Based on 97 polymers:

| Property | Mean | Std Dev | Units |
|----------|------|---------|-------|
| Density | 1300 | 300 | kg/m³ |
| Heat capacity | 2300 | 800 | J/(kg·K) |
| Thermal conductivity | 0.24 | 0.10 | W/(m·K) |
| Absorption coefficient | 2.9 | 2.0 | m²/kg |
| Emissivity | 0.88 | 0.05 | - |

## B.2 Averaged Decomposition Kinetics (35 polymers)

| Parameter | Mean | Range | Units |
|-----------|------|-------|-------|
| Pre-exponential A | 1.6×10¹⁶ | 10⁸-10²⁰ | s⁻¹ |
| Activation energy E | 218 | 100-300 | kJ/mol |
| Heat of reaction | -1300 | -3000 to +500 | kJ/kg |
| Char yield | 0.1 | 0-0.5 | - |

## B.3 Example Materials

### PMMA (Poly(methyl methacrylate))
```julia
pmma = Material(
    components = (
        SolidComponent(:pmma, ρ=1190, c=LinearProperty(1420, 2.5), k=0.21, ε=0.86, α=1.5),
        GaseousComponent(:MMA, M=0.100, c=1500, k=0.02, λ=1e-5)
    ),
    reactions = (
        Reaction(:depolymerization, :pmma => :MMA,
            A=8.5e12, E=188e3, h=ConstantProperty(870e3)),  # positive = endothermic
    ),
    β = 0.5
)
```

### Wood (Simplified 3-Component)
```julia
wood = Material(
    components = (
        SolidComponent(:cellulose, ρ=1500, c=1500, k=0.25),
        SolidComponent(:char, ρ=200, c=1000, k=0.1),
        GaseousComponent(:volatiles, M=0.044, c=1000, k=0.025, λ=1e-5),
        GaseousComponent(:water, M=0.018, c=2000, k=0.025, λ=2e-5)
    ),
    reactions = (
        Reaction(:drying, :cellulose => :water, A=1e13, E=100e3, h=2.4e6),  # positive = endothermic
        Reaction(:pyrolysis, :cellulose => (:char, :volatiles), 
            stoich=(1.0, 0.0, 0.2, 0.8), A=1e15, E=200e3, h=500e3),  # positive = endothermic
    )
)
```

---

# Appendix C: Error Messages and Troubleshooting

## C.1 Common Errors

### "Negative concentration detected"
- **Cause**: Numerical instability, often from too large timestep
- **Fix**: Reduce dt, use implicit solver, or add concentration limiters

### "Mesh quality degraded below threshold"
- **Cause**: Large deformations causing cell distortion
- **Fix**: Increase mesh smoothing weight, reduce Lagrangian weight, or trigger remeshing

### "GCL violation exceeds tolerance"
- **Cause**: Inconsistent mesh motion and flux computation
- **Fix**: Use GCL-compliant discretization, reduce dt

### "Linear solver failed to converge"
- **Cause**: Ill-conditioned Jacobian, often from very stiff reactions
- **Fix**: Use better preconditioner, tighter tolerances, or analytical Jacobian

### "Timestep became too small"
- **Cause**: Stiffness from fast reactions or sharp gradients
- **Fix**: Use stiffer solver (BDF), local timestepping, or adaptive mesh to resolve gradients

## C.2 Performance Issues

### "Simulation running slowly"
1. Check type stability with `@code_warntype`
2. Check allocations with `@allocated`
3. Profile with `@profile` and `ProfileView`
4. Consider coarser mesh or larger dt

### "Memory usage growing"
1. Clear solution history if not needed
2. Use callbacks instead of storing all timesteps
3. Check for array growth in custom code

---

# Appendix D: References

## D.1 ThermaKin Papers

1. Stoliarov, S.I., Lyon, R.E. (2008). "Thermo-kinetic model of burning for pyrolyzing materials." Fire Safety Science, 9, 1141-1152.

2. Stoliarov, S.I., Crowley, S., Lyon, R.E., Linteris, G.T. (2009). "Prediction of the burning rates of non-charring polymers." Combustion and Flame, 156(5), 1068-1083.

3. Stoliarov, S.I., Crowley, S., Walters, R.N., Lyon, R.E. (2010). "Prediction of the burning rates of charring polymers." Combustion and Flame, 157(11), 2024-2034.

4. Li, J., Stoliarov, S.I. (2013). "Measurement of kinetics and thermodynamics of the thermal degradation for charring polymers." Polymer Degradation and Stability, 98(12), 2484-2494.

5. Li, J., Gong, J., Stoliarov, S.I. (2014). "Gasification experiments for pyrolysis model parameterization and validation." International Journal of Heat and Mass Transfer, 77, 738-744.

## D.2 ALE Methods

6. Donea, J., Huerta, A., Ponthot, J.P., Rodríguez-Ferran, A. (2004). "Arbitrary Lagrangian-Eulerian methods." Encyclopedia of Computational Mechanics.

7. Farhat, C., Geuzaine, P., Grandmont, C. (2001). "The discrete geometric conservation law and the nonlinear stability of ALE schemes for the solution of flow problems on moving grids." Journal of Computational Physics, 174(2), 669-694.

## D.3 Adaptive Mesh Refinement

8. Berger, M.J., Colella, P. (1989). "Local adaptive mesh refinement for shock hydrodynamics." Journal of Computational Physics, 82(1), 64-84.

9. Huang, W., Russell, R.D. (2011). "Adaptive moving mesh methods." Springer.

## D.4 Finite Volume Methods

10. LeVeque, R.J. (2002). "Finite volume methods for hyperbolic problems." Cambridge University Press.

11. Versteeg, H.K., Malalasekera, W. (2007). "An introduction to computational fluid dynamics: The finite volume method." Pearson.

---

# Appendix E: Glossary

| Term | Definition |
|------|------------|
| ALE | Arbitrary Lagrangian-Eulerian: mesh moves independently of material |
| AMR | Adaptive Mesh Refinement: automatic mesh resolution adjustment |
| CFL | Courant-Friedrichs-Lewy: stability condition for explicit schemes |
| FVM | Finite Volume Method: spatial discretization conserving flux |
| GCL | Geometric Conservation Law: constraint on mesh motion |
| h-refinement | Changing mesh topology by splitting/merging cells |
| IMEX | Implicit-Explicit: time integration treating terms differently |
| MLR | Mass Loss Rate: key pyrolysis output [kg/(m²·s)] |
| MMS | Method of Manufactured Solutions: verification technique |
| MOL | Method of Lines: semi-discretization approach |
| r-refinement | Moving mesh nodes without changing topology |
| TVD | Total Variation Diminishing: property of non-oscillatory schemes |

---

**END OF DESIGN DOCUMENT**

*Document version: 2.0*
*Last updated: [Current Date]*
*Target implementation: Julia 1.9+*
