# Analytical Jacobian Performance Optimization Analysis

**Date:** 2026-01-17
**Target File:** [`src/TimeIntegration/analytical_jacobian.jl`](../src/TimeIntegration/analytical_jacobian.jl)
**Problem Size:** 400 cells, 13 components, 5,601 DOFs

---

## Executive Summary

After analyzing the analytical Jacobian implementation, I have identified **six critical bottlenecks** that consume most of the solver time. For a problem with 400 cells and 13 components (5,601 DOFs), the Jacobian computation involves approximately:

- 400 cell iterations for local terms (reactions, properties)
- 399 interior face iterations for conduction, gas transport, and convective heat flux
- Each face loop iterates over 13 components, with nested loops over 13 components for cross-derivatives

This results in roughly **200,000+ derivative computations** per Jacobian evaluation.

**Expected Total Improvement:** 55-75% reduction in Jacobian computation time

---

## Table of Contents

1. [Critical Priority Optimizations (50-60% speedup)](#critical-priority-high-impact-optimizations-50-60-speedup)
   - [1. Eliminate Redundant Property Caching](#1-eliminate-redundant-property-caching-in-face-loops-25-30-improvement)
   - [2. Hoist Common Derivative Computations](#2-hoist-common-derivative-computations-out-of-nested-loops-20-25-improvement)
   - [3. Fuse Duplicate Face Loop Functions](#3-fuse-duplicate-face-loop-functions-25-30-improvement)
2. [Medium Priority Optimizations (10-20% speedup)](#medium-priority-additional-optimizations-10-20-improvement)
   - [4. Batch Sparse Matrix Writes](#4-batch-sparse-matrix-writes-10-15-improvement)
   - [5. Exploit Radiation Absorption Sparsity](#5-exploit-radiation-absorption-sparsity-30-50-when-radiation-active)
   - [6. SIMD Optimization](#6-simd-optimization-for-volume-fraction-computation-2-5-improvement)
3. [Performance Impact Summary](#performance-impact-summary)
4. [Implementation Roadmap](#recommended-implementation-order)
5. [Validation & Benchmarking](#additional-recommendations)

---

## 🔴 Critical Priority: High-Impact Optimizations (50-60% speedup)

### **1. Eliminate Redundant Property Caching in Face Loops (25-30% improvement)**

#### Problem Location
- Lines 1385-1400: `compute_conduction_jacobian!`
- Lines 1470-1540: `compute_gas_transport_jacobian!`
- Lines 1680-1750: `compute_convective_heat_flux_jacobian!`

#### Issue
Each face loop calls `cache_all_component_properties!` twice per face (left and right cells). For 399 faces, this results in **~800 property evaluations**, but each cell is evaluated **multiple times** as it appears in different faces.

**Current Code Pattern:**
```julia
for f in eachindex(mesh.faces)
    # Left cell
    cache_all_component_properties!(base_cache, T_L, material)  # 91 function calls
    # ... computation ...

    # Right cell
    cache_all_component_properties!(base_cache, T_R, material)  # 91 function calls
end
```

For 400 cells: ~800 caching operations with 13 components × 7 properties = **72,618 redundant evaluations**

#### Impact
**Estimated 15-25% of total Jacobian time**

#### Optimization

Add per-cell property caches to the `JacobianCache` structure:

```julia
# Add to JacobianCache structure (after line 178):
comp_densities_per_cell::Matrix{Float64}        # [NC × n_cells]
comp_heat_capacities_per_cell::Matrix{Float64}  # [NC × n_cells]
comp_conductivities_per_cell::Matrix{Float64}   # [NC × n_cells]
comp_lambdas_per_cell::Matrix{Float64}          # [NC × n_cells]
comp_dρ_dT_per_cell::Matrix{Float64}            # [NC × n_cells]
comp_dcp_dT_per_cell::Matrix{Float64}           # [NC × n_cells]
comp_dk_dT_per_cell::Matrix{Float64}            # [NC × n_cells]
```

Create a new pre-computation function:

```julia
@generated function precompute_all_cell_properties!(cache, mesh, material::Material{NC}) where NC
    loop_body = Expr[]
    for j in 1:NC
        push!(loop_body, quote
            comp = material.components[$j]
            @inbounds for i in 1:length(mesh.cells)
                T = mesh.cell_states[i].T
                cache.comp_densities_per_cell[$j, i] = comp.density(T)
                cache.comp_heat_capacities_per_cell[$j, i] = comp.heat_capacity(T)
                cache.comp_conductivities_per_cell[$j, i] = comp.conductivity(T)
                cache.comp_lambdas_per_cell[$j, i] = comp.gas_transfer_coeff(T)
                cache.comp_dρ_dT_per_cell[$j, i] = derivative(comp.density, T)
                cache.comp_dcp_dT_per_cell[$j, i] = derivative(comp.heat_capacity, T)
                cache.comp_dk_dT_per_cell[$j, i] = derivative(comp.conductivity, T)
            end
        end)
    end

    return quote
        $(loop_body...)
        nothing
    end
end
```

Then in face loops, replace caching with lookups:

```julia
# Instead of: cache_all_component_properties!(base_cache, T_L, material)
# Use direct column views (zero-cost abstraction):
densities_L = @view cache.comp_densities_per_cell[:, left_id]
dk_dT_L_cached = @view cache.comp_dk_dT_per_cell[:, left_id]
# etc.
```

#### Benefits
- **Compute properties once per cell** instead of 2× per face
- Better **cache locality** - all properties for a cell stored contiguously
- Enables further optimizations by having all cell data pre-computed

---

### **2. Hoist Common Derivative Computations Out of Nested Loops (20-25% improvement)**

#### Problem Location
- Lines 1558-1615: `compute_gas_transport_jacobian!`
- Lines 1803-1885: `compute_convective_heat_flux_jacobian!`

#### Issue
In the gas transport Jacobian, for each gas component `j`, we loop over ALL components `k` and recompute `dρ_dξ` and `dλ_dξ`:

```julia
for j in 1:NC  # For each gas (typically 6)
    if !is_gas[j]; continue; end
    for k in 1:NC  # For each component (13)
        dρ_L_dξ = compute_dρ_dξ_cached(ξ_L, k, ...)  # RECOMPUTED for every j!
        dλ_L_dξ = compute_dλ_dξ_cached(ξ_L, k, ...)  # RECOMPUTED for every j!
        # ...
    end
end
```

This means: **399 faces × 6 gases × 13 components × 2 cells = 62,244 calls** to these functions, when we only need **399 × 13 × 2 = 10,374 calls**.

#### Impact
**Estimated 20-30% of Jacobian time**

#### Optimization

Add work arrays to `JacobianCache`:

```julia
# Add to JacobianCache (around line 164):
work_NC_5::Vector{Float64}  # dρ_dξ for left cell
work_NC_6::Vector{Float64}  # dρ_dξ for right cell
work_NC_7::Vector{Float64}  # dλ_dξ for left cell
work_NC_8::Vector{Float64}  # dλ_dξ for right cell
```

Restructure the gas transport Jacobian computation:

```julia
# In compute_gas_transport_jacobian! (replace lines 1558-1615):

# Pre-compute derivatives ONCE per face/cell
dρ_L_dξ_array = base_cache.work_NC_5
dλ_L_dξ_array = base_cache.work_NC_7
dρ_R_dξ_array = base_cache.work_NC_6
dλ_R_dξ_array = base_cache.work_NC_8

@inbounds for k in 1:NC
    dρ_L_dξ_array[k] = compute_dρ_dξ_cached(ξ_L, k, densities_L, swelling_factors, is_gas)
    dλ_L_dξ_array[k] = compute_dλ_dξ_cached(ξ_L, k, densities_L, lambdas_L, v_L, V_total_L, mixing, β)
    dρ_R_dξ_array[k] = compute_dρ_dξ_cached(ξ_R, k, densities_R, swelling_factors, is_gas)
    dλ_R_dξ_array[k] = compute_dλ_dξ_cached(ξ_R, k, densities_R, lambdas_R, v_R, V_total_R, mixing, β)
end

# NOW loop over gas components using pre-computed values
for j in 1:NC
    if !is_gas[j]; continue; end

    # ... λρ_face computation ...

    # Left cell derivatives
    for k in 1:NC
        ξ_col = n + (k-1)*n + left_id

        dρ_L_dξ = dρ_L_dξ_array[k]  # O(1) lookup!
        dλ_L_dξ = dλ_L_dξ_array[k]  # O(1) lookup!

        # ... rest of computation ...
    end

    # Right cell derivatives
    for k in 1:NC
        dρ_R_dξ = dρ_R_dξ_array[k]  # O(1) lookup!
        dλ_R_dξ = dλ_R_dξ_array[k]  # O(1) lookup!
        # ... rest of computation ...
    end
end
```

#### Additional Sparsity Optimization

Skip components with zero absorption coefficient:

```julia
# After computing dλ_L_dξ_array, identify contributing components:
n_contributing = 0
contributing_indices = base_cache.work_NR  # Reuse as Int buffer

for k in 1:NC
    if abs(dλ_L_dξ_array[k]) > JAC_ZERO_TOL || abs(dλ_R_dξ_array[k]) > JAC_ZERO_TOL
        n_contributing += 1
        contributing_indices[n_contributing] = k
    end
end

# Then only iterate contributing components:
for idx in 1:n_contributing
    k = Int(contributing_indices[idx])
    # ... computation with k ...
end
```

#### Benefits
- **6× reduction in derivative computations** (from 62,244 to 10,374)
- Enables sparsity exploitation to skip zero-contribution terms
- Same pattern can be applied to convective heat flux Jacobian

---

### **3. Fuse Duplicate Face Loop Functions (25-30% improvement)**

#### Problem Location
Three nearly identical face loops at:
- Lines 1320-1431: `compute_conduction_jacobian!`
- Lines 1453-1620: `compute_gas_transport_jacobian!`
- Lines 1673-1930: `compute_convective_heat_flux_jacobian!`

#### Issue
Each function computes similar quantities (mass fractions, density derivatives, harmonic means) but processes faces **three separate times**. This leads to:
- **Redundant property lookups** (3× per face)
- **Poor cache locality** (data evicted between loops)
- **Repeated derivative computations** (especially `dρ_dξ` used in both gas transport and convective heat)

#### Impact
**Estimated 25-35% of total Jacobian time due to triplication**

#### Optimization Strategy

Create a single fused face loop that computes shared quantities once and applies them to all three physics:

```julia
function compute_all_face_jacobians_fused!(J, smap, mesh, prop_cache, material,
                                            J_gas, n, Val(NC), jac_cache) where NC
    base_cache = jac_cache isa JacobianCacheALE ? jac_cache.base : jac_cache

    # Pre-cache ALL cell properties ONCE
    precompute_all_cell_properties!(base_cache, mesh, material)

    # Get mixing parameters (shared across physics)
    k_mixing = material.conductivity_mixing
    β_k = material.β
    λ_mixing = material.gas_transfer_mixing
    β_λ = material.β_λ
    is_gas = base_cache.is_gas
    swelling = base_cache.swelling_factors

    # Work arrays for derivatives (computed once per face)
    dρ_L_dξ = base_cache.work_NC_5
    dρ_R_dξ = base_cache.work_NC_6
    dλ_L_dξ = base_cache.work_NC_7
    dλ_R_dξ = base_cache.work_NC_8

    @inbounds for f in eachindex(mesh.faces)
        face = mesh.faces[f]
        left_id, right_id = face.cell_ids
        left_id == 0 || right_id == 0 && continue

        # Get cell states
        T_L = mesh.cell_states[left_id].T
        T_R = mesh.cell_states[right_id].T
        ξ_L = mesh.cell_states[left_id].ξ
        ξ_R = mesh.cell_states[right_id].ξ

        # Get pre-cached properties (column views - zero cost)
        densities_L = @view base_cache.comp_densities_per_cell[:, left_id]
        densities_R = @view base_cache.comp_densities_per_cell[:, right_id]
        conductivities_L = @view base_cache.comp_conductivities_per_cell[:, left_id]
        conductivities_R = @view base_cache.comp_conductivities_per_cell[:, right_id]
        dk_dT_L = @view base_cache.comp_dk_dT_per_cell[:, left_id]
        dk_dT_R = @view base_cache.comp_dk_dT_per_cell[:, right_id]
        lambdas_L = @view base_cache.comp_lambdas_per_cell[:, left_id]
        lambdas_R = @view base_cache.comp_lambdas_per_cell[:, right_id]

        # Compute shared geometry
        d_LR = compute_face_distance(face, mesh)
        d_L = mesh.cells[left_id].volume / face.area
        d_R = mesh.cells[right_id].volume / face.area

        # Compute shared mixing quantities
        k_L = conductivity_mixing(conductivities_L, ξ_L, k_mixing, β_k)
        k_R = conductivity_mixing(conductivities_R, ξ_R, k_mixing, β_k)
        k_face = harmonic_mean_weighted(k_L, k_R, d_L, d_R)

        # Compute bulk densities and volume totals (needed for multiple physics)
        V_total_L = compute_volume_total(ξ_L, densities_L)
        V_total_R = compute_volume_total(ξ_R, densities_R)
        ρ_L = compute_bulk_density(ξ_L, V_total_L)
        ρ_R = compute_bulk_density(ξ_R, V_total_R)

        # =================================================================
        # CONDUCTION JACOBIAN (inline, no function call)
        # =================================================================
        ΔT = T_R - T_L
        if abs(ΔT) > JAC_ZERO_TOL
            # Compute derivatives as in original function
            dk_dkL, dk_dkR = harmonic_mean_weighted_derivatives(k_L, k_R, d_L, d_R)

            # Temperature derivatives
            dk_L_dξ = base_cache.work_NC
            dk_R_dξ = base_cache.work_NC_2

            dk_dT_L_val = conductivity_derivatives_cached!(
                dk_L_dξ, ξ_L, densities_L, conductivities_L, dk_dT_L,
                is_gas, swelling, k_mixing, β_k, base_cache.ξ_work)
            dk_dT_R_val = conductivity_derivatives_cached!(
                dk_R_dξ, ξ_R, densities_R, conductivities_R, dk_dT_R,
                is_gas, swelling, k_mixing, β_k, base_cache.ξ_work)

            # Apply to Jacobian
            factor = -1.0 / d_LR
            coeff_L = factor / (prop_cache.ρcp[left_id] * mesh.cells[left_id].volume)
            coeff_R = factor / (prop_cache.ρcp[right_id] * mesh.cells[right_id].volume)

            # ... (rest of conduction Jacobian application)
        end

        # =================================================================
        # GAS TRANSPORT JACOBIAN (inline, reuse computations)
        # =================================================================
        # Pre-compute derivatives ONCE for this face (optimization #2)
        @inbounds for k in 1:NC
            dρ_L_dξ[k] = compute_dρ_dξ_cached(ξ_L, k, densities_L, swelling, is_gas)
            dρ_R_dξ[k] = compute_dρ_dξ_cached(ξ_R, k, densities_R, swelling, is_gas)
            dλ_L_dξ[k] = compute_dλ_dξ_cached(ξ_L, k, densities_L, lambdas_L,
                                              V_total_L, λ_mixing, β_λ)
            dλ_R_dξ[k] = compute_dλ_dξ_cached(ξ_R, k, densities_R, lambdas_R,
                                              V_total_R, λ_mixing, β_λ)
        end

        # Apply to each gas component
        for j in 1:NC
            if !is_gas[j]; continue; end

            # ... (gas transport Jacobian using pre-computed dρ_dξ, dλ_dξ)
        end

        # =================================================================
        # CONVECTIVE HEAT FLUX JACOBIAN (inline, REUSE dρ_dξ, dλ_dξ!)
        # =================================================================
        # The dρ_L_dξ and dλ_L_dξ arrays are already computed above!
        # No need to recompute them - this is a major win!

        for j in 1:NC
            if !is_gas[j]; continue; end

            # ... (convective heat Jacobian using SAME dρ_dξ, dλ_dξ arrays)
        end
    end

    return nothing
end
```

#### Benefits
- **Single property cache pass** instead of 3×
- **Single face iteration** instead of 3×
- **Reuse computed derivatives** across physics (e.g., `dρ_dξ` used in both gas transport and convective heat)
- **Better cache locality** - all face data accessed together, reducing cache misses
- **Easier to maintain** - all face-based physics in one place

#### Implementation Notes
- This is the most complex optimization but has the highest impact
- Requires careful testing to ensure all three physics are correctly implemented
- Consider implementing incrementally: start by fusing two functions, then add the third
- Maintain separate functions initially for testing, then switch to fused version

---

## 🟡 Medium Priority: Additional Optimizations (10-20% improvement)

### **4. Batch Sparse Matrix Writes (10-15% improvement)**

#### Problem Location
Throughout the Jacobian computation (lines 1414, 1416, 1422, 1424, 1583, 1584, etc.)

#### Issue
Every `fast_sparse_add!` call performs:
```julia
@inline function fast_sparse_add!(J, map::FastSparsityMap, row, col, val)
    idx = row + (col - 1) * map.n_rows  # Arithmetic
    nz_idx = @inbounds map.indices[idx]  # Memory access
    if nz_idx > 0                         # Branch
        @inbounds J.nzval[nz_idx] += val  # Memory write
    end
end
```

For ~200,000 Jacobian entries, this overhead (index calculation, map lookup, branch) adds up to **10-15% of Jacobian time**.

#### Optimization

Add batched sparse write functions:

```julia
# Add to sparsity map utilities:
@inline function fast_sparse_add_batch!(J, map::FastSparsityMap, row::Int,
                                        cols::NTuple{N,Int},
                                        vals::NTuple{N,Float64}) where N
    @inbounds for i in 1:N
        idx = row + (cols[i] - 1) * map.n_rows
        nz_idx = map.indices[idx]
        if nz_idx > 0
            J.nzval[nz_idx] += vals[i]
        end
    end
end

# For subtraction:
@inline function fast_sparse_sub_batch!(J, map::FastSparsityMap, row::Int,
                                        cols::NTuple{N,Int},
                                        vals::NTuple{N,Float64}) where N
    @inbounds for i in 1:N
        idx = row + (cols[i] - 1) * map.n_rows
        nz_idx = map.indices[idx]
        if nz_idx > 0
            J.nzval[nz_idx] -= vals[i]
        end
    end
end
```

Usage example (replace lines 1414-1425 in conduction Jacobian):

```julia
# Instead of 4 separate calls:
fast_sparse_sub!(J, smap, left_id, ξ_col_L, coeff_L * dq_dξ_jL)
fast_sparse_add!(J, smap, right_id, ξ_col_L, coeff_R * dq_dξ_jL)
fast_sparse_sub!(J, smap, left_id, ξ_col_R, coeff_L * dq_dξ_jR)
fast_sparse_add!(J, smap, right_id, ξ_col_R, coeff_R * dq_dξ_jR)

# Use two batched calls (one per row):
fast_sparse_sub_batch!(J, smap, left_id,
    (ξ_col_L, ξ_col_R),
    (coeff_L * dq_dξ_jL, coeff_L * dq_dξ_jR))

fast_sparse_add_batch!(J, smap, right_id,
    (ξ_col_L, ξ_col_R),
    (coeff_R * dq_dξ_jL, coeff_R * dq_dξ_jR))
```

#### Benefits
- **Reduced function call overhead** (2 calls instead of 4)
- **Better instruction pipelining** - compiler can optimize the loop
- **Fewer branches** - amortized across batch
- Works especially well when writing multiple entries to the same row

---

### **5. Exploit Radiation Absorption Sparsity (30-50% when radiation active)**

#### Problem Location
Lines 4010-4025: Beer-Lambert radiation Jacobian loop

#### Issue
The radiation Jacobian has an O(n²) loop over all cells:

```julia
for i in 1:n              # 400 iterations
    for k in (i+1):n      # ~200 average
        for j in 1:NC     # 13 iterations
            # Compute derivative
            # Total: 400 × 200 × 13 = 1,040,000 operations
        end
    end
end
```

However, **typically only 1-3 components have non-zero absorption coefficients**. We're iterating over all 13 components when most contribute nothing.

#### Impact
**When radiation is active, this can dominate (30-50% of Jacobian time)**

#### Optimization

Pre-identify absorbing components (do once at initialization):

```julia
# Add to JacobianCache or as a cached property:
struct RadiationSparsityInfo
    absorbing_components::Vector{Int}
    n_absorbing::Int
end

function identify_absorbing_components(material::Material{NC}) where NC
    absorbing = Int[]
    for j in 1:NC
        if material.components[j].absorption_coeff > 1e-12  # JAC_ZERO_TOL
            push!(absorbing, j)
        end
    end
    RadiationSparsityInfo(absorbing, length(absorbing))
end
```

Use in radiation Jacobian:

```julia
# In radiation Jacobian computation:
absorbing = radiation_sparsity.absorbing_components

for i in 1:n
    for k in (i+1):n
        for idx in 1:radiation_sparsity.n_absorbing
            j = absorbing[idx]  # Only iterate non-zero absorption components!
            # Computation
        end
    end
end

# New total: 400 × 200 × 3 = 240,000 operations (76% reduction if 3 absorbing components)
```

#### Benefits
- **Massive reduction** when few components absorb radiation (typical case)
- **Zero overhead** if all components absorb (just indirect indexing)
- Simple to implement with large impact

---

### **6. SIMD Optimization for Volume Fraction Computation (2-5% improvement)**

#### Problem Location
Lines 2274-2288: Volume fraction computation in various derivative functions

#### Current Code
```julia
@inbounds @simd for j in 1:NC
    V_total += ξ[j] / max(densities[j], DENSITY_MIN)
end
```

The `@simd` macro provides hints to the compiler but doesn't guarantee vectorization.

#### Optimization

Use `@turbo` from LoopVectorization.jl for explicit SIMD:

```julia
using LoopVectorization

@turbo for j in 1:NC
    V_total += ξ[j] / max(densities[j], DENSITY_MIN)
end
```

Apply to similar loops throughout the Jacobian computation:

```julia
# Example: density computation
@turbo for j in 1:NC
    ρ_total += ξ[j]
end

# Example: property mixing
@turbo for j in 1:NC
    k_total += v[j] * conductivities[j]
end
```

#### Benefits
- **Guaranteed vectorization** on compatible hardware
- **2-5% improvement** in innermost loops
- **Low implementation cost** - just add package and change macro
- Works best for simple arithmetic loops without complex branching

#### Notes
- LoopVectorization.jl works best with Float64 arrays
- May not vectorize loops with complex control flow
- Test on your target hardware to verify speedup

---

## Performance Impact Summary

| Optimization | Impact | Complexity | Priority | Lines Affected |
|-------------|--------|------------|----------|----------------|
| **1. Pre-cache cell properties** | 25-30% | Medium | 🔴 Critical | 1385-1400, 1470-1540, 1680-1750 |
| **2. Hoist dρ/dξ, dλ/dξ** | 20-25% | Low | 🔴 Critical | 1558-1615, 1803-1885 |
| **3. Fuse face loops** | 25-30% | High | 🔴 Critical | 1320-1931 (entire section) |
| **4. Batch sparse writes** | 10-15% | Medium | 🟡 Medium | Throughout (200+ locations) |
| **5. Radiation sparsity** | 30-50%* | Low | 🟡 Medium | 4010-4025 |
| **6. SIMD in loops** | 2-5% | Low | 🟢 Low | 2274-2288, others |

*When radiation model is active

### Expected Cumulative Improvement

- **Optimizations #1 + #2 alone:** 40-50% faster Jacobian
- **All critical optimizations (#1-#3):** 55-70% faster Jacobian
- **Including medium priority (#1-#5):** 60-75% faster overall
- **All optimizations:** 65-80% faster Jacobian computation

### Impact on Overall Solver Time

From the benchmark results in `comprehensive_solver_comparison.jl`, the analytical Jacobian is called dozens of times per solve. If Jacobian computation represents 60-70% of solver time:

- **40% Jacobian speedup** → 24-28% total solver speedup
- **60% Jacobian speedup** → 36-42% total solver speedup
- **75% Jacobian speedup** → 45-53% total solver speedup

---

## Recommended Implementation Order

### **Quick Wins (1-2 days)**

**Priority:** Get immediate 20-35% improvement with low risk

1. **Implement Optimization #2** (hoist derivatives)
   - **Why first:** Easiest to implement, no structural changes
   - **Expected gain:** 20-25%
   - **Files to modify:** `analytical_jacobian.jl` lines 1558-1615
   - **Testing:** Run existing test suite, verify Jacobian accuracy

2. **Implement Optimization #5** (radiation sparsity)
   - **Why second:** Trivial if radiation is used, independent of other changes
   - **Expected gain:** 30-50% (when radiation active)
   - **Files to modify:** `analytical_jacobian.jl` lines 4010-4025
   - **Testing:** Run radiation test cases

**Milestone:** ~35-45% improvement in typical cases

---

### **Major Refactor (3-5 days)**

**Priority:** Achieve maximum performance with structural improvements

3. **Implement Optimization #1** (pre-cache cell properties)
   - **Why third:** Requires cache structure changes but enables optimization #3
   - **Expected gain:** 25-30%
   - **Files to modify:**
     - `analytical_jacobian.jl` lines 129-230 (JacobianCache struct)
     - Lines 1385-1400, 1470-1540, 1680-1750 (face loops)
   - **Testing:** Comprehensive test suite, check allocations

4. **Implement Optimization #3** (fuse face loops)
   - **Why fourth:** Most complex but highest single-optimization gain
   - **Expected gain:** 25-30%
   - **Files to modify:** `analytical_jacobian.jl` lines 1320-1931
   - **Testing:** Extensive testing, compare results with original functions
   - **Strategy:**
     - Create new fused function alongside existing ones
     - Test extensively with flag to switch between implementations
     - Once validated, remove old functions

**Milestone:** ~55-70% improvement achieved

---

### **Polish (1-2 days)**

**Priority:** Extract final performance with finishing touches

5. **Implement Optimization #4** (batch sparse writes)
   - **Why fifth:** Requires many small changes throughout codebase
   - **Expected gain:** 10-15%
   - **Files to modify:** Create new utility functions, update ~200 call sites
   - **Testing:** Verify sparse matrix structure unchanged

6. **Implement Optimization #6** (SIMD)
   - **Why last:** Smallest gain, easy to add
   - **Expected gain:** 2-5%
   - **Files to modify:** Add LoopVectorization.jl, update inner loops
   - **Testing:** Benchmark on target hardware

**Milestone:** ~60-75% total improvement

---

### Implementation Strategy Notes

#### Test-Driven Development
For each optimization:
1. **Save reference results** before changes
2. **Implement optimization** with feature flag to toggle
3. **Run full test suite** comparing old vs new
4. **Benchmark performance** to verify improvement
5. **Remove old code** only after validation

#### Validation Checklist
- [ ] All unit tests pass
- [ ] Jacobian matches finite difference (within tolerance)
- [ ] Solver converges to same solution
- [ ] No new allocations introduced
- [ ] Performance improves as expected
- [ ] Code is documented and maintainable

#### Rollback Plan
- Keep git commits atomic (one optimization per commit)
- Tag working versions before major refactors
- Maintain feature flags for A/B testing
- Document any behavioral changes

---

## Additional Recommendations

### **Profile to Validate Assumptions**

Before implementing optimizations, confirm hotspots with Julia's profiler:

```julia
using Profile, ProfileView

# Setup problem
prob = ODEProblem(ode_fn_analytical_sparse, u0, (0.0, 1.0))

# Warmup
analytical_jac_ale!(J, u0, nothing, 0.0)

# Profile Jacobian computation
@profile for i in 1:100
    fill!(J.nzval, 0.0)
    analytical_jac_ale!(J, u0, nothing, 0.0)
end

# Visualize flamegraph
ProfileView.view()
```

**Look for:**
- Time spent in `cache_all_component_properties!`
- Time in nested loops (gas transport, convective heat)
- Time in `compute_dρ_dξ_cached`, `compute_dλ_dξ_cached`
- Time in sparse matrix operations
- Any unexpected allocations (red in flamegraph)

---

### **Benchmark Individual Optimizations**

After each optimization, measure impact:

```julia
using BenchmarkTools

# Before optimization
baseline = @benchmark begin
    fill!($J.nzval, 0.0)
    analytical_jac_ale!($J, $u0, nothing, 0.0)
end samples=50

println("Baseline: ", minimum(baseline.times) / 1e6, " ms")

# After optimization
optimized = @benchmark begin
    fill!($J.nzval, 0.0)
    analytical_jac_ale_optimized!($J, $u0, nothing, 0.0)
end samples=50

println("Optimized: ", minimum(optimized.times) / 1e6, " ms")
println("Speedup: ", minimum(baseline.times) / minimum(optimized.times), "x")
```

---

### **Memory Allocation Check**

Monitor allocations (should be near-zero after optimization):

```julia
# Should be ~0 bytes after warmup
allocs = @allocated begin
    fill!(J.nzval, 0.0)
    analytical_jac_ale!(J, u0, nothing, 0.0)
end

println("Allocations: ", allocs, " bytes")

if allocs > 1000
    # Use Cthulhu.jl to find allocation sources
    using Cthulhu
    @descend analytical_jac_ale!(J, u0, nothing, 0.0)
end
```

---

### **Type Stability Check**

Ensure all functions are type-stable:

```julia
using JET

# Check for type instabilities
@report_opt analytical_jac_ale!(J, u0, nothing, 0.0)

# Check for potential runtime errors
@report_call analytical_jac_ale!(J, u0, nothing, 0.0)
```

---

### **Solver-Level Integration Testing**

After Jacobian optimizations, verify solver performance:

```julia
using OrdinaryDiffEq, BenchmarkTools

# Baseline solve time
baseline_solve = @benchmark solve($prob, Rodas5P(linsolve=KLUFactorization());
    abstol=1e-8, reltol=1e-6, save_everystep=false) samples=10

# After Jacobian optimization
optimized_solve = @benchmark solve($prob_optimized, Rodas5P(linsolve=KLUFactorization());
    abstol=1e-8, reltol=1e-6, save_everystep=false) samples=10

println("Solver speedup: ",
    minimum(baseline_solve.times) / minimum(optimized_solve.times), "x")
```

---

### **Continuous Benchmarking**

Set up automated performance tracking:

```julia
# Save in test/benchmarks/jacobian_benchmark.jl
using BenchmarkTools, JSON

function run_jacobian_benchmarks()
    results = Dict()

    # Setup
    # ... (problem setup code)

    # Benchmark Jacobian
    b = @benchmark begin
        fill!($J.nzval, 0.0)
        analytical_jac_ale!($J, $u0, nothing, 0.0)
    end samples=50

    results["jacobian_time_ms"] = minimum(b.times) / 1e6
    results["jacobian_allocs_bytes"] = b.allocs
    results["jacobian_memory_mb"] = b.memory / 1e6

    # Save results
    open("benchmark_results.json", "w") do f
        JSON.print(f, results, 4)
    end

    return results
end
```

Run in CI to track performance regressions.

---

## Additional Potential Optimizations (Future Work)

### **7. Parallel Jacobian Computation**

For very large problems (n > 1000), consider parallelizing independent computations:

```julia
using Base.Threads

# Cell-local terms (reactions, sources) are independent
@threads for i in 1:n_cells
    compute_cell_jacobian_terms!(J, i, ...)
end

# Face terms have dependencies but can use task-based parallelism
@sync for f in eachindex(mesh.faces)
    @spawn compute_face_jacobian_terms!(J, f, ...)
end
```

**Expected impact:** 1.5-3× on multi-core systems for large problems
**Complexity:** High (requires thread-safe sparse matrix writes)

---

### **8. GPU Acceleration**

For extremely large problems (n > 5000), offload Jacobian to GPU:

```julia
using CUDA, CUDA.CUSPARSE

# Transfer state to GPU
u_gpu = CuArray(u)
J_gpu = CuSparseMatrixCSC(J)

# Compute Jacobian on GPU
analytical_jac_gpu!(J_gpu, u_gpu, ...)

# Transfer back
J .= Array(J_gpu)
```

**Expected impact:** 5-10× for very large problems
**Complexity:** Very high (complete GPU kernel rewrite)

---

### **9. Automatic Differentiation Validation**

Use AD to validate analytical Jacobian and find bugs:

```julia
using ForwardDiff, SparseDiffTools

# Finite difference Jacobian
J_fd = similar(J)
auto_jacvec!(J_fd, rhs!, u0, colorvec, jac_cache)

# Analytical Jacobian
J_analytical = similar(J)
analytical_jac_ale!(J_analytical, u0, nothing, 0.0)

# Compare
max_error = maximum(abs.(J_fd .- J_analytical))
println("Max Jacobian error: ", max_error)

if max_error > 1e-6
    # Find problematic entries
    errors = abs.(J_fd .- J_analytical)
    idx = argmax(errors)
    println("Largest error at: ", idx)
end
```

---

### **10. Adaptive Jacobian Recomputation**

For some problems, recompute Jacobian only when state changes significantly:

```julia
mutable struct AdaptiveJacobianCache
    last_u::Vector{Float64}
    last_J::SparseMatrixCSC{Float64}
    recompute_threshold::Float64
end

function analytical_jac_adaptive!(J, u, p, t, cache)
    # Check if state changed significantly
    Δu = norm(u - cache.last_u) / norm(u)

    if Δu > cache.recompute_threshold
        # Recompute Jacobian
        analytical_jac_ale!(J, u, p, t)
        cache.last_u .= u
        cache.last_J .= J
    else
        # Reuse previous Jacobian
        J .= cache.last_J
    end
end
```

**Expected impact:** Variable (problem-dependent)
**Complexity:** Medium (requires tuning threshold)

---

## Conclusion

The analytical Jacobian implementation has significant optimization opportunities. By implementing the **three critical optimizations** (pre-cache properties, hoist derivatives, fuse loops), you can achieve **55-70% reduction in Jacobian computation time**, translating to **35-50% faster overall solves**.

The optimizations are listed in order of decreasing impact and increasing complexity. I recommend starting with **optimization #2** (easiest, 20-25% gain) to build confidence, then proceeding to the more complex optimizations.

### Next Steps

1. **Run profiling** to validate assumptions about hotspots
2. **Implement optimization #2** as a proof of concept
3. **Measure impact** and verify correctness
4. **Proceed with remaining optimizations** based on measured impact
5. **Set up continuous benchmarking** to prevent regressions

Would you like me to implement any of these optimizations? I recommend starting with **#2 (hoist derivatives)** for immediate results.

---

**Document Version:** 1.0
**Last Updated:** 2026-01-17
**Analyzed Code:** [`src/TimeIntegration/analytical_jacobian.jl`](../src/TimeIntegration/analytical_jacobian.jl) (86,032 tokens)
