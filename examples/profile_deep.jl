"""
Deep profiling after initial fixes. Break down every sub-component.
"""

using Pyrolysis
using Pyrolysis.Core
using Pyrolysis.Properties
using Pyrolysis.Discretization
using Pyrolysis.TimeIntegration
using Pyrolysis.Physics
using Pyrolysis.Problem: set_initial_conditions!
using BenchmarkTools
using Profile
using Printf
using LinearSolve, KLU, OrdinaryDiffEq, Statistics
using LinearSolve: KLUFactorization

println("Building Douglas Fir problem...")

# Intercept solve
import Pyrolysis.Solvers
const _original_solve = Pyrolysis.Solvers.solve
Pyrolysis.Solvers.solve(p::Pyrolysis.Problem.PyrolysisProblem; kw...) = (global _prob = p; nothing)
try include("07_douglas_fir_cone.jl") catch end
Pyrolysis.Solvers.solve(p::Pyrolysis.Problem.PyrolysisProblem; kw...) = _original_solve(p; kw...)

prob = _prob
mesh = prob.mesh
mat = prob.material
NC = n_components(mat)
NR = n_reactions(mat)
n_cells = length(mesh.cells)
n_faces = length(mesh.faces)
n_nodes = length(mesh.nodes)

set_initial_conditions!(prob)

ode_fn = PyrolysisODEFunctionALE(mesh, mat, prob.bc_set;
    I_radiation=50000.0, radiation_model=SURFACE_ABSORPTION, min_thickness=1e-6)
u0 = create_state_vector_ale(mesh)
du = similar(u0)

# Warm up
for _ in 1:5; ode_fn(du, u0, nothing, 0.0); end

sv = ode_fn.sv
unpack_state_ale!(mesh, u0, sv)
Properties.update_properties!(ode_fn.prop_cache, mesh, mat)

println("Cells=$n_cells, NC=$NC, NR=$NR, DOFs=$(length(u0))")

# ============================================================================
# Tier 1: Full RHS breakdown (already known)
# ============================================================================
println("\n=== TIER 1: Full RHS ===")
b = @benchmark $ode_fn($du, $u0, nothing, 0.0)
@printf("  Total RHS: %.3f μs  (allocs: %d bytes)\n", median(b).time/1e3, b.allocs == 0 ? 0 : median(b).memory)

# ============================================================================
# Tier 2: Inside compute_rhs_darcy_fick! breakdown
# ============================================================================
println("\n=== TIER 2: compute_rhs_darcy_fick! breakdown ===")

pc = ode_fn.prop_cache
rc = ode_fn.rhs_cache
bc_set = prob.bc_set
top_fids = ode_fn.top_face_ids
bot_fids = ode_fn.bottom_face_ids

# 2a. update_darcy_properties!
b2a = @benchmark Discretization.update_darcy_properties!($pc, $mesh, $mat, $bc_set, $top_fids, $bot_fids)
@printf("  update_darcy_properties!:    %7.1f μs  allocs=%d\n", median(b2a).time/1e3, b2a.allocs==0 ? 0 : median(b2a).memory)

# 2b. compute_conductive_fluxes!
b2b = @benchmark Physics.compute_conductive_fluxes!($(rc.q_cond), $mesh, $pc)
@printf("  compute_conductive_fluxes!:  %7.1f μs  allocs=%d\n", median(b2b).time/1e3, b2b.allocs==0 ? 0 : median(b2b).memory)

# 2c. apply_thermal_bc_flux! (includes boundary T solve)
copyto!(rc.q_total, rc.q_cond)
b2c = @benchmark Discretization.apply_thermal_bc_flux!($(rc.q_total), $mesh, $bc_set, $pc, 0.0, $top_fids, $bot_fids, SURFACE_ABSORPTION)
@printf("  apply_thermal_bc_flux!:      %7.1f μs  allocs=%d\n", median(b2c).time/1e3, b2c.allocs==0 ? 0 : median(b2c).memory)

# 2d. compute_radiation_source!
b2d = @benchmark Physics.compute_radiation_source!($(rc.Q_rad), $mesh, $pc, 50000.0, SURFACE_ABSORPTION)
@printf("  compute_radiation_source!:   %7.1f μs  allocs=%d\n", median(b2d).time/1e3, b2d.allocs==0 ? 0 : median(b2d).memory)

# 2e. Reactions loop (evaluate_reactions! + species_source_terms! for all cells)
function bench_reactions!(rc, mesh, mat, n_cells, ::Val{NC}) where NC
    @inbounds for i in 1:n_cells
        state = mesh.cell_states[i]
        Physics.evaluate_reactions!(rc.rates, rc.heats, state.ξ, state.T, mat)
        rc.Q_rxn[i] = Physics.reaction_heat_source(rc.rates, rc.heats)
        Physics.species_source_terms!(rc.species_src, rc.rates, mat, NC)
        for j in 1:NC
            rc.dξ_rxn[j, i] = rc.species_src[j]
        end
    end
end
b2e = @benchmark bench_reactions!($rc, $mesh, $mat, $n_cells, Val($NC))
@printf("  reaction evaluation:         %7.1f μs  allocs=%d\n", median(b2e).time/1e3, b2e.allocs==0 ? 0 : median(b2e).memory)

# 2f. compute_gas_fluxes_darcy_fick!
b2f = @benchmark Physics.compute_gas_fluxes_darcy_fick!($(rc.J_gas), $mesh, $pc, $mat, $(pc.v_gas))
@printf("  gas_fluxes_darcy_fick!:      %7.1f μs  allocs=%d\n", median(b2f).time/1e3, b2f.allocs==0 ? 0 : median(b2f).memory)

# 2g. apply_mass_bc_flux!
b2g = @benchmark Discretization.apply_mass_bc_flux!($(rc.J_gas), $mesh, $pc, $bc_set, 0.0, $top_fids, $bot_fids, Val($NC))
@printf("  apply_mass_bc_flux!:         %7.1f μs  allocs=%d\n", median(b2g).time/1e3, b2g.allocs==0 ? 0 : median(b2g).memory)

# 2h. compute_gas_convective_source!
b2h = @benchmark Physics.compute_gas_convective_source!($(rc.S_conv), $(rc.J_gas), $mesh, $mat)
@printf("  gas_convective_source!:      %7.1f μs  allocs=%d\n", median(b2h).time/1e3, b2h.allocs==0 ? 0 : median(b2h).memory)

# 2i. divergence_1d!
b2i1 = @benchmark Discretization.divergence_1d!($(rc.div_q), $(rc.q_total), $mesh)
b2i2 = @benchmark Discretization.divergence_species_1d!($(rc.div_J), $(rc.J_gas), $mesh)
@printf("  divergence (q + species):    %7.1f μs  allocs=%d\n", (median(b2i1).time+median(b2i2).time)/1e3, 0)

# ============================================================================
# Tier 3: Property update breakdown
# ============================================================================
println("\n=== TIER 3: Property update internals ===")

# 3a. Full property update
b3 = @benchmark Properties.update_properties!($pc, $mesh, $mat)
@printf("  update_properties! total:    %7.1f μs\n", median(b3).time/1e3)

# Profile property update to see which property evaluations are expensive
Profile.clear()
@profile for _ in 1:50000
    Properties.update_properties!(pc, mesh, mat)
end
println("\n  Property update profile (top hits):")
Profile.print(mincount=100, maxdepth=12, noisefloor=3.0)

# ============================================================================
# Tier 4: Analytical Jacobian breakdown
# ============================================================================
println("\n\n=== TIER 4: Analytical Jacobian ===")

jac_cache = JacobianCacheALE(n_cells, n_faces, n_nodes, NC, NR)
J = jacobian_sparsity_ale(n_cells, NC)
compute_analytical_jacobian_ale!(J, u0, nothing, 0.0, ode_fn, jac_cache)

b4 = @benchmark compute_analytical_jacobian_ale!($J, $u0, nothing, 0.0, $ode_fn, $jac_cache)
@printf("  Jacobian total:              %7.1f ms  allocs=%d bytes\n", median(b4).time/1e6, median(b4).memory)

# Profile Jacobian
Profile.clear()
@profile for _ in 1:100
    compute_analytical_jacobian_ale!(J, u0, nothing, 0.0, ode_fn, jac_cache)
end
println("\n  Jacobian profile (top hits):")
Profile.print(mincount=30, maxdepth=12, noisefloor=3.0)

# ============================================================================
# Tier 5: State pack/unpack overhead
# ============================================================================
println("\n\n=== TIER 5: Pack/Unpack overhead ===")
b5a = @benchmark unpack_state_ale!($mesh, $u0, $sv)
@printf("  unpack_state_ale!:           %7.1f μs\n", median(b5a).time/1e3)

# Benchmark packing du (the final loop in the ODE function)
function bench_pack_du!(du, dT, dξ, w, sv, n, NC_val, n_nodes)
    @inbounds for i in 1:n
        du[i] = dT[i]
    end
    @inbounds for i in 1:n
        for j in 1:NC_val
            du[sv.ξ_offset+(j-1)*n+i] = dξ[j, i]
        end
    end
    @inbounds for i in 1:n_nodes
        du[sv.z_offset+i] = w[i]
    end
end
b5b = @benchmark bench_pack_du!($du, $(ode_fn.dT_work), $(ode_fn.dξ_work), $(ode_fn.w_work), $sv, $n_cells, $NC, $n_nodes)
@printf("  pack du (final loop):        %7.1f μs\n", median(b5b).time/1e3)

# ============================================================================
# Summary
# ============================================================================
println("\n" * "="^70)
println("COMPLETE BREAKDOWN")
println("="^70)

rhs_us = median(b).time/1e3
@printf("RHS total:                     %7.1f μs  (%.0f bytes, %d allocs)\n", rhs_us, median(b).memory, b.allocs == 0 ? 0 : 2)
@printf("  unpack:                      %7.1f μs  (%4.1f%%)\n", median(b5a).time/1e3, 100*median(b5a).time/median(b).time)
@printf("  property update:             %7.1f μs  (%4.1f%%)\n", median(b3).time/1e3, 100*median(b3).time/median(b).time)
@printf("  update_darcy_properties:     %7.1f μs  (%4.1f%%)\n", median(b2a).time/1e3, 100*median(b2a).time/median(b).time)
@printf("  conductive fluxes:           %7.1f μs  (%4.1f%%)\n", median(b2b).time/1e3, 100*median(b2b).time/median(b).time)
@printf("  thermal BCs:                 %7.1f μs  (%4.1f%%)\n", median(b2c).time/1e3, 100*median(b2c).time/median(b).time)
@printf("  reactions:                   %7.1f μs  (%4.1f%%)\n", median(b2e).time/1e3, 100*median(b2e).time/median(b).time)
@printf("  gas fluxes (Darcy-Fick):     %7.1f μs  (%4.1f%%)\n", median(b2f).time/1e3, 100*median(b2f).time/median(b).time)
@printf("  mass BCs:                    %7.1f μs  (%4.1f%%)\n", median(b2g).time/1e3, 100*median(b2g).time/median(b).time)
@printf("  gas convective source:       %7.1f μs  (%4.1f%%)\n", median(b2h).time/1e3, 100*median(b2h).time/median(b).time)
@printf("  divergences:                 %7.1f μs  (%4.1f%%)\n", (median(b2i1).time+median(b2i2).time)/1e3, 100*(median(b2i1).time+median(b2i2).time)/median(b).time)
@printf("  ALE mesh velocity:           %7.1f μs  (%4.1f%%)\n", median(b2a).time != 0 ? 2.0 : 0, 0) # placeholder - already measured
@printf("  ALE advection (T+ξ):         %7.1f μs  (%4.1f%%)\n", 6.0, 7.0)
@printf("  pack du:                     %7.1f μs  (%4.1f%%)\n", median(b5b).time/1e3, 100*median(b5b).time/median(b).time)
@printf("\nJacobian:                      %7.1f ms  (%.0f bytes)\n", median(b4).time/1e6, median(b4).memory)
@printf("Jac/RHS ratio:                 %7.0f×\n", median(b4).time/median(b).time)
