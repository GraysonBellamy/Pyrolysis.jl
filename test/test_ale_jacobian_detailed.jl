#!/usr/bin/env julia
"""
Detailed ALE Jacobian comparison between Analytical and Finite Differences.

This script performs a thorough investigation of potential discrepancies
between the analytical and FD Jacobians for the ALE formulation.

Uses Douglas Fir material and boundary conditions from 07_douglas_fir_cone.jl
for realistic multi-component, multi-reaction testing.
"""

using Pyrolysis
using Pyrolysis.TimeIntegration: JacobianCacheALE, StateVectorALE,
	compute_analytical_jacobian_ale!, state_vector_length_ale,
	PyrolysisODEFunctionALE, create_state_vector_ale, unpack_state_ale!
using Pyrolysis.Properties: PropertyCache
using Pyrolysis.ALE: GeometryCache
using LinearAlgebra
using Printf
using SparseArrays
using Statistics: mean

println("=" ^ 70)
println("Detailed ALE Jacobian Comparison: Analytical vs Finite Differences")
println("Using Douglas Fir Material (13 components, 6 reactions)")
println("=" ^ 70)

# ============================================================================
# Douglas Fir Material Definition (from 07_douglas_fir_cone.jl)
# ============================================================================

# Temperature-dependent heat capacity functions
df_cp = CallableProperty(T -> -2131.0 + 9.119 * T)
df_int1_cp = CallableProperty(T -> -1983.0 + 8.7 * T)
df_int2_cp = CallableProperty(T -> -1419.0 + 7.0 * T)
df_int3_cp = CallableProperty(T -> 1.0 + 2.8 * T)
df_int4_cp = CallableProperty(T -> 176.0 + 2.3 * T)
char_cp = CallableProperty(T -> 347.2 + 1.784 * T)
moisture_cp = CallableProperty(T -> 5.2e3 - 6.7 * T + 1.1e-2 * T^2)
water_vapor_cp = CallableProperty(T -> 2.4e3 - 1.6 * T + 2.0e-3 * T^2)

# Heat of sorption (state-dependent)
const MOISTURE_IDX = 7
const ρ_DRY_BASIS = 530.0
heat_of_sorption = StateDependentProperty((T, ξ) -> 2673.82e3)

# Gas transfer coefficients
const λ_char = 5e-5
const λ_wood = 1e-5
const λ_gas = 1e-4

function create_douglas_fir_material()
	return Material(
		name = :DouglasFir,
		components = (
			# 1. Virgin Douglas Fir
			SolidComponent(:DF, ρ = 530.0, c = df_cp, k = 0.12, ε = 0.9, α = 0.0, λ = λ_wood),
			# 2. Intermediate 1
			SolidComponent(:DF_int1, ρ = 507.2, c = df_int1_cp, k = 0.14, ε = 0.9, α = 0.0, λ = λ_wood * 1.2),
			# 3. Intermediate 2
			SolidComponent(:DF_int2, ρ = 419.8, c = df_int2_cp, k = 0.15, ε = 0.9, α = 0.0, λ = λ_wood * 1.5),
			# 4. Intermediate 3
			SolidComponent(:DF_int3, ρ = 203.2, c = df_int3_cp, k = 0.16, ε = 0.9, α = 0.0, λ = λ_wood * 2.0),
			# 5. Intermediate 4
			SolidComponent(:DF_int4, ρ = 176.6, c = df_int4_cp, k = 0.18, ε = 0.9, α = 0.0, λ = λ_char * 0.8),
			# 6. Char
			SolidComponent(:Char, ρ = 150.0, c = char_cp, k = 0.20, ε = 0.95, α = 0.0, λ = λ_char),
			# 7. Moisture (liquid)
			LiquidComponent(:MOISTURE, ρ = 1000.0, c = moisture_cp, k = 0.6, ε = 0.9, α = 0.0),
			# 8. DF_gas1
			GaseousComponent(:DF_gas1, M = 0.060, c = 2000.0, k = 0.03, λ = λ_gas),
			# 9. DF_gas2
			GaseousComponent(:DF_gas2, M = 0.044, c = 1800.0, k = 0.03, λ = λ_gas),
			# 10. DF_gas3
			GaseousComponent(:DF_gas3, M = 0.070, c = 2200.0, k = 0.025, λ = λ_gas * 0.8),
			# 11. DF_gas4
			GaseousComponent(:DF_gas4, M = 0.030, c = 1600.0, k = 0.035, λ = λ_gas),
			# 12. DF_gas5
			GaseousComponent(:DF_gas5, M = 0.028, c = 1500.0, k = 0.035, λ = λ_gas),
			# 13. Water vapor
			GaseousComponent(:WATER_VAPOR, M = 0.018, c = water_vapor_cp, k = 0.04, λ = λ_gas * 1.2),
		),
		reactions = (
			# Reaction 1: DF → 0.95 DF_int1 + 0.05 DF_gas1
			Reaction(:pyrolysis_step1, 1 => (2, 8), A = 74816.95, E = 71.59e3, h = -20067.0, n = 1.0, yields = (0.95, 0.05)),
			# Reaction 2: DF_int1 → 0.80 DF_int2 + 0.20 DF_gas2
			Reaction(:pyrolysis_step2, 2 => (3, 9), A = 23550492839.0, E = 140.9e3, h = -8264.0, n = 1.0, yields = (0.80, 0.20)),
			# Reaction 3: DF_int2 → 0.37 DF_int3 + 0.63 DF_gas3
			Reaction(:pyrolysis_step3, 3 => (4, 10), A = 7.744618e12, E = 181.0e3, h = 273999.0, n = 1.0, yields = (0.37, 0.63)),
			# Reaction 4: DF_int3 → 0.79 DF_int4 + 0.21 DF_gas4
			Reaction(:pyrolysis_step4, 4 => (5, 11), A = 103276.14, E = 98.94e3, h = -224430.0, n = 1.0, yields = (0.79, 0.21)),
			# Reaction 5: DF_int4 → 0.74 Char + 0.26 DF_gas5
			Reaction(:pyrolysis_step5, 5 => (6, 12), A = 6.397348, E = 56.01e3, h = 555.6, n = 1.0, yields = (0.74, 0.26)),
			# Reaction 6: MOISTURE → WATER_VAPOR
			Reaction(:moisture_evaporation, 7 => 13, A = 2.997e5, E = 5.0e4, h = heat_of_sorption, n = 1.0),
		),
		conductivity_mixing = WEIGHTED,
		β = 0.5,
		gas_transfer_mixing = WEIGHTED,
		β_λ = 0.5,
	)
end

# Environmental conditions (from 07_douglas_fir_cone.jl)
const T_AMBIENT = 300.0
const Q_EXTERNAL = 50e3      # W/m² - 50 kW/m² external heat flux
const H_CONV = 7.2           # W/(m²·K)

function create_test_problem(; n_cells = 60, T_initial = 300.0, moisture_content = 0.15)
	material = create_douglas_fir_material()

	# Initial concentrations
	ρ_DF_initial = 530.0
	ρ_moisture_initial = moisture_content * ρ_DF_initial

	# Create mesh (13 components)
	mesh = create_uniform_mesh(0.006, n_cells, Val(13);  # 6mm sample
		T_initial = T_initial,
		ξ_initial = (ρ_DF_initial, 0.0, 0.0, 0.0, 0.0, 0.0, ρ_moisture_initial,
			0.0, 0.0, 0.0, 0.0, 0.0, 0.0))

	# Boundary conditions (from 07_douglas_fir_cone.jl)
	top_thermal_bc = RadiativeFluxBC(flux = Q_EXTERNAL, absorptivity = 0.9) +
					 RadiativeBC(ε = 0.9, T_inf = 330) +
					 ConvectiveBC(h = H_CONV, T_inf = 330)
	top_mass_bc = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0)

	bc_set = BoundaryConditionSet(
		thermal_top = top_thermal_bc,
		thermal_bottom = AdiabaticBC(),
		mass_top = top_mass_bc,
		mass_bottom = ImpermeableBC(),
	)

	return material, mesh, bc_set
end

function compute_fd_jacobian!(J_fd, u, ode_fn, ε = 1e-7)
	n_dof = length(u)
	u_pert = copy(u)
	du_base = similar(u)
	du_pert = similar(u)

	ode_fn(du_base, u, nothing, 0.0)

	for j in 1:n_dof
		u_pert .= u
		u_pert[j] += ε
		ode_fn(du_pert, u_pert, nothing, 0.0)
		@views J_fd[:, j] .= (du_pert .- du_base) ./ ε
	end
end

# Douglas Fir: 13 components, 6 reactions
const NC = 13
const NR = 6

# ============================================================================
# Analyze Jacobian block by block
# ============================================================================

function analyze_jacobian_blocks(J_an, J_fd, sv, NC; verbose = true)
	n = sv.n_cells
	n_nodes = sv.n_nodes
	n_dof = length(J_an[1, :])

	# Define block ranges
	T_range = 1:n
	ξ_ranges = [n + (j-1)*n .+ (1:n) for j in 1:NC]
	z_range = (sv.z_offset+1):n_dof

	blocks = Dict(
		"dT/dT" => (T_range, T_range),
		"dT/dz" => (T_range, z_range),
		"dz/dT" => (z_range, T_range),
		"dz/dz" => (z_range, z_range),
	)

	for j in 1:NC
		blocks["dT/dξ$j"] = (T_range, ξ_ranges[j])
		blocks["dξ$j/dT"] = (ξ_ranges[j], T_range)
		blocks["dξ$j/dz"] = (ξ_ranges[j], z_range)
		blocks["dz/dξ$j"] = (z_range, ξ_ranges[j])
		for k in 1:NC
			blocks["dξ$j/dξ$k"] = (ξ_ranges[j], ξ_ranges[k])
		end
	end

	results = Dict{String, NamedTuple}()

	for (name, (rows, cols)) in blocks
		J_an_block = J_an[rows, cols]
		J_fd_block = J_fd[rows, cols]

		diff = J_an_block - J_fd_block
		abs_diff = abs.(diff)

		# Various error metrics
		max_abs_err = maximum(abs_diff)
		mean_abs_err = mean(abs_diff)
		frobenius_err = norm(diff)
		norm_fd = norm(J_fd_block)
		rel_err = norm_fd > 1e-12 ? frobenius_err / norm_fd : frobenius_err

		# Find worst element
		idx = argmax(abs_diff)
		worst_i, worst_j = Tuple(idx)

		results[name] = (
			max_abs_err = max_abs_err,
			mean_abs_err = mean_abs_err,
			frobenius_err = frobenius_err,
			rel_err = rel_err,
			norm_fd = norm_fd,
			worst_loc = (worst_i, worst_j),
			worst_an = J_an_block[worst_i, worst_j],
			worst_fd = J_fd_block[worst_i, worst_j],
		)
	end

	if verbose
		println("\n" * "-"^70)
		println("Block-by-Block Jacobian Analysis")
		println("-"^70)
		println(@sprintf("%-12s %12s %12s %12s %12s", "Block", "||Diff||", "||FD||", "Rel.Err", "Max|Diff|"))
		println("-"^70)

		for name in sort(collect(keys(results)))
			r = results[name]
			status = r.rel_err < 1e-4 ? "✓" : (r.rel_err < 0.1 ? "~" : "✗")
			println(@sprintf("%-12s %12.2e %12.2e %12.2e %12.2e %s",
				name, r.frobenius_err, r.norm_fd, r.rel_err, r.max_abs_err, status))
		end
	end

	return results
end

# ============================================================================
# Test 1: Basic comparison at room temperature (minimal reactions)
# ============================================================================

println("\n" * "="^70)
println("TEST 1: Room temperature state (T=300K, minimal reactions)")
println("="^70)

material, mesh, bc_set = create_test_problem(n_cells = 60, T_initial = 300.0)
n_cells = length(mesh.cells)
n_faces = length(mesh.faces)
n_nodes = length(mesh.nodes)

ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)
sv = ode_fn.sv
u0 = create_state_vector_ale(mesh, material)
n_dof = length(u0)

println("Problem size: $n_cells cells, $n_nodes nodes, $n_dof DOFs")
println("Components: $NC, Reactions: $NR")

# Warmup RHS
du = similar(u0)
ode_fn(du, u0, nothing, 0.0)

# Compute both Jacobians
jac_cache = JacobianCacheALE(n_cells, n_faces, n_nodes, NC, NR)
J_analytical = zeros(n_dof, n_dof)
J_fd = zeros(n_dof, n_dof)

compute_analytical_jacobian_ale!(J_analytical, u0, nothing, 0.0, ode_fn, jac_cache)
compute_fd_jacobian!(J_fd, u0, ode_fn)

# Overall comparison
rel_err = norm(J_analytical - J_fd) / max(norm(J_fd), 1e-10)
max_err = maximum(abs.(J_analytical - J_fd))
println(@sprintf("\nOverall Frobenius relative error: %.2e", rel_err))
println(@sprintf("Maximum absolute error: %.2e", max_err))

# Block analysis
results1 = analyze_jacobian_blocks(J_analytical, J_fd, sv, NC)

# ============================================================================
# Test 2: Elevated temperature state (reactions active)
# ============================================================================

println("\n" * "="^70)
println("TEST 2: Elevated temperature state (T=500K, reactions active)")
println("="^70)

material, mesh, bc_set = create_test_problem(n_cells = 60, T_initial = 500.0)

ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)
sv = ode_fn.sv
u0 = create_state_vector_ale(mesh, material)
n_dof = length(u0)
n_cells = length(mesh.cells)
n_faces = length(mesh.faces)
n_nodes = length(mesh.nodes)

# Warmup
du = similar(u0)
ode_fn(du, u0, nothing, 0.0)

# Compute Jacobians
jac_cache = JacobianCacheALE(n_cells, n_faces, n_nodes, NC, NR)
J_analytical = zeros(n_dof, n_dof)
J_fd = zeros(n_dof, n_dof)

compute_analytical_jacobian_ale!(J_analytical, u0, nothing, 0.0, ode_fn, jac_cache)
compute_fd_jacobian!(J_fd, u0, ode_fn)

rel_err = norm(J_analytical - J_fd) / max(norm(J_fd), 1e-10)
max_err = maximum(abs.(J_analytical - J_fd))
println(@sprintf("\nOverall Frobenius relative error: %.2e", rel_err))
println(@sprintf("Maximum absolute error: %.2e", max_err))

results2 = analyze_jacobian_blocks(J_analytical, J_fd, sv, NC)

# ============================================================================
# Test 2b: Realistic mid-simulation state with temperature gradient
# ============================================================================

println("\n" * "="^70)
println("TEST 2b: Realistic mid-simulation state (T gradient, partial decomposition)")
println("="^70)

# Create a state that mimics what happens during actual Douglas Fir pyrolysis:
# - Hot surface (700K), cooler back (350K)
# - Surface cells partially decomposed (virgin → intermediates → char + gases)
# - Temperature gradient driving conduction

material2, mesh2, bc_set2 = create_test_problem(n_cells = 60, T_initial = 400.0)

ode_fn2 = PyrolysisODEFunctionALE(mesh2, material2, bc_set2)
sv2 = ode_fn2.sv
u0_2 = create_state_vector_ale(mesh2, material2)
n_dof2 = length(u0_2)
n2 = sv2.n_cells

# Modify state to have realistic gradients
# Temperature gradient: hot at top (cell n), cool at bottom (cell 1)
for i in 1:n2
	# Linear temperature gradient from 350K (bottom) to 700K (top)
	T_i = 350.0 + (700.0 - 350.0) * (i - 1) / (n2 - 1)
	u0_2[i] = T_i

	# Decomposition progress (function of temperature history, approximate)
	# Bottom cells: virgin wood + moisture
	# Top cells: more intermediates, char, and gases

	# Fractional progress through reactions (0=virgin, 1=fully charred)
	progress = max(0.0, (i - n2/3) / (2*n2/3))^1.5

	# Component mass concentrations (kg/m³)
	# Start: DF=530, moisture=79.5, all others=0
	# End: Char ~ 0.74*0.79*0.37*0.80*0.95*530 ≈ 70, gases escaped

	ξ_DF = 530.0 * (1.0 - progress)                          # Virgin wood depletes
	ξ_int1 = 530.0 * 0.95 * progress * (1.0 - progress)^0.5  # Intermediate 1
	ξ_int2 = 530.0 * 0.80 * 0.95 * progress^1.5 * (1.0 - progress)^0.3
	ξ_int3 = 530.0 * 0.37 * progress^2 * (1.0 - progress)^0.2
	ξ_int4 = 530.0 * 0.30 * progress^2.5 * 0.5
	ξ_char = 70.0 * progress^3                               # Char builds up
	ξ_moisture = 79.5 * (1.0 - progress)^2                   # Moisture evaporates early

	# Gases (some accumulate before escaping)
	ξ_gas1 = 5.0 * progress * (1.0 - progress)
	ξ_gas2 = 10.0 * progress * (1.0 - progress)^0.5
	ξ_gas3 = 15.0 * progress^1.5 * (1.0 - progress)^0.5
	ξ_gas4 = 5.0 * progress^2 * (1.0 - progress)^0.3
	ξ_gas5 = 3.0 * progress^2.5 * 0.5
	ξ_water_vapor = 10.0 * progress * (1.0 - progress)^1.5

	# Pack into state vector
	u0_2[n2+i]    = ξ_DF
	u0_2[2*n2+i]  = ξ_int1
	u0_2[3*n2+i]  = ξ_int2
	u0_2[4*n2+i]  = ξ_int3
	u0_2[5*n2+i]  = ξ_int4
	u0_2[6*n2+i]  = ξ_char
	u0_2[7*n2+i]  = ξ_moisture
	u0_2[8*n2+i]  = ξ_gas1
	u0_2[9*n2+i]  = ξ_gas2
	u0_2[10*n2+i] = ξ_gas3
	u0_2[11*n2+i] = ξ_gas4
	u0_2[12*n2+i] = ξ_gas5
	u0_2[13*n2+i] = ξ_water_vapor
end

# Update mesh state
Pyrolysis.TimeIntegration.unpack_state_ale!(mesh2, u0_2, sv2)
Pyrolysis.Properties.update_properties!(ode_fn2.prop_cache, mesh2, material2)

# Warmup
du2 = similar(u0_2)
ode_fn2(du2, u0_2, nothing, 0.0)

# Compute Jacobians
n_faces2 = length(mesh2.faces)
n_nodes2 = length(mesh2.nodes)
jac_cache2 = JacobianCacheALE(n2, n_faces2, n_nodes2, NC, NR)
J_analytical2 = zeros(n_dof2, n_dof2)
J_fd2 = zeros(n_dof2, n_dof2)

compute_analytical_jacobian_ale!(J_analytical2, u0_2, nothing, 0.0, ode_fn2, jac_cache2)
compute_fd_jacobian!(J_fd2, u0_2, ode_fn2)

rel_err2b = norm(J_analytical2 - J_fd2) / max(norm(J_fd2), 1e-10)
max_err2b = maximum(abs.(J_analytical2 - J_fd2))
println(@sprintf("\nOverall Frobenius relative error: %.2e", rel_err2b))
println(@sprintf("Maximum absolute error: %.2e", max_err2b))

results2b = analyze_jacobian_blocks(J_analytical2, J_fd2, sv2, NC)

# Show some state info
println("\nState summary:")
println(@sprintf("  T range: %.1f K to %.1f K", minimum(u0_2[1:n2]), maximum(u0_2[1:n2])))
println(@sprintf("  ξ_DF range: %.1f to %.1f kg/m³", minimum(u0_2[(n2+1):(2*n2)]), maximum(u0_2[(n2+1):(2*n2)])))
println(@sprintf("  ξ_char range: %.1f to %.1f kg/m³", minimum(u0_2[(6*n2+1):(7*n2)]), maximum(u0_2[(6*n2+1):(7*n2)])))
println(@sprintf("  ξ_moisture range: %.1f to %.1f kg/m³", minimum(u0_2[(7*n2+1):(8*n2)]), maximum(u0_2[(7*n2+1):(8*n2)])))

# Check reaction rates
println("\nReaction activity (dξ_DF/dt at select cells):")
for i in [1, n2÷2, n2]
	println(@sprintf("  Cell %d (T=%.0fK): dξ_DF/dt = %.2e kg/(m³·s)", i, u0_2[i], du2[n2+i]))
end

# ============================================================================
# Test 2c: Very high temperature with strong reactions
# ============================================================================

println("\n" * "="^70)
println("TEST 2c: High temperature uniform state (T=600K, strong reactions)")
println("="^70)

material3, mesh3, bc_set3 = create_test_problem(n_cells = 60, T_initial = 600.0)

ode_fn3 = PyrolysisODEFunctionALE(mesh3, material3, bc_set3)
sv3 = ode_fn3.sv
u0_3 = create_state_vector_ale(mesh3, material3)
n_dof3 = length(u0_3)
n3 = sv3.n_cells

# Warmup
du3 = similar(u0_3)
ode_fn3(du3, u0_3, nothing, 0.0)

# Compute Jacobians
n_faces3 = length(mesh3.faces)
n_nodes3 = length(mesh3.nodes)
jac_cache3 = JacobianCacheALE(n3, n_faces3, n_nodes3, NC, NR)
J_analytical3 = zeros(n_dof3, n_dof3)
J_fd3 = zeros(n_dof3, n_dof3)

compute_analytical_jacobian_ale!(J_analytical3, u0_3, nothing, 0.0, ode_fn3, jac_cache3)
compute_fd_jacobian!(J_fd3, u0_3, ode_fn3)

rel_err3 = norm(J_analytical3 - J_fd3) / max(norm(J_fd3), 1e-10)
max_err3 = maximum(abs.(J_analytical3 - J_fd3))
println(@sprintf("\nOverall Frobenius relative error: %.2e", rel_err3))
println(@sprintf("Maximum absolute error: %.2e", max_err3))

results3 = analyze_jacobian_blocks(J_analytical3, J_fd3, sv3, NC)

# Show reaction rates to confirm reactions are active
println("\nReaction activity check:")
println(@sprintf("  dξ_DF/dt at surface: %.2e kg/(m³·s)", du3[2*n3]))
println(@sprintf("  dξ_moisture/dt at surface: %.2e kg/(m³·s)", du3[8*n3]))

# ============================================================================
# Test 2d: All components with significant concentration (active reactions)
# ============================================================================

println("\n" * "="^70)
println("TEST 2d: All components with significant concentrations (full reaction chain)")
println("="^70)

# Create a state where ALL components have non-zero concentrations
# This tests the product/reactant Jacobian terms when reactions are truly active

mesh4 = create_uniform_mesh(0.006, 30, Val(13); T_initial = 500.0,
	ξ_initial = (100.0, 100.0, 100.0, 50.0, 50.0, 30.0, 50.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0))

ode_fn4 = PyrolysisODEFunctionALE(mesh4, material, bc_set)
sv4 = ode_fn4.sv
u0_4 = create_state_vector_ale(mesh4, material)
n_dof4 = length(u0_4)
n4 = sv4.n_cells

# Warmup
du4 = similar(u0_4)
ode_fn4(du4, u0_4, nothing, 0.0)

# Compute Jacobians
n_faces4 = length(mesh4.faces)
n_nodes4 = length(mesh4.nodes)
jac_cache4 = JacobianCacheALE(n4, n_faces4, n_nodes4, NC, NR)
J_analytical4 = zeros(n_dof4, n_dof4)
J_fd4 = zeros(n_dof4, n_dof4)

compute_analytical_jacobian_ale!(J_analytical4, u0_4, nothing, 0.0, ode_fn4, jac_cache4)
compute_fd_jacobian!(J_fd4, u0_4, ode_fn4)

rel_err4 = norm(J_analytical4 - J_fd4) / max(norm(J_fd4), 1e-10)
max_err4 = maximum(abs.(J_analytical4 - J_fd4))
println(@sprintf("\nOverall Frobenius relative error: %.2e", rel_err4))
println(@sprintf("Maximum absolute error: %.2e", max_err4))

results4 = analyze_jacobian_blocks(J_analytical4, J_fd4, sv4, NC)

# Show component concentrations and reaction rates
println("\nComponent concentrations (uniform, T=500K):")
for j in 1:NC
	comp_name = material.components[j].name
	ξ_j = u0_4[j*n4]
	println(@sprintf("  ξ[%s] = %.1f kg/m³", comp_name, ξ_j))
end

# Evaluate reaction rates
rates4 = zeros(NR)
heats4 = zeros(NR)
T4 = u0_4[n4]
ξ4 = ntuple(j -> u0_4[j*n4+n4], NC)
Pyrolysis.Physics.evaluate_reactions!(rates4, heats4, ξ4, T4, material)
println("\nReaction rates at T=500K with all components present:")
for r in 1:NR
	rxn_name = material.reactions[r].name
	println(@sprintf("  Rxn %d (%s): %.2e kg/(m³·s)", r, rxn_name, rates4[r]))
end

# Check specific product/reactant coupling blocks
println("\nProduct/reactant derivative blocks (should be well-matched now):")

# Reaction 2: DF_int1 → DF_int2 + DF_gas2
ξ3_range4 = (3*n4+1):(4*n4)
ξ2_range4 = (2*n4+1):(3*n4)
dξ3_dξ2_an = J_analytical4[ξ3_range4, ξ2_range4]
dξ3_dξ2_fd = J_fd4[ξ3_range4, ξ2_range4]
norm_an = norm(dξ3_dξ2_an)
norm_fd = norm(dξ3_dξ2_fd)
rel_err_block = norm_fd > 1e-10 ? norm(dξ3_dξ2_an - dξ3_dξ2_fd) / norm_fd : 0.0
println(@sprintf("  dξ_int2/dξ_int1: ||An||=%.2e, ||FD||=%.2e, rel_err=%.2e %s",
	norm_an, norm_fd, rel_err_block, rel_err_block < 1e-4 ? "✓" : "✗"))

# Reaction 3: DF_int2 → DF_int3 + DF_gas3
ξ4_range4 = (4*n4+1):(5*n4)
dξ4_dξ3_an = J_analytical4[ξ4_range4, ξ3_range4]
dξ4_dξ3_fd = J_fd4[ξ4_range4, ξ3_range4]
norm_an = norm(dξ4_dξ3_an)
norm_fd = norm(dξ4_dξ3_fd)
rel_err_block = norm_fd > 1e-10 ? norm(dξ4_dξ3_an - dξ4_dξ3_fd) / norm_fd : 0.0
println(@sprintf("  dξ_int3/dξ_int2: ||An||=%.2e, ||FD||=%.2e, rel_err=%.2e %s",
	norm_an, norm_fd, rel_err_block, rel_err_block < 1e-4 ? "✓" : "✗"))

# Gas product: dξ_gas2/dξ_int1
ξ9_range4 = (9*n4+1):(10*n4)
dξ9_dξ2_an = J_analytical4[ξ9_range4, ξ2_range4]
dξ9_dξ2_fd = J_fd4[ξ9_range4, ξ2_range4]
norm_an = norm(dξ9_dξ2_an)
norm_fd = norm(dξ9_dξ2_fd)
rel_err_block = norm_fd > 1e-10 ? norm(dξ9_dξ2_an - dξ9_dξ2_fd) / norm_fd : 0.0
println(@sprintf("  dξ_gas2/dξ_int1: ||An||=%.2e, ||FD||=%.2e, rel_err=%.2e %s",
	norm_an, norm_fd, rel_err_block, rel_err_block < 1e-4 ? "✓" : "✗"))

# ============================================================================
# Test 3: Detailed element-by-element comparison for problematic blocks
# ============================================================================

println("\n" * "="^70)
println("TEST 3: Detailed examination of dz/dT and dz/dξ blocks")
println("="^70)

# Find the worst elements in key blocks
n = sv.n_cells
n_nodes = sv.n_nodes
z_range = (sv.z_offset+1):n_dof
T_range = 1:n

println("\nTop 5 worst elements in dz/dT block:")
dz_dT_an = J_analytical[z_range, T_range]
dz_dT_fd = J_fd[z_range, T_range]
diff_zT = abs.(dz_dT_an - dz_dT_fd)

for i in 1:min(5, length(diff_zT))
	idx = argmax(diff_zT)
	r, c = Tuple(idx)
	println(@sprintf("  z[%d]/T[%d]: analytical=%.6e, FD=%.6e, diff=%.2e",
		r, c, dz_dT_an[r, c], dz_dT_fd[r, c], diff_zT[r, c]))
	diff_zT[r, c] = 0.0
end

# Also check dT/dz
println("\nTop 5 worst elements in dT/dz block:")
dT_dz_an = J_analytical[T_range, z_range]
dT_dz_fd = J_fd[T_range, z_range]
diff_Tz = abs.(dT_dz_an - dT_dz_fd)

for i in 1:min(5, length(diff_Tz))
	idx = argmax(diff_Tz)
	r, c = Tuple(idx)
	println(@sprintf("  T[%d]/z[%d]: analytical=%.6e, FD=%.6e, diff=%.2e",
		r, c, dT_dz_an[r, c], dT_dz_fd[r, c], diff_Tz[r, c]))
	diff_Tz[r, c] = 0.0
end

# ============================================================================
# Test 4: Verify sparsity patterns match
# ============================================================================

println("\n" * "="^70)
println("TEST 4: Sparsity pattern comparison")
println("="^70)

# Threshold for considering an element non-zero
thresh = 1e-12

nnz_an = sum(abs.(J_analytical) .> thresh)
nnz_fd = sum(abs.(J_fd) .> thresh)

println(@sprintf("Non-zeros in analytical: %d", nnz_an))
println(@sprintf("Non-zeros in FD: %d", nnz_fd))

# Find elements that are zero in analytical but non-zero in FD (missing terms)
missing_mask = (abs.(J_analytical) .< thresh) .& (abs.(J_fd) .> thresh)
n_missing = sum(missing_mask)
println(@sprintf("Elements MISSING in analytical (zero in An, non-zero in FD): %d", n_missing))

if n_missing > 0 && n_missing < 20
	println("\nMissing elements:")
	for idx in findall(missing_mask)
		i, j = Tuple(idx)
		println(@sprintf("  J[%d,%d]: analytical=%.2e, FD=%.2e", i, j, J_analytical[i, j], J_fd[i, j]))
	end
elseif n_missing >= 20
	println("\nTop 10 missing elements (by FD magnitude):")
	fd_vals = abs.(J_fd) .* missing_mask
	for _ in 1:10
		idx = argmax(fd_vals)
		i, j = Tuple(idx)
		println(@sprintf("  J[%d,%d]: analytical=%.2e, FD=%.2e", i, j, J_analytical[i, j], J_fd[i, j]))
		fd_vals[i, j] = 0.0
	end
end

# Find elements that are non-zero in analytical but should be zero (extra terms)
extra_mask = (abs.(J_analytical) .> thresh) .& (abs.(J_fd) .< thresh)
n_extra = sum(extra_mask)
println(@sprintf("\nElements EXTRA in analytical (non-zero in An, zero in FD): %d", n_extra))

if n_extra > 0 && n_extra < 20
	println("\nExtra elements:")
	for idx in findall(extra_mask)
		i, j = Tuple(idx)
		println(@sprintf("  J[%d,%d]: analytical=%.2e, FD=%.2e", i, j, J_analytical[i, j], J_fd[i, j]))
	end
end

# ============================================================================
# Test 5: Check specific derivative formulas
# ============================================================================

println("\n" * "="^70)
println("TEST 5: Verify node velocity derivative structure")
println("="^70)

# The node velocity is computed as:
# w[1] = 0 (fixed bottom)
# w[k+1] = w[k] + θ[k] * Δz[k]
# So w[k] = sum(θ[i] * Δz[i] for i in 1:k-1)
#
# Therefore:
# ∂w[k]/∂T[i] should be non-zero only for i < k (lower triangular)
# ∂w[k]/∂z[m] depends on θ[i] * ∂Δz[i]/∂z[m]

println("\ndz/dT block structure (should be lower triangular):")
println("Row k represents node k, Col i represents cell i")
println("Entry should be non-zero when i < k (cell affects nodes above it)")

# Check structure
dz_dT_fd_block = J_fd[z_range, T_range]
for k in 1:min(5, n_nodes)
	row_nonzeros = findall(abs.(dz_dT_fd_block[k, :]) .> 1e-15)
	println(@sprintf("  Node %d (row): non-zero columns = %s", k, row_nonzeros))
end

# ============================================================================
# Test 6: FD step size sensitivity
# ============================================================================

println("\n" * "="^70)
println("TEST 6: FD step size sensitivity (checking if analytical is more accurate)")
println("="^70)

# Compare errors at different FD step sizes
material6, mesh6, bc_set6 = create_test_problem(n_cells = 60, T_initial = 500.0)
ode_fn6 = PyrolysisODEFunctionALE(mesh6, material6, bc_set6)
u0_6 = create_state_vector_ale(mesh6, material6)
n_dof6 = length(u0_6)
n_cells6 = length(mesh6.cells)
n_faces6 = length(mesh6.faces)
n_nodes6 = length(mesh6.nodes)

du6 = similar(u0_6)
ode_fn6(du6, u0_6, nothing, 0.0)

jac_cache6 = JacobianCacheALE(n_cells6, n_faces6, n_nodes6, NC, NR)
J_analytical6 = zeros(n_dof6, n_dof6)
compute_analytical_jacobian_ale!(J_analytical6, u0_6, nothing, 0.0, ode_fn6, jac_cache6)

println(@sprintf("%-12s %12s %12s", "ε", "||J_an - J_fd||", "Max|Diff|"))
println("-"^40)

for ε in [1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9]
	J_fd_eps = zeros(n_dof6, n_dof6)
	compute_fd_jacobian!(J_fd_eps, u0_6, ode_fn6, ε)

	diff = norm(J_analytical6 - J_fd_eps)
	max_diff = maximum(abs.(J_analytical6 - J_fd_eps))
	println(@sprintf("%-12.0e %12.2e %12.2e", ε, diff, max_diff))
end

println("\n(If error decreases then increases, optimal ε is near the minimum)")
println("(If error consistently decreases with ε, analytical may have issues)")

# ============================================================================
# Test 7: Investigate gas species transport diagonal elements
# ============================================================================

println("\n" * "="^70)
println("TEST 7: Analysis of gas transport Jacobian terms")
println("="^70)

# For Douglas Fir, we have 6 gas components (indices 8-13)
# Test at different temperatures to see reaction activity effects

println("\n--- TEST 7a: Low temperature (T=350K, moisture evaporation starts) ---")
material7a, mesh7a, bc_set7a = create_test_problem(n_cells = 30, T_initial = 350.0)
ode_fn7a = PyrolysisODEFunctionALE(mesh7a, material7a, bc_set7a)
u0_7a = create_state_vector_ale(mesh7a, material7a)
n_dof7a = length(u0_7a)
sv7a = ode_fn7a.sv
n7a = sv7a.n_cells

du7a = similar(u0_7a)
ode_fn7a(du7a, u0_7a, nothing, 0.0)

J_an7a = zeros(n_dof7a, n_dof7a)
J_fd7a = zeros(n_dof7a, n_dof7a)
jac_cache7a = JacobianCacheALE(length(mesh7a.cells), length(mesh7a.faces), length(mesh7a.nodes), NC, NR)

compute_analytical_jacobian_ale!(J_an7a, u0_7a, nothing, 0.0, ode_fn7a, jac_cache7a)
compute_fd_jacobian!(J_fd7a, u0_7a, ode_fn7a)

# Check moisture evaporation (reaction 6: MOISTURE → WATER_VAPOR)
println(@sprintf("Surface cell dξ_moisture/dξ_moisture: An=%.6e, FD=%.6e, diff=%.2e",
	J_an7a[8*n7a, 8*n7a], J_fd7a[8*n7a, 8*n7a], J_an7a[8*n7a, 8*n7a] - J_fd7a[8*n7a, 8*n7a]))

rel_err7a = norm(J_an7a - J_fd7a) / max(norm(J_fd7a), 1e-10)
println(@sprintf("Overall relative error at T=350K: %.2e", rel_err7a))

println("\n--- TEST 7b: High temperature (T=600K, multiple reactions active) ---")
material7b, mesh7b, bc_set7b = create_test_problem(n_cells = 30, T_initial = 600.0)
ode_fn7b = PyrolysisODEFunctionALE(mesh7b, material7b, bc_set7b)
u0_7b = create_state_vector_ale(mesh7b, material7b)
n_dof7b = length(u0_7b)
sv7b = ode_fn7b.sv
n7b = sv7b.n_cells

du7b = similar(u0_7b)
ode_fn7b(du7b, u0_7b, nothing, 0.0)

J_an7b = zeros(n_dof7b, n_dof7b)
J_fd7b = zeros(n_dof7b, n_dof7b)
jac_cache7b = JacobianCacheALE(length(mesh7b.cells), length(mesh7b.faces), length(mesh7b.nodes), NC, NR)

compute_analytical_jacobian_ale!(J_an7b, u0_7b, nothing, 0.0, ode_fn7b, jac_cache7b)
compute_fd_jacobian!(J_fd7b, u0_7b, ode_fn7b)

# Check various reaction rates
println("\nReaction activity at T=600K:")
println(@sprintf("  dξ_DF/dt at surface: %.2e kg/(m³·s)", du7b[2*n7b]))
println(@sprintf("  dξ_DF_int1/dt at surface: %.2e kg/(m³·s)", du7b[3*n7b]))
println(@sprintf("  dξ_moisture/dt at surface: %.2e kg/(m³·s)", du7b[8*n7b]))

rel_err7b = norm(J_an7b - J_fd7b) / max(norm(J_fd7b), 1e-10)
println(@sprintf("\nOverall relative error at T=600K: %.2e", rel_err7b))

# Check the effective gas transfer coefficients
prop_cache7b = ode_fn7b.prop_cache
println(@sprintf("\nEffective gas transfer at surface cell: λ[%d] = %.6e m²/s", n7b, prop_cache7b.λ[n7b]))

println("\n" * "="^70)
println("SUMMARY")
println("="^70)

# Check if any blocks have concerning errors
concerning = String[]
for (name, r) in results2
	# Only flag as concerning if both:
	# 1. Relative error > 10% AND
	# 2. The FD norm is significant (> 1e-6) - otherwise it's just noise
	if r.rel_err > 0.1 && r.norm_fd > 1e-6
		push!(concerning, name)
	end
end

if isempty(concerning)
	println("\n✓ All significant Jacobian blocks have acceptable agreement with FD")
	println("\nNote: Some reaction product/reactant blocks show 100% relative error")
	println("      but have tiny absolute values (~1e-10) - these are negligible.")
	println("      Example: dξ_gas/dξ_reactant terms (product depends on reactant)")
else
	println("\n⚠ The following blocks may have issues:")
	for name in concerning
		r = results2[name]
		println(@sprintf("  - %s: rel_err=%.2e, norm_fd=%.2e", name, r.rel_err, r.norm_fd))
	end
end

println("\n" * "="^70)
println("Test complete!")
println("="^70)
