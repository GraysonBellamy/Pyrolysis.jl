#!/usr/bin/env julia
"""
Comprehensive Jacobian Verification Test Suite

This script performs exhaustive comparison of the analytical Jacobian against
finite differences across a wide range of conditions. It tests:

1. Multiple temperature ranges (300K to 800K)
2. Various decomposition states (virgin → partial → fully charred)
3. Presence of gases in the solid matrix
4. Darcy flow (pressure-driven transport)
5. Different radiation models (NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT)
6. ALE mesh deformation
7. Temperature gradients and concentration gradients

Usage:
	julia --project=. test/test_jacobian_comprehensive.jl

The script reports per-block and per-element errors, identifying any terms
where the analytical Jacobian disagrees with finite differences.
"""

using Pyrolysis
using Pyrolysis.Core
using Pyrolysis.Geometry
using Pyrolysis.Properties
using Pyrolysis.Physics
using Pyrolysis.TimeIntegration
using Pyrolysis.TimeIntegration:
	JacobianCacheALE, JacobianCache, StateVectorALE, StateVector,
	compute_analytical_jacobian_ale!, compute_analytical_jacobian!,
	PyrolysisODEFunctionALE, PyrolysisODEFunction,
	create_state_vector_ale, create_state_vector, unpack_state_ale!
using LinearAlgebra
using Printf
using SparseArrays
using Statistics: mean, std

# ============================================================================
# Configuration
# ============================================================================

const VERBOSE = true                    # Print detailed output
const TOLERANCE_REL = 1e-3              # Relative error threshold for "pass"
const TOLERANCE_ABS = 1e-10             # Absolute error threshold for small values
const FD_STEP_SIZES = [1e-4, 1e-6, 1e-7, 1e-8]  # FD step sizes to test convergence
const DEFAULT_FD_STEP = 1e-7            # Default FD step for main comparisons

# ============================================================================
# Test Results Tracking
# ============================================================================

mutable struct TestResult
	name::String
	passed::Bool
	rel_error::Float64
	abs_error::Float64
	norm_fd::Float64
	worst_element::Tuple{Int, Int}
	worst_an::Float64
	worst_fd::Float64
	details::String
end

const ALL_RESULTS = TestResult[]

function report_test(name, passed, rel_err, abs_err, norm_fd, worst_loc, worst_an, worst_fd; details = "")
	result = TestResult(name, passed, rel_err, abs_err, norm_fd, worst_loc, worst_an, worst_fd, details)
	push!(ALL_RESULTS, result)

	status = passed ? "✓ PASS" : "✗ FAIL"
	if VERBOSE
		println(@sprintf("  %-50s %s (rel=%.2e, abs=%.2e)", name, status, rel_err, abs_err))
	end
	return passed
end

# ============================================================================
# Material Definitions (Douglas Fir base from 07_douglas_fir_cone.jl)
# ============================================================================

# Temperature-dependent heat capacity functions
const df_cp = CallableProperty(T -> -2131.0 + 9.119 * T)
const df_int1_cp = CallableProperty(T -> -1983.0 + 8.7 * T)
const df_int2_cp = CallableProperty(T -> -1419.0 + 7.0 * T)
const df_int3_cp = CallableProperty(T -> 1.0 + 2.8 * T)
const df_int4_cp = CallableProperty(T -> 176.0 + 2.3 * T)
const char_cp = CallableProperty(T -> 347.2 + 1.784 * T)
const moisture_cp = CallableProperty(T -> 5.2e3 - 6.7 * T + 1.1e-2 * T^2)
const water_vapor_cp = CallableProperty(T -> 2.4e3 - 1.6 * T + 2.0e-3 * T^2)

# Gas transfer and permeability coefficients
const λ_char = 5e-5
const λ_wood = 1e-5
const λ_gas = 1e-4
const κ_virgin = 1e-15
const κ_char = 1e-12

# Heat of sorption (state-dependent)
const heat_of_sorption = StateDependentProperty((T, ξ) -> 2673.82e3)

"""
Create Douglas Fir material without Darcy flow (diffusion-only transport).
"""
function create_douglas_fir_diffusion_only()
	return Material(
		name = :DouglasFir_DiffusionOnly,
		components = (
			SolidComponent(:DF, ρ = 530.0, c = df_cp, k = 0.12, ε = 0.9, α = 0.0, λ = λ_wood, κ = 0.0),
			SolidComponent(:DF_int1, ρ = 507.2, c = df_int1_cp, k = 0.14, ε = 0.9, α = 0.0, λ = λ_wood * 1.2, κ = 0.0),
			SolidComponent(:DF_int2, ρ = 419.8, c = df_int2_cp, k = 0.15, ε = 0.9, α = 0.0, λ = λ_wood * 1.5, κ = 0.0),
			SolidComponent(:DF_int3, ρ = 203.2, c = df_int3_cp, k = 0.16, ε = 0.9, α = 0.0, λ = λ_wood * 2.0, κ = 0.0),
			SolidComponent(:DF_int4, ρ = 176.6, c = df_int4_cp, k = 0.18, ε = 0.9, α = 0.0, λ = λ_char * 0.8, κ = 0.0),
			SolidComponent(:Char, ρ = 150.0, c = char_cp, k = 0.20, ε = 0.95, α = 0.0, λ = λ_char, κ = 0.0),
			LiquidComponent(:MOISTURE, ρ = 1000.0, c = moisture_cp, k = 0.6, ε = 0.9, α = 0.0),
			GaseousComponent(:DF_gas1, M = 0.060, c = 2000.0, k = 0.03, λ = λ_gas),
			GaseousComponent(:DF_gas2, M = 0.044, c = 1800.0, k = 0.03, λ = λ_gas),
			GaseousComponent(:DF_gas3, M = 0.070, c = 2200.0, k = 0.025, λ = λ_gas * 0.8),
			GaseousComponent(:DF_gas4, M = 0.030, c = 1600.0, k = 0.035, λ = λ_gas),
			GaseousComponent(:DF_gas5, M = 0.028, c = 1500.0, k = 0.035, λ = λ_gas),
			GaseousComponent(:WATER_VAPOR, M = 0.018, c = water_vapor_cp, k = 0.04, λ = λ_gas * 1.2),
		),
		reactions = (
			Reaction(:pyrolysis_step1, 1 => (2, 8), A = 74816.95, E = 71.59e3, h = -20067.0, n = 1.0, yields = (0.95, 0.05)),
			Reaction(:pyrolysis_step2, 2 => (3, 9), A = 23550492839.0, E = 140.9e3, h = -8264.0, n = 1.0, yields = (0.80, 0.20)),
			Reaction(:pyrolysis_step3, 3 => (4, 10), A = 7.744618e12, E = 181.0e3, h = 273999.0, n = 1.0, yields = (0.37, 0.63)),
			Reaction(:pyrolysis_step4, 4 => (5, 11), A = 103276.14, E = 98.94e3, h = -224430.0, n = 1.0, yields = (0.79, 0.21)),
			Reaction(:pyrolysis_step5, 5 => (6, 12), A = 6.397348, E = 56.01e3, h = 555.6, n = 1.0, yields = (0.74, 0.26)),
			Reaction(:moisture_evaporation, 7 => 13, A = 2.997e5, E = 5.0e4, h = heat_of_sorption, n = 1.0),
		),
		conductivity_mixing = WEIGHTED,
		β = 0.5,
		gas_transfer_mixing = WEIGHTED,
		β_λ = 0.5,
	)
end

"""
Create Douglas Fir material with Darcy flow (pressure-driven transport).
"""
function create_douglas_fir_with_darcy()
	return Material(
		name = :DouglasFir_Darcy,
		components = (
			SolidComponent(:DF, ρ = 530.0, c = df_cp, k = 0.12, ε = 0.9, α = 0.0, λ = λ_wood, κ = κ_virgin),
			SolidComponent(:DF_int1, ρ = 507.2, c = df_int1_cp, k = 0.14, ε = 0.9, α = 0.0, λ = λ_wood * 1.2, κ = κ_virgin * 2),
			SolidComponent(:DF_int2, ρ = 419.8, c = df_int2_cp, k = 0.15, ε = 0.9, α = 0.0, λ = λ_wood * 1.5, κ = κ_virgin * 5),
			SolidComponent(:DF_int3, ρ = 203.2, c = df_int3_cp, k = 0.16, ε = 0.9, α = 0.0, λ = λ_wood * 2.0, κ = κ_virgin * 10),
			SolidComponent(:DF_int4, ρ = 176.6, c = df_int4_cp, k = 0.18, ε = 0.9, α = 0.0, λ = λ_char * 0.8, κ = κ_char * 0.5),
			SolidComponent(:Char, ρ = 150.0, c = char_cp, k = 0.20, ε = 0.95, α = 0.0, λ = λ_char, κ = κ_char),
			LiquidComponent(:MOISTURE, ρ = 1000.0, c = moisture_cp, k = 0.6, ε = 0.9, α = 0.0),
			GaseousComponent(:DF_gas1, M = 0.060, c = 2000.0, k = 0.03, λ = λ_gas),
			GaseousComponent(:DF_gas2, M = 0.044, c = 1800.0, k = 0.03, λ = λ_gas),
			GaseousComponent(:DF_gas3, M = 0.070, c = 2200.0, k = 0.025, λ = λ_gas * 0.8),
			GaseousComponent(:DF_gas4, M = 0.030, c = 1600.0, k = 0.035, λ = λ_gas),
			GaseousComponent(:DF_gas5, M = 0.028, c = 1500.0, k = 0.035, λ = λ_gas),
			GaseousComponent(:WATER_VAPOR, M = 0.018, c = water_vapor_cp, k = 0.04, λ = λ_gas * 1.2),
		),
		reactions = (
			Reaction(:pyrolysis_step1, 1 => (2, 8), A = 74816.95, E = 71.59e3, h = -20067.0, n = 1.0, yields = (0.95, 0.05)),
			Reaction(:pyrolysis_step2, 2 => (3, 9), A = 23550492839.0, E = 140.9e3, h = -8264.0, n = 1.0, yields = (0.80, 0.20)),
			Reaction(:pyrolysis_step3, 3 => (4, 10), A = 7.744618e12, E = 181.0e3, h = 273999.0, n = 1.0, yields = (0.37, 0.63)),
			Reaction(:pyrolysis_step4, 4 => (5, 11), A = 103276.14, E = 98.94e3, h = -224430.0, n = 1.0, yields = (0.79, 0.21)),
			Reaction(:pyrolysis_step5, 5 => (6, 12), A = 6.397348, E = 56.01e3, h = 555.6, n = 1.0, yields = (0.74, 0.26)),
			Reaction(:moisture_evaporation, 7 => 13, A = 2.997e5, E = 5.0e4, h = heat_of_sorption, n = 1.0),
		),
		conductivity_mixing = WEIGHTED,
		β = 0.5,
		gas_transfer_mixing = WEIGHTED,
		β_λ = 0.5,
		permeability_mixing = WEIGHTED,
		β_κ = 0.5,
	)
end

"""
Create Douglas Fir material with Beer-Lambert radiation absorption.
"""
function create_douglas_fir_beer_lambert()
	return Material(
		name = :DouglasFir_BeerLambert,
		components = (
			SolidComponent(:DF, ρ = 530.0, c = df_cp, k = 0.12, ε = 0.9, α = 5000.0, λ = λ_wood),
			SolidComponent(:DF_int1, ρ = 507.2, c = df_int1_cp, k = 0.14, ε = 0.9, α = 4000.0, λ = λ_wood * 1.2),
			SolidComponent(:DF_int2, ρ = 419.8, c = df_int2_cp, k = 0.15, ε = 0.9, α = 3000.0, λ = λ_wood * 1.5),
			SolidComponent(:DF_int3, ρ = 203.2, c = df_int3_cp, k = 0.16, ε = 0.9, α = 2000.0, λ = λ_wood * 2.0),
			SolidComponent(:DF_int4, ρ = 176.6, c = df_int4_cp, k = 0.18, ε = 0.9, α = 1500.0, λ = λ_char * 0.8),
			SolidComponent(:Char, ρ = 150.0, c = char_cp, k = 0.20, ε = 0.95, α = 10000.0, λ = λ_char),
			LiquidComponent(:MOISTURE, ρ = 1000.0, c = moisture_cp, k = 0.6, ε = 0.9, α = 100.0),
			GaseousComponent(:DF_gas1, M = 0.060, c = 2000.0, k = 0.03, λ = λ_gas),
			GaseousComponent(:DF_gas2, M = 0.044, c = 1800.0, k = 0.03, λ = λ_gas),
			GaseousComponent(:DF_gas3, M = 0.070, c = 2200.0, k = 0.025, λ = λ_gas * 0.8),
			GaseousComponent(:DF_gas4, M = 0.030, c = 1600.0, k = 0.035, λ = λ_gas),
			GaseousComponent(:DF_gas5, M = 0.028, c = 1500.0, k = 0.035, λ = λ_gas),
			GaseousComponent(:WATER_VAPOR, M = 0.018, c = water_vapor_cp, k = 0.04, λ = λ_gas * 1.2),
		),
		reactions = (
			Reaction(:pyrolysis_step1, 1 => (2, 8), A = 74816.95, E = 71.59e3, h = -20067.0, n = 1.0, yields = (0.95, 0.05)),
			Reaction(:pyrolysis_step2, 2 => (3, 9), A = 23550492839.0, E = 140.9e3, h = -8264.0, n = 1.0, yields = (0.80, 0.20)),
			Reaction(:pyrolysis_step3, 3 => (4, 10), A = 7.744618e12, E = 181.0e3, h = 273999.0, n = 1.0, yields = (0.37, 0.63)),
			Reaction(:pyrolysis_step4, 4 => (5, 11), A = 103276.14, E = 98.94e3, h = -224430.0, n = 1.0, yields = (0.79, 0.21)),
			Reaction(:pyrolysis_step5, 5 => (6, 12), A = 6.397348, E = 56.01e3, h = 555.6, n = 1.0, yields = (0.74, 0.26)),
			Reaction(:moisture_evaporation, 7 => 13, A = 2.997e5, E = 5.0e4, h = heat_of_sorption, n = 1.0),
		),
		conductivity_mixing = WEIGHTED,
		β = 0.5,
		gas_transfer_mixing = WEIGHTED,
		β_λ = 0.5,
	)
end

"""
Create a simpler 3-component material for focused testing.
"""
function create_simple_material()
	return Material(
		name = :Simple3,
		components = (
			SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2, ε = 0.9, α = 0.0, λ = 1e-5, κ = 0.0),
			SolidComponent(:char, ρ = 200.0, c = 1000.0, k = 0.1, ε = 0.95, α = 0.0, λ = 5e-5, κ = 0.0),
			GaseousComponent(:gas, M = 0.044, c = 1500.0, k = 0.03, λ = 1e-4),
		),
		reactions = (
			Reaction(:decomposition, 1 => (2, 3), A = 1e12, E = 150e3, h = 500e3, n = 1.0, yields = (0.2, 0.8)),
		),
		conductivity_mixing = WEIGHTED,
		β = 0.5,
		gas_transfer_mixing = WEIGHTED,
		β_λ = 0.5,
	)
end

# ============================================================================
# Finite Difference Jacobian Computation
# ============================================================================

"""
Compute finite difference Jacobian with specified step size.
"""
function compute_fd_jacobian!(J_fd, u, ode_fn, ε = DEFAULT_FD_STEP; debug_col = 0)
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

		# Debug specific column
		if j == debug_col && debug_col > 0
			println("FD DEBUG col $j: u[$j]=$(u[j]), u_pert[$j]=$(u_pert[j])")
			println("FD DEBUG: du_base[40]=$(du_base[40]), du_pert[40]=$(du_pert[40])")
			println("FD DEBUG: J_fd[40,$j] = $(J_fd[40, j])")
		end
	end
end

"""
Compute central difference Jacobian for higher accuracy.
"""
function compute_fd_jacobian_central!(J_fd, u, ode_fn, ε = DEFAULT_FD_STEP)
	n_dof = length(u)
	u_plus = copy(u)
	u_minus = copy(u)
	du_plus = similar(u)
	du_minus = similar(u)

	for j in 1:n_dof
		u_plus .= u
		u_minus .= u
		u_plus[j] += ε
		u_minus[j] -= ε
		ode_fn(du_plus, u_plus, nothing, 0.0)
		ode_fn(du_minus, u_minus, nothing, 0.0)
		@views J_fd[:, j] .= (du_plus .- du_minus) ./ (2ε)
	end
end

# ============================================================================
# Block Analysis Functions
# ============================================================================

"""
Define block ranges for the Jacobian given state vector structure.
"""
function get_block_ranges(sv, NC, n_dof; is_ale = true)
	n = sv.n_cells

	blocks = Dict{String, Tuple{UnitRange{Int}, UnitRange{Int}}}()

	# Temperature block
	T_range = 1:n
	blocks["dT/dT"] = (T_range, T_range)

	# Species blocks
	ξ_ranges = [n + (j-1)*n .+ (1:n) for j in 1:NC]

	for j in 1:NC
		blocks["dT/dξ$j"] = (T_range, ξ_ranges[j])
		blocks["dξ$j/dT"] = (ξ_ranges[j], T_range)
		for k in 1:NC
			blocks["dξ$j/dξ$k"] = (ξ_ranges[j], ξ_ranges[k])
		end
	end

	# ALE blocks (node positions)
	if is_ale && hasfield(typeof(sv), :z_offset)
		z_range = (sv.z_offset+1):n_dof
		blocks["dT/dz"] = (T_range, z_range)
		blocks["dz/dT"] = (z_range, T_range)
		blocks["dz/dz"] = (z_range, z_range)

		for j in 1:NC
			blocks["dξ$j/dz"] = (ξ_ranges[j], z_range)
			blocks["dz/dξ$j"] = (z_range, ξ_ranges[j])
		end
	end

	return blocks
end

"""
Analyze Jacobian block-by-block and report errors.
"""
function analyze_jacobian_blocks(J_an, J_fd, blocks; verbose = true, test_name = "")
	results = Dict{String, NamedTuple}()
	all_passed = true

	for (name, (rows, cols)) in blocks
		J_an_block = J_an[rows, cols]
		J_fd_block = J_fd[rows, cols]

		diff = J_an_block - J_fd_block
		abs_diff = abs.(diff)

		# Error metrics
		max_abs_err = maximum(abs_diff)
		mean_abs_err = mean(abs_diff)
		frobenius_err = norm(diff)
		norm_fd = norm(J_fd_block)
		norm_an = norm(J_an_block)
		rel_err = norm_fd > TOLERANCE_ABS ? frobenius_err / norm_fd : 0.0

		# Find worst element
		idx = argmax(abs_diff)
		worst_i, worst_j = Tuple(idx)

		# Determine pass/fail
		# Pass if any of:
		#   1. Relative error is small
		#   2. Both FD norm and max absolute error are tiny (both matrices near zero)
		#   3. Max absolute error is tiny (covers case where analytical is zero but FD shows FD noise)
		#      Use a slightly higher threshold (1e-6) for this case since FD has numerical noise
		passed = (rel_err < TOLERANCE_REL) ||
				 (norm_fd < TOLERANCE_ABS && max_abs_err < TOLERANCE_ABS) ||
				 (max_abs_err < 1e-6)  # Absolute tolerance for near-zero elements

		if !passed
			all_passed = false
		end

		results[name] = (
			passed = passed,
			max_abs_err = max_abs_err,
			mean_abs_err = mean_abs_err,
			frobenius_err = frobenius_err,
			rel_err = rel_err,
			norm_fd = norm_fd,
			norm_an = norm_an,
			worst_loc = (worst_i, worst_j),
			worst_an = J_an_block[worst_i, worst_j],
			worst_fd = J_fd_block[worst_i, worst_j],
		)

		# Record result
		full_name = isempty(test_name) ? name : "$test_name :: $name"
		report_test(full_name, passed, rel_err, max_abs_err, norm_fd,
			(rows[worst_i], cols[worst_j]),
			J_an_block[worst_i, worst_j],
			J_fd_block[worst_i, worst_j])
	end

	return results, all_passed
end

"""
Find elements with largest discrepancy.
"""
function find_worst_elements(J_an, J_fd; n_worst = 10)
	diff = abs.(J_an - J_fd)
	worst = []

	for _ in 1:n_worst
		idx = argmax(diff)
		i, j = Tuple(idx)
		if diff[i, j] > 0
			push!(worst, (i = i, j = j, an = J_an[i, j], fd = J_fd[i, j], diff = diff[i, j]))
			diff[i, j] = 0.0
		end
	end

	return worst
end

# ============================================================================
# State Modification Functions
# ============================================================================
# Note: State is stored in the state vector u0, not directly in mesh.
# The pattern is:
#   u0[i] = T_i for i = 1:n_cells
#   u0[n_cells + (j-1)*n_cells + i] = ξ_j,i for component j at cell i
# For ALE, node positions are at the end.

"""
Modify state vector to have uniform temperature and virgin wood composition.
Returns a function that modifies the state vector u0 after it's created.
"""
function modify_state_virgin(u0, sv, NC, T_uniform, material)
	n = sv.n_cells

	# Set temperatures
	for i in 1:n
		u0[i] = T_uniform
	end

	# Set virgin composition (only first component non-zero)
	ρ_virgin = material.components[1].density(T_uniform)
	for j in 1:NC
		for i in 1:n
			u0[n+(j-1)*n+i] = (j == 1) ? ρ_virgin : 0.0
		end
	end

	return u0
end

"""
Modify state vector to have partial decomposition with temperature gradient.
"""
function modify_state_partial_decomposition(u0, sv, NC, T_bottom, T_top, material)
	n = sv.n_cells

	for i in 1:n
		# Temperature gradient: cool at bottom (cell 1), hot at top (cell n)
		progress = (i - 1) / (n - 1)
		T_i = T_bottom + (T_top - T_bottom) * progress
		u0[i] = T_i

		# Decomposition progress (more at top where it's hotter)
		decomp = max(0.0, (i - n/3) / (2*n/3))^1.5

		# Component concentrations
		ξ_DF = 530.0 * (1.0 - decomp)
		ξ_int1 = 530.0 * 0.95 * decomp * (1.0 - decomp)^0.5
		ξ_int2 = 530.0 * 0.80 * 0.95 * decomp^1.5 * (1.0 - decomp)^0.3
		ξ_int3 = 530.0 * 0.37 * decomp^2 * (1.0 - decomp)^0.2
		ξ_int4 = 530.0 * 0.30 * decomp^2.5 * 0.5
		ξ_char = 70.0 * decomp^3
		ξ_moisture = 79.5 * (1.0 - decomp)^2
		ξ_gas1 = 5.0 * decomp * (1.0 - decomp)
		ξ_gas2 = 10.0 * decomp * (1.0 - decomp)^0.5
		ξ_gas3 = 15.0 * decomp^1.5 * (1.0 - decomp)^0.5
		ξ_gas4 = 5.0 * decomp^2 * (1.0 - decomp)^0.3
		ξ_gas5 = 3.0 * decomp^2.5 * 0.5
		ξ_water = 10.0 * decomp * (1.0 - decomp)^1.5

		if NC == 13
			u0[n+i]    = ξ_DF
			u0[2*n+i]  = ξ_int1
			u0[3*n+i]  = ξ_int2
			u0[4*n+i]  = ξ_int3
			u0[5*n+i]  = ξ_int4
			u0[6*n+i]  = ξ_char
			u0[7*n+i]  = ξ_moisture
			u0[8*n+i]  = ξ_gas1
			u0[9*n+i]  = ξ_gas2
			u0[10*n+i] = ξ_gas3
			u0[11*n+i] = ξ_gas4
			u0[12*n+i] = ξ_gas5
			u0[13*n+i] = ξ_water
		elseif NC == 3
			u0[n+i]   = 800.0 * (1.0 - decomp)
			u0[2*n+i] = 100.0 * decomp
			u0[3*n+i] = 10.0 * decomp * (1.0 - decomp)
		end
	end

	return u0
end

"""
Modify state vector to have all components with non-zero concentrations.
"""
function modify_state_all_components(u0, sv, NC, T_uniform)
	n = sv.n_cells

	for i in 1:n
		u0[i] = T_uniform
	end

	if NC == 13
		concs = (100.0, 100.0, 100.0, 50.0, 50.0, 30.0, 50.0,
			10.0, 10.0, 10.0, 10.0, 10.0, 10.0)
		for j in 1:NC
			for i in 1:n
				u0[n+(j-1)*n+i] = concs[j]
			end
		end
	elseif NC == 3
		concs = (500.0, 100.0, 5.0)
		for j in 1:NC
			for i in 1:n
				u0[n+(j-1)*n+i] = concs[j]
			end
		end
	else
		for j in 1:NC
			for i in 1:n
				u0[n+(j-1)*n+i] = j <= NC÷2 ? 100.0 : 10.0
			end
		end
	end

	return u0
end

"""
Modify state vector to have high gas concentrations.
"""
function modify_state_high_gas(u0, sv, NC, T_uniform)
	n = sv.n_cells

	for i in 1:n
		u0[i] = T_uniform
	end

	if NC == 13
		# High gas concentrations
		concs = (50.0, 20.0, 20.0, 10.0, 10.0, 50.0, 10.0,
			50.0, 50.0, 50.0, 30.0, 30.0, 40.0)
		for j in 1:NC
			for i in 1:n
				u0[n+(j-1)*n+i] = concs[j]
			end
		end
	elseif NC == 3
		concs = (200.0, 100.0, 50.0)
		for j in 1:NC
			for i in 1:n
				u0[n+(j-1)*n+i] = concs[j]
			end
		end
	else
		for j in 1:NC
			for i in 1:n
				u0[n+(j-1)*n+i] = j > NC÷2 ? 50.0 : 100.0
			end
		end
	end

	return u0
end

"""
Modify mesh node positions for ALE deformation testing.
Also modifies state vector for partial decomposition.
"""
function modify_mesh_deformed(mesh, u0, sv, NC, T_bottom, T_top, shrinkage, material)
	n = sv.n_cells
	n_nodes = length(mesh.nodes)

	# Apply mesh deformation (shrinkage increases toward top)
	for (k, node) in enumerate(mesh.nodes)
		progress = (k - 1) / (n_nodes - 1)
		z_new = node.z * (1.0 - shrinkage * progress)
		mesh.nodes[k] = Node1D(node.id, z_new, node.z_ref, node.is_boundary, node.boundary_id)
	end

	# Update geometric quantities
	Pyrolysis.Geometry.update_geometry!(mesh)

	# Update node positions in state vector (for ALE)
	if hasfield(typeof(sv), :z_offset)
		for k in 1:n_nodes
			u0[sv.z_offset+k] = mesh.nodes[k].z
		end
	end

	# Set temperature gradient and partial decomposition
	modify_state_partial_decomposition(u0, sv, NC, T_bottom, T_top, material)

	return mesh, u0
end

# ============================================================================
# Core Test Functions
# ============================================================================

"""
Run comprehensive Jacobian test for a given configuration.
"""
function test_jacobian_configuration(;
	name::String,
	material,
	bc_set,
	state_modifier,  # Function (u0, sv, NC, material) -> u0
	use_ale::Bool = true,
	radiation_model::RadiationModel = SURFACE_ABSORPTION,
	n_cells::Int = 30,
	mesh_modifier = nothing,  # Optional: (mesh, u0, sv, NC, material) -> (mesh, u0)
)
	NC = n_components(material)
	NR = n_reactions(material)

	println("\n" * "-"^70)
	println("Testing: $name")
	println("-"^70)

	# Create fresh mesh
	mesh = create_uniform_mesh(0.006, n_cells, Val(NC);
		T_initial = 300.0,
		ξ_initial = ntuple(j -> j == 1 ? 530.0 : (j == 7 && NC == 13 ? 80.0 : 0.0), NC))

	if use_ale
		# ALE mode
		ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set; radiation_model = radiation_model)
		sv = ode_fn.sv
		u0 = create_state_vector_ale(mesh, material)
		n_dof = length(u0)
		n_faces = length(mesh.faces)
		n_nodes = length(mesh.nodes)

		# Apply mesh modification if provided
		if mesh_modifier !== nothing
			mesh, u0 = mesh_modifier(mesh, u0, sv, NC, material)
		end

		# Apply state modification
		u0 = state_modifier(u0, sv, NC, material)

		# Sync state to mesh
		unpack_state_ale!(mesh, u0, sv)

		# Update property cache
		Pyrolysis.Properties.update_properties!(ode_fn.prop_cache, mesh, material)

		# Warmup RHS
		du = similar(u0)
		ode_fn(du, u0, nothing, 0.0)

		# Compute Jacobians
		jac_cache = JacobianCacheALE(n_cells, n_faces, n_nodes, NC, NR)
		J_an = zeros(n_dof, n_dof)
		J_fd = zeros(n_dof, n_dof)

		compute_analytical_jacobian_ale!(J_an, u0, nothing, 0.0, ode_fn, jac_cache)
		compute_fd_jacobian!(J_fd, u0, ode_fn)

		blocks = get_block_ranges(sv, NC, n_dof; is_ale = true)
	else
		# Non-ALE mode
		ode_fn = PyrolysisODEFunction(mesh, material, bc_set; radiation_model = radiation_model)
		sv = ode_fn.sv
		u0 = create_state_vector(mesh, NC)
		n_dof = length(u0)
		n_faces = length(mesh.faces)

		# Apply state modification
		u0 = state_modifier(u0, sv, NC, material)

		# Sync state to mesh
		Pyrolysis.TimeIntegration.unpack_state!(mesh, u0, sv)

		# Update property cache
		Pyrolysis.Properties.update_properties!(ode_fn.prop_cache, mesh, material)

		# Warmup RHS
		du = similar(u0)
		ode_fn(du, u0, nothing, 0.0)

		# Compute Jacobians
		jac_cache = JacobianCache(n_cells, n_faces, NC, NR)
		J_an = zeros(n_dof, n_dof)
		J_fd = zeros(n_dof, n_dof)

		compute_analytical_jacobian!(J_an, u0, nothing, 0.0, ode_fn, jac_cache)
		compute_fd_jacobian!(J_fd, u0, ode_fn)

		blocks = get_block_ranges(sv, NC, n_dof; is_ale = false)
	end

	# Overall metrics
	rel_err = norm(J_an - J_fd) / max(norm(J_fd), 1e-10)
	max_err = maximum(abs.(J_an - J_fd))

	println(@sprintf("Overall: ||J_an - J_fd|| / ||J_fd|| = %.2e, max|diff| = %.2e", rel_err, max_err))

	# Block-by-block analysis
	results, all_passed = analyze_jacobian_blocks(J_an, J_fd, blocks; verbose = VERBOSE, test_name = name)

	# Report worst elements if significant error
	if !all_passed || rel_err > TOLERANCE_REL
		println("\nWorst 5 elements:")
		worst = find_worst_elements(J_an, J_fd; n_worst = 5)
		for w in worst
			println(@sprintf("  J[%d,%d]: An=%.6e, FD=%.6e, diff=%.2e",
				w.i, w.j, w.an, w.fd, w.diff))
		end
	end

	return all_passed, rel_err, max_err
end

"""
Test FD step size sensitivity.
"""
function test_fd_step_sensitivity(material, bc_set, state_modifier; n_cells = 30)
	NC = n_components(material)
	NR = n_reactions(material)

	mesh = create_uniform_mesh(0.006, n_cells, Val(NC);
		T_initial = 500.0,
		ξ_initial = ntuple(j -> j == 1 ? 530.0 : 0.0, NC))

	ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)
	sv = ode_fn.sv
	u0 = create_state_vector_ale(mesh, material)

	# Apply state modification
	u0 = state_modifier(u0, sv, NC, material)
	unpack_state_ale!(mesh, u0, sv)
	Pyrolysis.Properties.update_properties!(ode_fn.prop_cache, mesh, material)

	n_dof = length(u0)
	n_faces = length(mesh.faces)
	n_nodes = length(mesh.nodes)

	du = similar(u0)
	ode_fn(du, u0, nothing, 0.0)

	jac_cache = JacobianCacheALE(n_cells, n_faces, n_nodes, NC, NR)
	J_an = zeros(n_dof, n_dof)
	compute_analytical_jacobian_ale!(J_an, u0, nothing, 0.0, ode_fn, jac_cache)

	println("\nFD Step Size Sensitivity Analysis:")
	println(@sprintf("%-12s %15s %15s", "ε", "||J_an - J_fd||", "Max|Diff|"))
	println("-"^45)

	errors = Float64[]
	for ε in FD_STEP_SIZES
		J_fd = zeros(n_dof, n_dof)
		compute_fd_jacobian!(J_fd, u0, ode_fn, ε)
		diff = norm(J_an - J_fd)
		max_diff = maximum(abs.(J_an - J_fd))
		push!(errors, diff)
		println(@sprintf("%-12.0e %15.2e %15.2e", ε, diff, max_diff))
	end

	# Check if analytical matches FD better at smaller ε (expected if FD limited)
	# or if error stays constant (analytical may have issues)
	if length(errors) >= 3
		if errors[end] < errors[1] * 0.1
			println("✓ Error decreases with ε → Analytical Jacobian is likely correct")
		elseif errors[end] > errors[1] * 0.5
			println("⚠ Error does NOT decrease significantly → Check analytical Jacobian")
		end
	end

	return errors
end

# ============================================================================
# Main Test Suite
# ============================================================================

function run_all_tests()
	println("="^70)
	println("COMPREHENSIVE JACOBIAN VERIFICATION TEST SUITE")
	println("="^70)
	println("Comparing Analytical Jacobian vs Finite Differences")
	println("Douglas Fir material (13 components, 6 reactions)")
	println("="^70)

	# Standard boundary conditions
	Q_EXTERNAL = 50e3
	H_CONV = 7.2
	T_AMBIENT = 300.0

	bc_set = BoundaryConditionSet(
		thermal_top = RadiativeFluxBC(flux = Q_EXTERNAL, absorptivity = 0.9) +
					  RadiativeBC(ε = 0.9, T_inf = 330) +
					  ConvectiveBC(h = H_CONV, T_inf = 330),
		thermal_bottom = AdiabaticBC(),
		mass_top = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0),
		mass_bottom = ImpermeableBC(),
	)

	bc_set_no_radiation = BoundaryConditionSet(
		thermal_top = ConvectiveBC(h = H_CONV, T_inf = 400),
		thermal_bottom = AdiabaticBC(),
		mass_top = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0),
		mass_bottom = ImpermeableBC(),
	)

	n_passed = 0
	n_total = 0

	# ========================================================================
	# GROUP 1: Temperature Variation Tests (Diffusion-Only)
	# ========================================================================
	println("\n" * "="^70)
	println("GROUP 1: Temperature Variation Tests (Diffusion-Only Transport)")
	println("="^70)

	material_diff = create_douglas_fir_diffusion_only()

	for T in [300.0, 400.0, 500.0, 600.0, 700.0, 800.0]
		passed, _, _ = test_jacobian_configuration(
			name = "Uniform T=$(Int(T))K, Virgin",
			material = material_diff,
			bc_set = bc_set,
			state_modifier = (u0, sv, NC, mat) -> modify_state_virgin(u0, sv, NC, T, mat),
			use_ale = true,
			n_cells = 40,
		)
		n_passed += passed ? 1 : 0
		n_total += 1
	end

	# ========================================================================
	# GROUP 2: Decomposition State Variation
	# ========================================================================
	println("\n" * "="^70)
	println("GROUP 2: Decomposition State Variation")
	println("="^70)

	# Partial decomposition with temperature gradient
	passed, _, _ = test_jacobian_configuration(
		name = "Partial Decomposition (T=350-700K gradient)",
		material = material_diff,
		bc_set = bc_set,
		state_modifier = (u0, sv, NC, mat) -> modify_state_partial_decomposition(u0, sv, NC, 350.0, 700.0, mat),
		use_ale = true,
		n_cells = 40,
	)
	n_passed += passed ? 1 : 0
	n_total += 1

	# All components present
	for T in [450.0, 550.0, 650.0]
		passed, _, _ = test_jacobian_configuration(
			name = "All Components Active, T=$(Int(T))K",
			material = material_diff,
			bc_set = bc_set,
			state_modifier = (u0, sv, NC, mat) -> modify_state_all_components(u0, sv, NC, T),
			use_ale = true,
			n_cells = 30,
		)
		n_passed += passed ? 1 : 0
		n_total += 1
	end

	# ========================================================================
	# GROUP 3: High Gas Concentration Tests
	# ========================================================================
	println("\n" * "="^70)
	println("GROUP 3: High Gas Concentration (Stress Gas Transport)")
	println("="^70)

	for T in [400.0, 500.0, 600.0]
		passed, _, _ = test_jacobian_configuration(
			name = "High Gas Content, T=$(Int(T))K",
			material = material_diff,
			bc_set = bc_set,
			state_modifier = (u0, sv, NC, mat) -> modify_state_high_gas(u0, sv, NC, T),
			use_ale = true,
			n_cells = 30,
		)
		n_passed += passed ? 1 : 0
		n_total += 1
	end

	# ========================================================================
	# GROUP 4: Darcy Flow Tests
	# ========================================================================
	println("\n" * "="^70)
	println("GROUP 4: Darcy Flow Tests (Pressure-Driven Transport)")
	println("="^70)

	material_darcy = create_douglas_fir_with_darcy()

	# Virgin with Darcy
	for T in [400.0, 500.0, 600.0]
		passed, _, _ = test_jacobian_configuration(
			name = "Darcy: Virgin T=$(Int(T))K",
			material = material_darcy,
			bc_set = bc_set,
			state_modifier = (u0, sv, NC, mat) -> modify_state_virgin(u0, sv, NC, T, mat),
			use_ale = true,
			n_cells = 30,
		)
		n_passed += passed ? 1 : 0
		n_total += 1
	end

	# Darcy with gas accumulation
	for T in [450.0, 550.0]
		passed, _, _ = test_jacobian_configuration(
			name = "Darcy: High Gas T=$(Int(T))K",
			material = material_darcy,
			bc_set = bc_set,
			state_modifier = (u0, sv, NC, mat) -> modify_state_high_gas(u0, sv, NC, T),
			use_ale = true,
			n_cells = 30,
		)
		n_passed += passed ? 1 : 0
		n_total += 1
	end

	# Darcy with partial decomposition
	passed, _, _ = test_jacobian_configuration(
		name = "Darcy: Partial Decomposition (350-650K)",
		material = material_darcy,
		bc_set = bc_set,
		state_modifier = (u0, sv, NC, mat) -> modify_state_partial_decomposition(u0, sv, NC, 350.0, 650.0, mat),
		use_ale = true,
		n_cells = 40,
	)
	n_passed += passed ? 1 : 0
	n_total += 1

	# ========================================================================
	# GROUP 5: ALE Mesh Deformation Tests
	# ========================================================================
	println("\n" * "="^70)
	println("GROUP 5: ALE Mesh Deformation Tests")
	println("="^70)

	for shrinkage in [0.05, 0.15, 0.30]
		passed, _, _ = test_jacobian_configuration(
			name = "Deformed Mesh ($(Int(shrinkage*100))% shrinkage)",
			material = material_diff,
			bc_set = bc_set,
			state_modifier = (u0, sv, NC, mat) -> modify_state_partial_decomposition(u0, sv, NC, 350.0, 650.0, mat),
			mesh_modifier = (mesh, u0, sv, NC, mat) -> modify_mesh_deformed(mesh, u0, sv, NC, 350.0, 650.0, shrinkage, mat),
			use_ale = true,
			n_cells = 30,
		)
		n_passed += passed ? 1 : 0
		n_total += 1
	end

	# ========================================================================
	# GROUP 6: Beer-Lambert Radiation Tests
	# ========================================================================
	println("\n" * "="^70)
	println("GROUP 6: Beer-Lambert In-Depth Radiation")
	println("="^70)

	material_bl = create_douglas_fir_beer_lambert()

	for T in [400.0, 500.0, 600.0]
		passed, _, _ = test_jacobian_configuration(
			name = "Beer-Lambert: T=$(Int(T))K",
			material = material_bl,
			bc_set = bc_set,
			state_modifier = (u0, sv, NC, mat) -> modify_state_virgin(u0, sv, NC, T, mat),
			use_ale = true,
			radiation_model = BEER_LAMBERT,
			n_cells = 30,
		)
		n_passed += passed ? 1 : 0
		n_total += 1
	end

	# Beer-Lambert with partial decomposition
	passed, _, _ = test_jacobian_configuration(
		name = "Beer-Lambert: Partial Decomposition",
		material = material_bl,
		bc_set = bc_set,
		state_modifier = (u0, sv, NC, mat) -> modify_state_partial_decomposition(u0, sv, NC, 350.0, 650.0, mat),
		use_ale = true,
		radiation_model = BEER_LAMBERT,
		n_cells = 40,
	)
	n_passed += passed ? 1 : 0
	n_total += 1

	# ========================================================================
	# GROUP 7: Non-ALE (Fixed Mesh) Tests
	# ========================================================================
	# GROUP 7: Non-ALE (Fixed Mesh) Tests
	# ========================================================================
	println("\n" * "="^70)
	println("GROUP 7: Non-ALE (Fixed Mesh) Tests")
	println("="^70)

	# NOTE: Use NO_RADIATION here because bc_set has RadiativeFluxBC
	# SURFACE_ABSORPTION would require setting I_radiation parameter instead
	for T in [400.0, 550.0, 700.0]
		passed, _, _ = test_jacobian_configuration(
			name = "Fixed Mesh: T=$(Int(T))K",
			material = material_diff,
			bc_set = bc_set,
			state_modifier = (u0, sv, NC, mat) -> modify_state_virgin(u0, sv, NC, T, mat),
			use_ale = false,
			radiation_model = NO_RADIATION,
			n_cells = 40,
		)
		n_passed += passed ? 1 : 0
		n_total += 1
	end

	# Fixed mesh with all components
	passed, _, _ = test_jacobian_configuration(
		name = "Fixed Mesh: All Components T=500K",
		material = material_diff,
		bc_set = bc_set,
		state_modifier = (u0, sv, NC, mat) -> modify_state_all_components(u0, sv, NC, 500.0),
		use_ale = false,
		radiation_model = NO_RADIATION,
		n_cells = 30,
	)
	n_passed += passed ? 1 : 0
	n_total += 1

	# ========================================================================
	# GROUP 8: Simple Material (Focused Testing)
	# ========================================================================
	println("\n" * "="^70)
	println("GROUP 8: Simple 3-Component Material (Focused Testing)")
	println("="^70)

	material_simple = create_simple_material()

	bc_set_simple = BoundaryConditionSet(
		thermal_top = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300),
		thermal_bottom = AdiabaticBC(),
		mass_top = DiffusiveMassBC(0.01),
		mass_bottom = ImpermeableBC(),
	)

	for T in [400.0, 550.0, 700.0]
		passed, _, _ = test_jacobian_configuration(
			name = "Simple Material: T=$(Int(T))K",
			material = material_simple,
			bc_set = bc_set_simple,
			state_modifier = (u0, sv, NC, mat) -> begin
				# Simple material: virgin=1, char=2, gas=3
				n = sv.n_cells
				ρ_virgin = mat.components[1].density(T)
				ρ_char = 150.0  # char density
				ρ_gas = 5.0     # gas concentration
				for i in 1:n
					u0[i] = T
					u0[n+i] = ρ_virgin * 0.9      # mostly virgin
					u0[n+n+i] = ρ_char * 0.05   # small amount char
					u0[n+2*n+i] = ρ_gas         # some gas
				end
				return u0  # Must return u0!
			end,
			use_ale = true,
			n_cells = 30,
		)
		n_passed += passed ? 1 : 0
		n_total += 1
	end

	# ========================================================================
	# GROUP 9: FD Step Size Sensitivity
	# ========================================================================
	println("\n" * "="^70)
	println("GROUP 9: FD Step Size Sensitivity Analysis")
	println("="^70)

	test_fd_step_sensitivity(
		material_diff, bc_set,
		(u0, sv, NC, mat) -> modify_state_all_components(u0, sv, NC, 500.0),
	)

	# ========================================================================
	# SUMMARY
	# ========================================================================
	println("\n" * "="^70)
	println("FINAL SUMMARY")
	println("="^70)

	println(@sprintf("\nTotal tests run: %d", n_total))
	println(@sprintf("Tests passed:    %d (%.1f%%)", n_passed, 100.0 * n_passed / n_total))
	println(@sprintf("Tests failed:    %d (%.1f%%)", n_total - n_passed, 100.0 * (n_total - n_passed) / n_total))

	# List any failures
	failures = filter(r -> !r.passed, ALL_RESULTS)
	if !isempty(failures)
		println("\n⚠ FAILED TESTS:")
		for f in failures
			println(@sprintf("  - %s", f.name))
			println(@sprintf("      rel_err=%.2e, abs_err=%.2e, norm_fd=%.2e",
				f.rel_error, f.abs_error, f.norm_fd))
			println(@sprintf("      worst: J[%d,%d] An=%.4e, FD=%.4e",
				f.worst_element[1], f.worst_element[2], f.worst_an, f.worst_fd))
		end
	else
		println("\n✓ ALL TESTS PASSED!")
	end

	println("\n" * "="^70)
	println("Test suite complete!")
	println("="^70)

	return n_passed == n_total
end

# ============================================================================
# Run Tests
# ============================================================================

if abspath(PROGRAM_FILE) == @__FILE__
	success = run_all_tests()
	exit(success ? 0 : 1)
end
