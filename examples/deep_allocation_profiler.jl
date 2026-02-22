"""
Deep Allocation Profiler for Douglas Fir Pyrolysis Simulation

This script performs an exhaustive, hierarchical profiling to pinpoint
EXACTLY where the massive allocations (410 GB for 50s simulation!) are occurring.

The approach is:
1. Start from the top (ODE solve) and measure at each level
2. Drill down into the component that allocates
3. Continue until we find the allocating function(s)

Key insight: With 40 cells × 14 DOFs = 560 DOFs, a dense Jacobian evaluation
requires 560 RHS calls per Jacobian. If each RHS call allocates even 100 bytes,
that's 56 KB per Jacobian. With ~millions of Jacobian evaluations, this explodes.

Usage:
	julia --project=. examples/deep_allocation_profiler.jl
"""

using Pyrolysis
using Pyrolysis.Core
using Pyrolysis.Core: R_GAS, n_components, n_reactions, has_darcy_flow, Phase, SOLID, LIQUID, GAS
using Pyrolysis.Properties
using Pyrolysis.Properties: PropertyCache, update_properties!, effective_permeability,
	volume_fractions, effective_conductivity, effective_heat_capacity,
	effective_emissivity, mixture_density, effective_gas_transfer,
	effective_absorption, harmonic_mean_weighted, harmonic_mean_product_weighted
using Pyrolysis.Physics
using Pyrolysis.Physics: ideal_gas_pressure, evaluate_reactions!, reaction_heat_source,
	species_source_terms!, compute_conductive_fluxes!, compute_gas_fluxes!,
	compute_radiation_source!, compute_convective_heat_fluxes_from_mass_fluxes!
using Pyrolysis.Discretization
using Pyrolysis.Discretization:
	compute_rhs!, compute_rhs_darcy_fick!, RHSCache,
	update_darcy_properties!, boundary_face_ids, divergence_1d!, divergence_species_1d!,
	apply_thermal_bc_flux!, apply_mass_bc_flux!, compute_ale_advection_term!,
	compute_ale_species_advection!
using Pyrolysis.TimeIntegration
using Pyrolysis.TimeIntegration: PyrolysisODEFunctionALE, StateVectorALE,
	create_state_vector_ale, unpack_state_ale!, JacobianCacheALE,
	compute_analytical_jacobian_ale!, jacobian_sparsity_ale
using Pyrolysis.Geometry: GeometryCache
using Printf
using InteractiveUtils
using SparseArrays
using OrdinaryDiffEq

# ============================================================================
# Configuration
# ============================================================================
const PROFILE_N_CELLS = 40
const N_ITERS = 100        # Iterations for timing
const N_ITERS_FAST = 1000  # For very fast functions

println("=" ^ 80)
println("DEEP ALLOCATION PROFILER - Douglas Fir Pyrolysis")
println("=" ^ 80)
println()

# ============================================================================
# Material Definition (same as main example)
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

# State-dependent heat of sorption
heat_of_sorption = StateDependentProperty((T, ξ) -> 2673.82e3)

# Gas transfer coefficients
λ_char = 5e-5
λ_wood = 1e-5
λ_gas = 1e-4
κ_virgin = 1e-15
κ_char = 1e-12

const douglas_fir = Material(
	name = :DouglasFir,
	components = (
		SolidComponent(:DF, ρ = 530.0, c = df_cp, k = 0.12, ε = 0.9, α = 0.0, λ = λ_wood, κ = κ_virgin),
		SolidComponent(:DF_int1, ρ = 507.2, c = df_int1_cp, k = 0.14, ε = 0.9, α = 0.0, λ = λ_wood*1.2, κ = 5e-15),
		SolidComponent(:DF_int2, ρ = 419.8, c = df_int2_cp, k = 0.15, ε = 0.9, α = 0.0, λ = λ_wood*1.5, κ = 1e-14),
		SolidComponent(:DF_int3, ρ = 203.2, c = df_int3_cp, k = 0.16, ε = 0.9, α = 0.0, λ = λ_wood*2.0, κ = 1e-13),
		SolidComponent(:DF_int4, ρ = 176.6, c = df_int4_cp, k = 0.18, ε = 0.9, α = 0.0, λ = λ_char*0.8, κ = 5e-13),
		SolidComponent(:Char, ρ = 150.0, c = char_cp, k = 0.20, ε = 0.95, α = 0.0, λ = λ_char, κ = κ_char),
		LiquidComponent(:MOISTURE, ρ = 1000.0, c = moisture_cp, k = 0.6, ε = 0.9, α = 0.0, κ = 0.0),
		GaseousComponent(:DF_gas1, M = 0.060, c = 2000.0, k = 0.03, λ = λ_gas),
		GaseousComponent(:DF_gas2, M = 0.044, c = 1800.0, k = 0.03, λ = λ_gas),
		GaseousComponent(:DF_gas3, M = 0.070, c = 2200.0, k = 0.025, λ = λ_gas*0.8),
		GaseousComponent(:DF_gas4, M = 0.030, c = 1600.0, k = 0.035, λ = λ_gas),
		GaseousComponent(:DF_gas5, M = 0.028, c = 1500.0, k = 0.035, λ = λ_gas),
		GaseousComponent(:WATER_VAPOR, M = 0.018, c = water_vapor_cp, k = 0.04, λ = λ_gas*1.2),
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

const NC = 13
const NR = 6

# ============================================================================
# Create Test Infrastructure
# ============================================================================

println("Creating test infrastructure...")

sample_thickness = 12.0e-3
T_ambient = 300.0
q_external = 50e3
moisture_content_db = 0.15
ρ_DF_initial = 530.0
ρ_moisture_initial = moisture_content_db * ρ_DF_initial

# Create mesh with some reacted material (more realistic state)
mesh = create_uniform_mesh(
	sample_thickness, PROFILE_N_CELLS, Val(NC);
	T_initial = 500.0,  # Hot enough for reactions
	ξ_initial = (400.0, 50.0, 30.0, 20.0, 10.0, 10.0, 30.0,
		5.0, 3.0, 3.0, 2.0, 2.0, 3.0),
)

n_cells = length(mesh.cells)
n_faces = length(mesh.faces)
n_nodes = length(mesh.nodes)

# Boundary conditions
top_thermal_bc = RadiativeFluxBC(flux = q_external, absorptivity = 0.9) +
				 RadiativeBC(ε = 0.9, T_inf = 330) +
				 ConvectiveBC(h = 7.2, T_inf = 330)
top_mass_bc = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0)
bottom_thermal_bc = AdiabaticBC()
bottom_mass_bc = ImpermeableBC()

bc_set = BoundaryConditionSet(
	thermal_top = top_thermal_bc,
	thermal_bottom = bottom_thermal_bc,
	mass_top = top_mass_bc,
	mass_bottom = bottom_mass_bc,
)

# Pre-cache boundary face IDs
top_face_ids = boundary_face_ids(mesh, 1)
bottom_face_ids = boundary_face_ids(mesh, 2)

# Create caches
prop_cache = PropertyCache(n_cells, n_faces)
rhs_cache = RHSCache(n_cells, n_faces, NC, NR)

# Output arrays
dT = zeros(n_cells)
dξ = zeros(NC, n_cells)

# Create ALE ODE function
ode_fn = PyrolysisODEFunctionALE(mesh, douglas_fir, bc_set;
	I_radiation = q_external * 0.9, radiation_model = SURFACE_ABSORPTION,
	min_thickness = 0.5e-3, min_density = 1.0)

# Create state vector
u0 = create_state_vector_ale(mesh, douglas_fir)
du = similar(u0)

# Create Jacobian cache
jac_cache = JacobianCacheALE(n_cells, n_faces, n_nodes, NC, NR)
jac_prototype = jacobian_sparsity_ale(n_cells, NC)
J = copy(jac_prototype)

println("  Mesh: $n_cells cells, $n_faces faces, $n_nodes nodes")
println("  DOFs: $(length(u0)) = $n_cells × (1 + $NC) + $n_nodes")
println("  Jacobian nnz: $(nnz(J))")
println()

# ============================================================================
# Helper Functions
# ============================================================================

function print_alloc(name::String, allocs::Int; indent::Int = 0, threshold::Int = 0)
	prefix = "  " ^ indent
	status = allocs <= threshold ? "✓" : "✗"
	color = allocs <= threshold ? "\e[32m" : "\e[31m"
	reset = "\e[0m"
	if allocs > 0
		@printf("%s%-50s %12d bytes  %s%s%s\n", prefix, name, allocs, color, status, reset)
	else
		@printf("%s%-50s %12d bytes  %s%s%s\n", prefix, name, allocs, color, status, reset)
	end
	return allocs
end

function test_alloc(f, name::String; n_warmup::Int = 2, indent::Int = 0, threshold::Int = 0)
	# Warmup
	for _ in 1:n_warmup
		f()
	end
	# Measure
	allocs = @allocated f()
	print_alloc(name, allocs; indent = indent, threshold = threshold)
	return allocs
end

# ============================================================================
# LEVEL 1: Top-Level Functions
# ============================================================================

println("=" ^ 80)
println("LEVEL 1: Top-Level Functions")
println("=" ^ 80)
println()
println("Testing the main entry points that are called by DifferentialEquations.jl")
println()

# 1.1 Full ODE RHS evaluation
println("-" ^ 60)
println("1.1 PyrolysisODEFunctionALE (du, u, p, t) - Full RHS")
println("-" ^ 60)

allocs_ode = test_alloc("ode_fn(du, u0, nothing, 0.0)", indent = 1) do
	ode_fn(du, u0, nothing, 0.0)
end

if allocs_ode > 0
	println("\n  ⚠ ODE function allocates! This is called every RHS evaluation.")
	println("    With stiff solvers, this can be called millions of times.")
	println("    Estimated impact: $(allocs_ode) bytes × ~1M calls = $(allocs_ode / 1e6) GB")
end
println()

# 1.2 Analytical Jacobian evaluation
println("-" ^ 60)
println("1.2 compute_analytical_jacobian_ale! - Full Jacobian")
println("-" ^ 60)

allocs_jac = test_alloc("compute_analytical_jacobian_ale!(J, u0, nothing, 0.0, ...)", indent = 1) do
	compute_analytical_jacobian_ale!(J, u0, nothing, 0.0, ode_fn, jac_cache)
end

if allocs_jac > 0
	println("\n  ⚠ Jacobian computation allocates! Called ~thousands of times.")
	println("    Estimated impact: $(allocs_jac) bytes × ~10K calls = $(allocs_jac / 1e3) MB")
end
println()

# ============================================================================
# LEVEL 2: ODE Function Components
# ============================================================================

println("=" ^ 80)
println("LEVEL 2: ODE Function Internal Components")
println("=" ^ 80)
println()
println("Breaking down the ODE function to find where allocations occur...")
println()

# 2.1 State unpacking
println("-" ^ 60)
println("2.1 State Unpacking")
println("-" ^ 60)

allocs_unpack = test_alloc("unpack_state_ale!(mesh, u0, sv)", indent = 1) do
	unpack_state_ale!(mesh, u0, ode_fn.sv)
end
println()

# 2.2 Property update
println("-" ^ 60)
println("2.2 Property Cache Update")
println("-" ^ 60)

allocs_props = test_alloc("update_properties!(prop_cache, mesh, material)", indent = 1) do
	Properties.update_properties!(prop_cache, mesh, douglas_fir)
end

if allocs_props > 0
	println("\n  ⚠ Property update allocates! Drilling down...")
end
println()

# 2.3 Full RHS computation (compute_rhs_darcy_fick!)
println("-" ^ 60)
println("2.3 RHS Computation (Darcy-Fick)")
println("-" ^ 60)

# First update properties
Properties.update_properties!(prop_cache, mesh, douglas_fir)

allocs_rhs = test_alloc("compute_rhs_darcy_fick!(dT, dξ, ...)", indent = 1) do
	compute_rhs_darcy_fick!(dT, dξ, mesh, prop_cache, rhs_cache,
		douglas_fir, bc_set, 0.0, top_face_ids, bottom_face_ids;
		I_radiation = q_external * 0.9, radiation_model = SURFACE_ABSORPTION)
end

if allocs_rhs > 0
	println("\n  ⚠ RHS computation allocates! This is a critical hot path.")
end
println()

# 2.4 Mesh velocity computation
println("-" ^ 60)
println("2.4 ALE Mesh Velocity Computation")
println("-" ^ 60)

# This tests the private function that computes mesh velocities
# We need to call it through the ODE function's internal method
w_work = zeros(n_nodes)
θ_work = zeros(n_cells)
rates_work = zeros(NR)
heats_work = zeros(NR)
dξ_dt_work = zeros(NC)

# Warmup and test
allocs_mesh_vel = test_alloc("_compute_ale_mesh_velocity!(...)", indent = 1) do
	Pyrolysis.TimeIntegration._compute_ale_mesh_velocity!(w_work, mesh, prop_cache, douglas_fir,
		0.0, 1.0;
		θ_out = θ_work,
		precomputed_dξ_rxn = rhs_cache.dξ_rxn,
		rates = rates_work,
		heats = heats_work,
		dξ_dt = dξ_dt_work)
end
println()

# 2.5 ALE Advection Terms
println("-" ^ 60)
println("2.5 ALE Advection Terms")
println("-" ^ 60)

ale_advection = zeros(n_cells)
ale_species_advection = zeros(NC, n_cells)

allocs_ale_adv = test_alloc("compute_ale_advection_term!(...)", indent = 1) do
	compute_ale_advection_term!(ale_advection, mesh, w_work)
end

allocs_ale_species = test_alloc("compute_ale_species_advection!(...)", indent = 1) do
	compute_ale_species_advection!(ale_species_advection, mesh, douglas_fir, w_work)
end
println()

# ============================================================================
# LEVEL 3: Property Evaluation Details
# ============================================================================

println("=" ^ 80)
println("LEVEL 3: Property Evaluation Details")
println("=" ^ 80)
println()

# Get a test cell state
test_state = mesh.cell_states[1]
test_ξ = test_state.ξ
test_T = test_state.T

println("-" ^ 60)
println("3.1 Per-Cell Property Functions")
println("-" ^ 60)

# Test each property function
allocs_vf = test_alloc("volume_fractions(ξ, T, material)", indent = 1) do
	volume_fractions(test_ξ, test_T, douglas_fir)
end

allocs_k = test_alloc("effective_conductivity(ξ, T, material)", indent = 1) do
	effective_conductivity(test_ξ, test_T, douglas_fir)
end

allocs_cp = test_alloc("effective_heat_capacity(ξ, T, material)", indent = 1) do
	effective_heat_capacity(test_ξ, test_T, douglas_fir)
end

allocs_rho = test_alloc("mixture_density(ξ, T, material)", indent = 1) do
	mixture_density(test_ξ, test_T, douglas_fir)
end

allocs_lambda = test_alloc("effective_gas_transfer(ξ, T, material)", indent = 1) do
	effective_gas_transfer(test_ξ, test_T, douglas_fir)
end

allocs_eps = test_alloc("effective_emissivity(ξ, T, material)", indent = 1) do
	effective_emissivity(test_ξ, test_T, douglas_fir)
end

allocs_kappa = test_alloc("effective_permeability(ξ, T, material)", indent = 1) do
	effective_permeability(test_ξ, test_T, douglas_fir)
end

allocs_P = test_alloc("ideal_gas_pressure(ξ, T, material)", indent = 1) do
	ideal_gas_pressure(test_ξ, test_T, douglas_fir)
end
println()

# ============================================================================
# LEVEL 3B: CallableProperty Allocations
# ============================================================================

println("-" ^ 60)
println("3.2 CallableProperty Evaluation (Temperature-Dependent cp)")
println("-" ^ 60)

# Test the CallableProperty functions directly
# These are used for heat capacities in Douglas Fir

allocs_cp_call_1 = test_alloc("df_cp(T) - CallableProperty", indent = 1) do
	df_cp(test_T)
end

allocs_cp_call_2 = test_alloc("moisture_cp(T) - CallableProperty", indent = 1) do
	moisture_cp(test_T)
end

allocs_cp_call_3 = test_alloc("water_vapor_cp(T) - CallableProperty", indent = 1) do
	water_vapor_cp(test_T)
end

# Test constant property for comparison
const_prop = ConstantProperty(1000.0)
allocs_const = test_alloc("ConstantProperty(T)", indent = 1) do
	const_prop(test_T)
end

# Test linear property
lin_prop = LinearProperty(1000.0, 2.0)
allocs_lin = test_alloc("LinearProperty(T)", indent = 1) do
	lin_prop(test_T)
end
println()

# ============================================================================
# LEVEL 3C: Component Property Access
# ============================================================================

println("-" ^ 60)
println("3.3 Component Property Access")
println("-" ^ 60)

# Test accessing component properties
comp1 = douglas_fir.components[1]
allocs_comp_cp = test_alloc("component.heat_capacity(T)", indent = 1) do
	comp1.heat_capacity(test_T)
end

allocs_comp_k = test_alloc("component.thermal_conductivity(T)", indent = 1) do
	comp1.thermal_conductivity(test_T)
end

# Gas component
gas_comp = douglas_fir.components[8]  # DF_gas1
allocs_gas_cp = test_alloc("gas_component.heat_capacity(T)", indent = 1) do
	gas_comp.heat_capacity(test_T)
end
println()

# ============================================================================
# LEVEL 4: Reaction Evaluation Details
# ============================================================================

println("=" ^ 80)
println("LEVEL 4: Reaction Evaluation")
println("=" ^ 80)
println()

rates = zeros(NR)
heats = zeros(NR)

println("-" ^ 60)
println("4.1 Full Reaction Evaluation")
println("-" ^ 60)

allocs_rxn = test_alloc("evaluate_reactions!(rates, heats, ξ, T, material)", indent = 1) do
	evaluate_reactions!(rates, heats, test_ξ, test_T, douglas_fir)
end
println()

println("-" ^ 60)
println("4.2 Individual Reaction Rate")
println("-" ^ 60)

rxn1 = douglas_fir.reactions[1]
allocs_rxn_rate = test_alloc("reaction_rate(rxn, ξ, T)", indent = 1) do
	Physics.reaction_rate(rxn1, test_ξ, test_T)
end
println()

println("-" ^ 60)
println("4.3 Reaction Heat Evaluation (State-Dependent)")
println("-" ^ 60)

# This is the one with StateDependentProperty for heat_of_sorption
rxn6 = douglas_fir.reactions[6]  # moisture_evaporation
allocs_rxn6_heat = test_alloc("rxn.heat(T, ξ) - StateDependentProperty", indent = 1) do
	rxn6.heat(test_T, test_ξ)
end

# Compare with constant heat
rxn1 = douglas_fir.reactions[1]
allocs_rxn1_heat = test_alloc("rxn.heat(T, ξ) - ConstantProperty", indent = 1) do
	rxn1.heat(test_T, test_ξ)
end
println()

# ============================================================================
# LEVEL 5: Physics Flux Computations
# ============================================================================

println("=" ^ 80)
println("LEVEL 5: Physics Flux Computations")
println("=" ^ 80)
println()

# Update properties first
Properties.update_properties!(prop_cache, mesh, douglas_fir)

println("-" ^ 60)
println("5.1 Conductive Heat Fluxes")
println("-" ^ 60)

allocs_cond = test_alloc("compute_conductive_fluxes!(q_cond, mesh, prop_cache)", indent = 1) do
	compute_conductive_fluxes!(rhs_cache.q_cond, mesh, prop_cache)
end
println()

println("-" ^ 60)
println("5.2 Gas Transport Fluxes")
println("-" ^ 60)

allocs_gas = test_alloc("compute_gas_fluxes!(J_gas, mesh, prop_cache, material)", indent = 1) do
	compute_gas_fluxes!(rhs_cache.J_gas, mesh, prop_cache, douglas_fir)
end
println()

println("-" ^ 60)
println("5.3 Radiation Source")
println("-" ^ 60)

allocs_rad = test_alloc("compute_radiation_source!(Q_rad, mesh, ...)", indent = 1) do
	compute_radiation_source!(rhs_cache.Q_rad, mesh, prop_cache, q_external * 0.9, SURFACE_ABSORPTION)
end
println()

println("-" ^ 60)
println("5.4 Convective Heat Flux from Mass Fluxes")
println("-" ^ 60)

allocs_conv = test_alloc("compute_convective_heat_fluxes_from_mass_fluxes!(...)", indent = 1) do
	compute_convective_heat_fluxes_from_mass_fluxes!(rhs_cache.q_conv, rhs_cache.J_gas, mesh, douglas_fir)
end
println()

# ============================================================================
# LEVEL 6: Boundary Condition Evaluation
# ============================================================================

println("=" ^ 80)
println("LEVEL 6: Boundary Condition Evaluation")
println("=" ^ 80)
println()

println("-" ^ 60)
println("6.1 Thermal BC Flux Application")
println("-" ^ 60)

q_total = zeros(n_faces)
allocs_thermal_bc = test_alloc("apply_thermal_bc_flux!(...)", indent = 1) do
	apply_thermal_bc_flux!(q_total, mesh, bc_set, prop_cache, 0.0,
		top_face_ids, bottom_face_ids, SURFACE_ABSORPTION)
end
println()

println("-" ^ 60)
println("6.2 Mass BC Flux Application")
println("-" ^ 60)

J_gas = zeros(NC, n_faces)
allocs_mass_bc = test_alloc("apply_mass_bc_flux!(...)", indent = 1) do
	apply_mass_bc_flux!(J_gas, mesh, prop_cache, bc_set, 0.0,
		top_face_ids, bottom_face_ids, Val(NC))
end
println()

# ============================================================================
# LEVEL 7: Divergence Operations
# ============================================================================

println("=" ^ 80)
println("LEVEL 7: Divergence Operations")
println("=" ^ 80)
println()

div_q = zeros(n_cells)
div_J = zeros(NC, n_cells)

allocs_div_q = test_alloc("divergence_1d!(div_q, fluxes, mesh)", indent = 1) do
	divergence_1d!(div_q, rhs_cache.q_total, mesh)
end

allocs_div_J = test_alloc("divergence_species_1d!(div_J, J, mesh)", indent = 1) do
	divergence_species_1d!(div_J, rhs_cache.J_gas, mesh)
end
println()

# ============================================================================
# LEVEL 8: Jacobian Components
# ============================================================================

println("=" ^ 80)
println("LEVEL 8: Jacobian Computation Details")
println("=" ^ 80)
println()

println("-" ^ 60)
println("8.1 Sparsity Map Creation (once per solve)")
println("-" ^ 60)

allocs_sparsity = test_alloc("FastSparsityMap(J) or SparsityMap(J)", indent = 1) do
	Pyrolysis.TimeIntegration.FastSparsityMap(J)
end
println()

# ============================================================================
# TYPE STABILITY CHECKS
# ============================================================================

println("=" ^ 80)
println("TYPE STABILITY CHECKS")
println("=" ^ 80)
println()
println("Running @code_warntype on critical functions...")
println("(Red/yellow text indicates type instability → allocations)")
println()

println("-" ^ 60)
println("evaluate_reactions! type stability:")
println("-" ^ 60)
@code_warntype evaluate_reactions!(rates, heats, test_ξ, test_T, douglas_fir)
println()

println("-" ^ 60)
println("effective_heat_capacity type stability:")
println("-" ^ 60)
@code_warntype effective_heat_capacity(test_ξ, test_T, douglas_fir)
println()

println("-" ^ 60)
println("CallableProperty call type stability:")
println("-" ^ 60)
@code_warntype df_cp(test_T)
println()

# ============================================================================
# SUMMARY
# ============================================================================

println()
println("=" ^ 80)
println("ALLOCATION SUMMARY")
println("=" ^ 80)
println()

total_per_rhs = allocs_ode
total_per_jac = allocs_jac

println("Per-RHS-call allocations:       $(total_per_rhs) bytes")
println("Per-Jacobian-call allocations:  $(total_per_jac) bytes")
println()

# Estimate total allocations for a solve
# Typical stiff solve: ~1000 steps, each step needs ~10 Jacobians, ~100 RHS evals
est_steps = 1000
est_jac_per_step = 10
est_rhs_per_step = 100
est_total_rhs = est_steps * est_rhs_per_step
est_total_jac = est_steps * est_jac_per_step

est_alloc_rhs = total_per_rhs * est_total_rhs
est_alloc_jac = total_per_jac * est_total_jac
est_total = est_alloc_rhs + est_alloc_jac

println("Estimated allocations for $(est_steps)-step solve:")
println("  RHS calls: $(est_total_rhs) × $(total_per_rhs) bytes = $(est_alloc_rhs / 1e9) GB")
println("  Jacobian calls: $(est_total_jac) × $(total_per_jac) bytes = $(est_alloc_jac / 1e9) GB")
println("  TOTAL: $(est_total / 1e9) GB")
println()

if total_per_rhs > 0
	println("⚠ CRITICAL: ODE RHS function allocates $(total_per_rhs) bytes per call")
	println("   This is the PRIMARY cause of performance issues!")
	println()
	println("   Functions that allocate (from above analysis):")
	allocs_ode > 0 && println("     - Full ODE function: $(allocs_ode) bytes")
	allocs_props > 0 && println("     - Property update: $(allocs_props) bytes")
	allocs_rhs > 0 && println("     - RHS computation: $(allocs_rhs) bytes")
	allocs_mesh_vel > 0 && println("     - Mesh velocity: $(allocs_mesh_vel) bytes")
end

println()
println("=" ^ 80)
println("LEVEL 9: Mixing Rule Functions")
println("=" ^ 80)
println()
println("Testing allocation-free mixing rule implementations...")
println()

println("-" ^ 60)
println("9.1 Harmonic Mean Mixing Rules (2-value)")
println("-" ^ 60)

# Test harmonic mean for face properties (uses 2 values)
k1 = douglas_fir.components[1].thermal_conductivity(test_T)
k2 = douglas_fir.components[2].thermal_conductivity(test_T)
d1 = 0.1
d2 = 0.2

allocs_harmonic_weighted = test_alloc("harmonic_mean_weighted(k1, k2, d1, d2)", indent = 1) do
	harmonic_mean_weighted(k1, k2, d1, d2)
end

vf1 = 0.5
vf2 = 0.5
allocs_harmonic_product = test_alloc("harmonic_mean_product_weighted(k1, k2, d1, d2, vf1, vf2)", indent = 1) do
	harmonic_mean_product_weighted(k1, k2, d1, d2, vf1, vf2)
end
println()

# ============================================================================
# LEVEL 10: Species Source Terms
# ============================================================================

println("=" ^ 80)
println("LEVEL 10: Species Source Terms")
println("=" ^ 80)
println()

dξ_dt = zeros(NC)
allocs_species_source = test_alloc("species_source_terms!(dξ_dt, rates, material, Val(NC))", indent = 1) do
	species_source_terms!(dξ_dt, rates, douglas_fir, Val(NC))
end
println()

# ============================================================================
# LEVEL 11: Darcy Flow Properties
# ============================================================================

println("=" ^ 80)
println("LEVEL 11: Darcy Flow Property Updates")
println("=" ^ 80)
println()

# Test Darcy-specific property updates
println("-" ^ 60)
println("11.1 Darcy Property Cache Update")
println("-" ^ 60)

allocs_darcy_props = test_alloc("update_darcy_properties!(prop_cache, mesh, material, bc_set, ...)", indent = 1) do
	update_darcy_properties!(prop_cache, mesh, douglas_fir, bc_set, top_face_ids, bottom_face_ids)
end
println()

# ============================================================================
# LEVEL 12: Full Component Property Loops
# ============================================================================

println("=" ^ 80)
println("LEVEL 12: Component Iteration Patterns")
println("=" ^ 80)
println()
println("Testing loops over all components to check for allocation patterns...")
println()

println("-" ^ 60)
println("12.1 Heat Capacity Sum Loop")
println("-" ^ 60)

allocs_cp_loop = test_alloc("sum over component heat capacities", indent = 1) do
	sum = 0.0
	for i in 1:NC
		comp = douglas_fir.components[i]
		sum += test_ξ[i] * comp.heat_capacity(test_T)
	end
	sum
end

println("-" ^ 60)
println("12.2 Volume Fraction Calculation Loop")
println("-" ^ 60)

allocs_vf_loop = test_alloc("volume fraction calculation loop", indent = 1) do
	vf = zeros(NC)
	for i in 1:NC
		comp = douglas_fir.components[i]
		if comp.phase == SOLID || comp.phase == LIQUID
			vf[i] = test_ξ[i] / comp.density(test_T)
		else
			vf[i] = 0.0
		end
	end
	vf
end
println()

# ============================================================================
# LEVEL 13: Mesh Access Patterns
# ============================================================================

println("=" ^ 80)
println("LEVEL 13: Mesh Data Access Patterns")
println("=" ^ 80)
println()
println("Testing mesh data structure access for allocations...")
println()

println("-" ^ 60)
println("13.1 Cell State Access")
println("-" ^ 60)

allocs_cell_access = test_alloc("mesh.cell_states[i] access", indent = 1, n_warmup = 5) do
	state = mesh.cell_states[1]
	state.T + sum(state.ξ)
end

println("-" ^ 60)
println("13.2 Face Neighbor Access")
println("-" ^ 60)

allocs_face_neighbors = test_alloc("mesh.faces[i].cell_ids access", indent = 1, n_warmup = 5) do
	face = mesh.faces[1]
	idx_L, idx_R = face.cell_ids
	idx_L + idx_R
end

println("-" ^ 60)
println("13.3 Cell Volume Access Loop")
println("-" ^ 60)

allocs_volume_loop = test_alloc("loop over cell volumes", indent = 1) do
	total = 0.0
	for i in 1:length(mesh.cells)
		total += mesh.cell_volumes[i]
	end
	total
end
println()

# ============================================================================
# LEVEL 14: Cache Array Access
# ============================================================================

println("=" ^ 80)
println("LEVEL 14: Cache Array Access Patterns")
println("=" ^ 80)
println()

println("-" ^ 60)
println("14.1 Property Cache Array Access")
println("-" ^ 60)

allocs_prop_array = test_alloc("prop_cache.k[i] access", indent = 1, n_warmup = 5) do
	sum(prop_cache.k)
end

println("-" ^ 60)
println("14.2 RHS Cache Matrix Access")
println("-" ^ 60)

allocs_rhs_matrix = test_alloc("rhs_cache.dξ_rxn[:, i] access", indent = 1, n_warmup = 5) do
	sum(rhs_cache.dξ_rxn)
end
println()

# ============================================================================
# LEVEL 15: Material Tuple Access
# ============================================================================

println("=" ^ 80)
println("LEVEL 15: Material Tuple Access")
println("=" ^ 80)
println()
println("Testing material.components and material.reactions tuple access...")
println()

println("-" ^ 60)
println("15.1 Component Tuple Iteration")
println("-" ^ 60)

allocs_comp_iter = test_alloc("for comp in material.components", indent = 1) do
	count = 0
	for comp in douglas_fir.components
		count += 1
	end
	count
end

println("-" ^ 60)
println("15.2 Reaction Tuple Iteration")
println("-" ^ 60)

allocs_rxn_iter = test_alloc("for rxn in material.reactions", indent = 1) do
	count = 0
	for rxn in douglas_fir.reactions
		count += 1
	end
	count
end

println("-" ^ 60)
println("15.3 Component Indexing")
println("-" ^ 60)

allocs_comp_index = test_alloc("material.components[i] in loop", indent = 1) do
	sum_ρ = 0.0
	for i in 1:NC
		comp = douglas_fir.components[i]
		sum_ρ += comp.density(test_T)
	end
	sum_ρ
end
println()

# ============================================================================
# LEVEL 16: Callable/StateDependentProperty Deep Dive
# ============================================================================

println("=" ^ 80)
println("LEVEL 16: Property Function Call Deep Dive")
println("=" ^ 80)
println()
println("Deep analysis of property function evaluation mechanisms...")
println()

println("-" ^ 60)
println("16.1 Multiple CallableProperty Calls")
println("-" ^ 60)

allocs_multi_cp = test_alloc("10 CallableProperty evaluations", indent = 1) do
	sum = 0.0
	for _ in 1:10
		sum += df_cp(test_T)
	end
	sum
end

println("-" ^ 60)
println("16.2 StateDependentProperty with State")
println("-" ^ 60)

allocs_state_dep = test_alloc("heat_of_sorption(T, ξ)", indent = 1) do
	heat_of_sorption(test_T, test_ξ)
end

println("-" ^ 60)
println("16.3 Mixed Property Calls")
println("-" ^ 60)

allocs_mixed_props = test_alloc("const + linear + callable properties", indent = 1) do
	c1 = const_prop(test_T)
	c2 = lin_prop(test_T)
	c3 = df_cp(test_T)
	c1 + c2 + c3
end
println()

# ============================================================================
# LEVEL 17: Array Operations in Physics
# ============================================================================

println("=" ^ 80)
println("LEVEL 17: Array Operation Patterns")
println("=" ^ 80)
println()

println("-" ^ 60)
println("17.1 In-Place Array Operations")
println("-" ^ 60)

test_arr = zeros(n_cells)
test_vals = rand(n_cells)

allocs_inplace_add = test_alloc("arr .+= vals in-place", indent = 1) do
	test_arr .+= test_vals
end

allocs_inplace_mul = test_alloc("arr .*= scalar in-place", indent = 1) do
	test_arr .*= 2.0
end

println("-" ^ 60)
println("17.2 View vs Copy Operations")
println("-" ^ 60)

test_matrix = zeros(NC, n_cells)

allocs_view = test_alloc("@view matrix[:, i]", indent = 1, n_warmup = 5) do
	v = @view test_matrix[:, 1]
	sum(v)
end

allocs_copy = test_alloc("matrix[:, i] (copy)", indent = 1, n_warmup = 5) do
	v = test_matrix[:, 1]
	sum(v)
end
println()

println()

# ============================================================================
# LEVEL 18: Profile-Guided Optimization
# ============================================================================

println("=" ^ 80)
println("LEVEL 18: Profile-Guided Optimization - CPU Time Hotspots")
println("=" ^ 80)
println()
println("Profiling execution time to identify computational bottlenecks...")
println("This complements allocation profiling by finding slow functions.")
println()

# Profile configuration
const PROFILE_SAMPLES = 1000  # Number of RHS evaluations to profile
const PROFILE_WARMUP = 10     # Warmup iterations

println("-" ^ 60)
println("18.1 Full RHS Execution Timing")
println("-" ^ 60)

# Warmup
for _ in 1:PROFILE_WARMUP
	ode_fn(du, u0, nothing, 0.0)
end

# Timing test
using BenchmarkTools
println("\nBenchmarking RHS evaluation:")
rhs_bench = @benchmark $ode_fn($du, $u0, nothing, 0.0)
println("  Median time: $(median(rhs_bench.times) / 1e6) ms")
println("  Min time:    $(minimum(rhs_bench.times) / 1e6) ms")
println("  Max time:    $(maximum(rhs_bench.times) / 1e6) ms")
println("  Allocations: $(rhs_bench.allocs)")
println("  Memory:      $(rhs_bench.memory) bytes")
println()

println("-" ^ 60)
println("18.2 CPU Profiling with Profile.jl")
println("-" ^ 60)

# Profile the RHS function
using Profile

# Clear any previous profiling data
Profile.clear()

# Run with profiling
println("\nCollecting profile samples ($(PROFILE_SAMPLES) iterations)...")
Profile.@profile for _ in 1:PROFILE_SAMPLES
	ode_fn(du, u0, nothing, 0.0)
end

println("Profile data collected: $(Profile.len_data()) samples")
println()

# Analyze profile data
println("Top functions by execution time:")
println("(Showing functions taking >1% of total time)")
println()

# Print flat profile
Profile.print(
	format = :flat,
	sortedby = :count,
	noisefloor = 2.0,  # Filter out functions with <2% of samples
	mincount = div(PROFILE_SAMPLES, 50)  # At least 2% of iterations
)

println()
println("-" ^ 60)
println("18.3 Per-Component Timing Breakdown")
println("-" ^ 60)

# Time individual components
components = [
	("State unpacking", () -> unpack_state_ale!(mesh, u0, ode_fn.sv)),
	("Property update", () -> Properties.update_properties!(prop_cache, mesh, douglas_fir)),
	("Darcy properties", () -> update_darcy_properties!(prop_cache, mesh, douglas_fir, bc_set, top_face_ids, bottom_face_ids)),
	("RHS computation", () -> compute_rhs_darcy_fick!(dT, dξ, mesh, prop_cache, rhs_cache,
		douglas_fir, bc_set, 0.0, top_face_ids, bottom_face_ids;
		I_radiation = q_external * 0.9, radiation_model = SURFACE_ABSORPTION)),
	("Mesh velocity", () -> Pyrolysis.TimeIntegration._compute_ale_mesh_velocity!(w_work, mesh, prop_cache, douglas_fir,
		0.0, 1.0; θ_out = θ_work, precomputed_dξ_rxn = rhs_cache.dξ_rxn,
		rates = rates_work, heats = heats_work, dξ_dt = dξ_dt_work)),
]

# Ensure properties are updated first
Properties.update_properties!(prop_cache, mesh, douglas_fir)
update_darcy_properties!(prop_cache, mesh, douglas_fir, bc_set, top_face_ids, bottom_face_ids)

timings = []
for (name, func) in components
	# Warmup
	for _ in 1:5
		func()
	end
	# Benchmark
	b = @benchmark $func()
	push!(timings, (name, median(b.times) / 1e6))
	println("  $(name): $(median(b.times) / 1e6) ms")
end

println()
total_time = sum(t[2] for t in timings)
println("Total component time: $(total_time) ms")
println()
println("Percentage breakdown:")
for (name, time) in timings
	pct = 100.0 * time / total_time
	bars = "█" ^ round(Int, pct / 2)
	@printf("  %-25s %5.1f%%  %s\n", name, pct, bars)
end
println()

println("-" ^ 60)
println("18.4 Jacobian Timing")
println("-" ^ 60)

# Warmup Jacobian
for _ in 1:5
	compute_analytical_jacobian_ale!(J, u0, nothing, 0.0, ode_fn, jac_cache)
end

jac_bench = @benchmark compute_analytical_jacobian_ale!($J, $u0, nothing, 0.0, $ode_fn, $jac_cache)
println("  Median time: $(median(jac_bench.times) / 1e6) ms")
println("  Min time:    $(minimum(jac_bench.times) / 1e6) ms")
println("  Allocations: $(jac_bench.allocs)")
println("  Memory:      $(jac_bench.memory) bytes")

# Ratio
jac_to_rhs_ratio = median(jac_bench.times) / median(rhs_bench.times)
println("\n  Jacobian/RHS ratio: $(round(jac_to_rhs_ratio, digits=1))×")
println("  (For reference: Finite differences would be ~$(n_cells * (1 + NC))× = $(n_cells * (1 + NC))×)")
println()

println("-" ^ 60)
println("18.5 Estimated Solve Time")
println("-" ^ 60)

# Estimate total solve time based on profiling
rhs_time_ms = median(rhs_bench.times) / 1e6
jac_time_ms = median(jac_bench.times) / 1e6

# Typical solve statistics (rough estimates for stiff problems)
est_steps = 1000
est_rhs_per_step = 10      # Newton iterations, error estimation
est_jac_per_step = 2       # Jacobian updates per step

total_rhs_calls = est_steps * est_rhs_per_step
total_jac_calls = est_steps * est_jac_per_step

est_rhs_time = total_rhs_calls * rhs_time_ms / 1000  # Convert to seconds
est_jac_time = total_jac_calls * jac_time_ms / 1000
est_total_time = est_rhs_time + est_jac_time

println("Estimated time for $(est_steps)-step solve:")
println("  RHS calls:      $(total_rhs_calls) × $(round(rhs_time_ms, digits=3)) ms = $(round(est_rhs_time, digits=1)) s")
println("  Jacobian calls: $(total_jac_calls) × $(round(jac_time_ms, digits=3)) ms = $(round(est_jac_time, digits=1)) s")
println("  Total:          ~$(round(est_total_time, digits=1)) s")
println()

# Speedup estimates
println("Potential speedups from optimizations:")
println("  - Reduce to 200 cells:     ~$(round((n_cells/200)^2, digits=1))× faster (RHS) + ~$(round((n_cells/200)^1.5, digits=1))× faster (Jacobian)")
println("  - Simplify to 5 components: ~$(round((NC/5)^1.5, digits=1))× faster")
println("  - Disable ALE mesh motion:  ~1.5× faster")
println("  - Use QNDF solver:          ~1.2-2× faster (problem dependent)")
println()

println("=" ^ 80)
println("OPTIMIZATION RECOMMENDATIONS")
println("=" ^ 80)
println()

# Analyze results and provide recommendations
rhs_pct = 100.0 * est_rhs_time / est_total_time
jac_pct = 100.0 * est_jac_time / est_total_time

println("Performance breakdown:")
println("  RHS computation:  $(round(rhs_pct, digits=1))% of solve time")
println("  Jacobian:         $(round(jac_pct, digits=1))% of solve time")
println()

if rhs_pct > 70
	println("⚠ RHS dominates computation time ($(round(rhs_pct, digits=1))%)")
	println("  Priority optimizations:")
	println("  1. Reduce mesh resolution (biggest impact)")
	println("  2. Simplify chemistry (reduce number of components)")
	println("  3. Profile RHS components to find specific bottlenecks")
elseif jac_pct > 70
	println("⚠ Jacobian dominates computation time ($(round(jac_pct, digits=1))%)")
	println("  Priority optimizations:")
	println("  1. Reduce mesh resolution (Jacobian scales as O(N²))")
	println("  2. Increase Jacobian reuse (larger dtmax, looser tolerances)")
	println("  3. Consider iterative linear solvers for very large problems")
else
	println("✓ Balanced RHS/Jacobian time - both should be optimized")
	println("  Priority optimizations:")
	println("  1. Reduce mesh resolution (affects both)")
	println("  2. Simplify chemistry if possible")
	println("  3. Optimize solver settings (tolerances, dtmax)")
end
println()

println("=" ^ 80)
println("NEXT STEPS")
println("=" ^ 80)
println()
println("1. Focus on functions with > 0 bytes allocation")
println("2. Focus on functions taking >5% of CPU time (see profile above)")
println("3. Use @code_warntype to identify type instabilities")
println("4. Check for:")
println("   - Array slicing (creates views that allocate)")
println("   - Tuple indexing with runtime indices (allocates)")
println("   - Closures capturing variables (can allocate)")
println("   - Type-unstable return types")
println("5. Pay special attention to:")
println("   - Property function calls (CallableProperty, StateDependentProperty)")
println("   - Component/reaction tuple access patterns")
println("   - Mixing rule calculations")
println("   - Cache array access patterns")
println("6. Experiment with solver configuration:")
println("   - Try different solvers (Rodas5P, QNDF, Rodas4P)")
println("   - Adjust tolerances (abstol, reltol)")
println("   - Limit timestep (dtmax) to reduce Jacobian updates")
println()
println("To save profile data for visualization:")
println("  using ProfileView")
println("  ProfileView.view()")
println()
println("To generate flamegraph:")
println("  using ProfileSVG")
println("  ProfileSVG.save(\"profile.svg\")")
println()
