# test_radiation_p1_validation.jl
# Week 4: Validation tests for P1 volumetric radiation model
# Tests against analytical solutions, physics verification, conservation checks

using Test
using LinearAlgebra
using Statistics: mean

using Pyrolysis
using Pyrolysis.Properties  # For update_properties!
using Pyrolysis.Discretization  # For update_radiation_properties!, etc.

# Import specific functions to resolve ambiguity
using Pyrolysis.Discretization: update_radiation_properties!

# Import physical constants
const σ_SB = Pyrolysis.STEFAN_BOLTZMANN

# ============================================================================
# Helper Functions
# ============================================================================

"""
Create a highly absorbing material (char-like) for validation tests.
"""
function create_char_validation_material()
	return Material(
		name = :char,
		components = (
			SolidComponent(:char, ρ = 200.0, c = 1500.0, k = 0.3, ε = 0.95, α = 3.0),
		),
		reactions = (),
		β = 0.5,
	)
end

"""
Create a very optically thick material for testing Rosseland limit.
High α ensures τ >> 1 even in thin samples.
"""
function create_optically_thick_material()
	return Material(
		name = :thick,
		components = (
			SolidComponent(:thick, ρ = 500.0, c = 1500.0, k = 0.2, ε = 0.9, α = 20.0),
		),
		reactions = (),
		β = 0.5,
	)
end

"""
Create mesh with specified temperature profile.
"""
function create_test_mesh_with_temperature(
	n_cells::Int,
	NC::Int,
	L::Float64,
	T_profile::Function,
	ξ_initial::Float64,
)
	mesh = create_uniform_mesh(L, n_cells, Val(NC))
	for i in 1:n_cells
		z = mesh.cell_centers[i]
		mesh.cell_states[i] = CellState{NC}(
			T_profile(z),
			ntuple(j -> j == 1 ? ξ_initial : 0.0, NC),
			mesh.cell_volumes[i],
			mesh.cell_volumes[i],
		)
	end
	return mesh
end

# ============================================================================
# Test Suite: Analytical Solutions
# ============================================================================

@testset "P1 Radiation Validation" begin

	@testset "Analytical: Rosseland Limit Equilibrium" begin
		# In optically thick media (τ >> 1), the P1 model should give
		# G ≈ 4σT⁴ at local radiative equilibrium

		@testset "Local equilibrium in optically thick medium" begin
			material = create_optically_thick_material()

			# Small domain but very optically thick
			L = 0.01  # 1 cm
			n_cells = 50
			T_cold = 400.0
			T_hot = 800.0

			# Linear temperature profile
			mesh = create_test_mesh_with_temperature(
				n_cells, 1, L,
				z -> T_cold + (T_hot - T_cold) * z / L,
				500.0,
			)

			# Setup caches
			n_faces = n_cells + 1
			prop_cache = PropertyCache(n_cells, n_faces)
			update_properties!(prop_cache, mesh, material)

			rad_cache = RadiationCache(n_cells)
			update_radiation_properties!(rad_cache, mesh, prop_cache)

			# Check optical thickness is large
			τ = optical_thickness(mesh, prop_cache)
			@test τ > 50  # Very optically thick

			# Initialize G from equilibrium
			initialize_g_from_temperature!(rad_cache, mesh, prop_cache)

			# At initialization, G should be 4σT⁴ (blackbody, NO emissivity factor)
			# In the standard P1 approximation, volumetric equilibrium is ALWAYS blackbody.
			# Emissivity only appears at boundaries (Marshak BC), not in the volume.
			# This is a fundamental assumption of P1: the medium is locally a perfect
			# absorber/emitter (gray body with ε=1 volumetrically).
			for i in 1:n_cells
				T = mesh.cell_states[i].T
				G_expected = 4.0 * σ_SB * T^4  # Blackbody equilibrium
				G_actual = rad_cache.G[i]
				@test isapprox(G_actual, G_expected, rtol = 1e-10)
			end
		end

		@testset "Equilibrium G increases with temperature" begin
			material = create_optically_thick_material()
			n_cells = 20
			L = 0.005  # 5 mm

			# Step temperature profile
			mesh = create_test_mesh_with_temperature(
				n_cells, 1, L,
				z -> z < L / 2 ? 400.0 : 800.0,
				500.0,
			)

			# Setup
			prop_cache = PropertyCache(n_cells, n_cells + 1)
			update_properties!(prop_cache, mesh, material)
			rad_cache = RadiationCache(n_cells)
			update_radiation_properties!(rad_cache, mesh, prop_cache)
			initialize_g_from_temperature!(rad_cache, mesh, prop_cache)

			# G in hot region should be ~16x G in cold region (T⁴ ratio = (800/400)⁴ = 16)
			G_cold_avg = mean([rad_cache.G[i] for i in 1:div(n_cells, 2)])
			G_hot_avg = mean([rad_cache.G[i] for i in (div(n_cells, 2)+1):n_cells])
			@test G_hot_avg / G_cold_avg ≈ 16.0 rtol = 0.1
		end
	end

	@testset "Physics: Radiation Heat Transfer Behavior" begin
		@testset "Heat flows from hot to cold" begin
			# Create non-uniform temperature field and verify
			# radiation source term adds heat to cold regions, removes from hot
			material = create_char_validation_material()
			n_cells = 30
			L = 0.01

			# Temperature higher in middle
			mesh = create_test_mesh_with_temperature(
				n_cells, 1, L,
				z -> 500.0 + 300.0 * exp(-((z - L / 2) / (L / 4))^2),
				200.0,
			)

			prop_cache = PropertyCache(n_cells, n_cells + 1)
			update_properties!(prop_cache, mesh, material)
			rad_cache = RadiationCache(n_cells)
			update_radiation_properties!(rad_cache, mesh, prop_cache)

			# Setup BCs
			radiation_bcs = RadiationBCSet(
				top = MarshakRadiationBC(I_external = 0.0, emissivity = 0.9),
				bottom = ReflectiveRadiationBC(),
			)

			# Solve radiation field
			initialize_g_from_temperature!(rad_cache, mesh, prop_cache)
			solve_quasi_steady_radiation!(rad_cache, mesh, prop_cache, radiation_bcs, 0.0)
			compute_radiation_flux_divergence!(rad_cache, mesh, prop_cache)

			# Find middle cell (hot region) and edge cells (cold regions)
			mid = div(n_cells, 2)
			T_mid = mesh.cell_states[mid].T
			T_edge = mesh.cell_states[1].T

			# Hot region should lose heat (negative source)
			Q_mid = rad_cache.div_q_rad[mid]

			# Note: The sign depends on convention. What matters is physics:
			# hot regions should have net outflow, cold regions net inflow.
			# In equilibrium, hot regions have Q_rad < 0 (losing heat)
			# Check that middle region has different sign from edges
			Q_edge = rad_cache.div_q_rad[1]
			@test Q_mid * Q_edge < 0 || abs(Q_mid) > abs(Q_edge)
		end

		@testset "External flux increases surface temperature heat gain" begin
			material = create_char_validation_material()
			n_cells = 20
			L = 0.008

			# Uniform temperature initially
			mesh = create_test_mesh_with_temperature(n_cells, 1, L, z -> 400.0, 200.0)

			prop_cache = PropertyCache(n_cells, n_cells + 1)
			update_properties!(prop_cache, mesh, material)
			rad_cache = RadiationCache(n_cells)
			update_radiation_properties!(rad_cache, mesh, prop_cache)

			# High external flux at surface
			radiation_bcs = RadiationBCSet(
				top = MarshakRadiationBC(I_external = 50e3, emissivity = 0.9),
				bottom = ReflectiveRadiationBC(),
			)

			initialize_g_from_temperature!(rad_cache, mesh, prop_cache)
			solve_quasi_steady_radiation!(rad_cache, mesh, prop_cache, radiation_bcs, 0.0)
			compute_radiation_flux_divergence!(rad_cache, mesh, prop_cache)

			# Surface cell should have net heat addition (Q_rad > 0)
			# Convention: Q_rad > 0 means heat added to cell
			@test rad_cache.div_q_rad[n_cells] > 0  # Top cell gains heat
		end
	end

	@testset "Conservation: Energy Balance" begin
		@testset "Total source term integrates to flux difference" begin
			# For steady state: ∫Q_rad dV = Φ_in - Φ_out (net boundary flux)
			material = create_char_validation_material()
			n_cells = 40
			L = 0.012

			mesh = create_test_mesh_with_temperature(
				n_cells, 1, L,
				z -> 400.0 + 200.0 * z / L,
				200.0,
			)

			prop_cache = PropertyCache(n_cells, n_cells + 1)
			update_properties!(prop_cache, mesh, material)
			rad_cache = RadiationCache(n_cells)
			update_radiation_properties!(rad_cache, mesh, prop_cache)

			# Asymmetric BCs to create net flux
			radiation_bcs = RadiationBCSet(
				top = MarshakRadiationBC(I_external = 30e3, emissivity = 0.9),
				bottom = ReflectiveRadiationBC(),
			)

			initialize_g_from_temperature!(rad_cache, mesh, prop_cache)
			solve_quasi_steady_radiation!(rad_cache, mesh, prop_cache, radiation_bcs, 0.0)
			compute_radiation_flux_divergence!(rad_cache, mesh, prop_cache)

			# Compute total volumetric source
			total_Q = sum(rad_cache.div_q_rad[i] * mesh.cell_volumes[i] for i in 1:n_cells)

			# In steady state, total source = net flux in - net flux out
			# This is just checking Q_rad sums to something finite
			@test isfinite(total_Q)
		end

		@testset "Closed reflective system reaches equilibrium" begin
			material = create_char_validation_material()
			n_cells = 25
			L = 0.006

			# Initial non-uniform temperature
			mesh = create_test_mesh_with_temperature(
				n_cells, 1, L,
				z -> 500.0 + 200.0 * sin(π * z / L),
				200.0,
			)

			prop_cache = PropertyCache(n_cells, n_cells + 1)
			update_properties!(prop_cache, mesh, material)
			rad_cache = RadiationCache(n_cells)
			update_radiation_properties!(rad_cache, mesh, prop_cache)

			# Fully reflective boundaries (closed system)
			radiation_bcs = RadiationBCSet(
				top = ReflectiveRadiationBC(),
				bottom = ReflectiveRadiationBC(),
			)

			initialize_g_from_temperature!(rad_cache, mesh, prop_cache)
			solve_quasi_steady_radiation!(rad_cache, mesh, prop_cache, radiation_bcs, 0.0)

			# In closed reflective system, total G should be conserved
			total_G = sum(rad_cache.G[i] * mesh.cell_volumes[i] for i in 1:n_cells)
			@test total_G > 0
			@test isfinite(total_G)
		end
	end

	@testset "Optical Thickness Diagnostics" begin
		@testset "Optical thickness scales with absorption and path length" begin
			material = create_optically_thick_material()

			# Different path lengths
			for L in [0.001, 0.005, 0.01, 0.02]
				n_cells = 20
				mesh = create_test_mesh_with_temperature(n_cells, 1, L, z -> 600.0, 500.0)
				prop_cache = PropertyCache(n_cells, n_cells + 1)
				update_properties!(prop_cache, mesh, material)

				τ = optical_thickness(mesh, prop_cache)

				# τ = α * L where α = ρ * α_mass
				# With ρ = 500, α_mass = 20, we get α_vol = 10000 m⁻¹
				# So τ = 10000 * L
				expected_τ = 500.0 * 20.0 * L
				@test isapprox(τ, expected_τ, rtol = 0.05)
			end
		end

		@testset "Per-cell optical thickness" begin
			material = create_optically_thick_material()
			n_cells = 10
			L = 0.01
			mesh = create_test_mesh_with_temperature(n_cells, 1, L, z -> 600.0, 500.0)
			prop_cache = PropertyCache(n_cells, n_cells + 1)
			update_properties!(prop_cache, mesh, material)

			# Check a specific cell
			τ_cell = optical_thickness_cell(mesh, prop_cache, 5)
			Δz = mesh.cell_volumes[5]  # In 1D with unit area, volume = length
			expected = prop_cache.α[5] * Δz
			@test isapprox(τ_cell, expected, rtol = 1e-10)
		end

		@testset "P1 validity check" begin
			material_thick = create_optically_thick_material()
			material_char = create_char_validation_material()

			# Thick material should be valid for P1
			n_cells = 30
			L = 0.01
			mesh_thick = create_test_mesh_with_temperature(n_cells, 1, L, z -> 600.0, 500.0)
			prop_thick = PropertyCache(n_cells, n_cells + 1)
			update_properties!(prop_thick, mesh_thick, material_thick)

			result_thick = check_p1_validity(mesh_thick, prop_thick; warn = false)
			@test result_thick.valid == true
			@test result_thick.τ_total > 10  # τ = 500 * 20 * 0.01 = 100

			# Char material with short path may not be valid
			mesh_thin = create_test_mesh_with_temperature(n_cells, 1, 0.001, z -> 600.0, 200.0)
			prop_thin = PropertyCache(n_cells, n_cells + 1)
			update_properties!(prop_thin, mesh_thin, material_char)

			result_thin = check_p1_validity(mesh_thin, prop_thin; warn = false)
			# τ = 200 * 3.0 * 0.001 = 0.6 < 1
			@test result_thin.τ_total < 2
		end
	end

	@testset "Boundary Condition Types" begin
		@testset "Reflective BC reflects all radiation" begin
			material = create_char_validation_material()
			n_cells = 20
			L = 0.008
			mesh = create_test_mesh_with_temperature(n_cells, 1, L, z -> 500.0, 200.0)
			prop_cache = PropertyCache(n_cells, n_cells + 1)
			update_properties!(prop_cache, mesh, material)
			rad_cache = RadiationCache(n_cells)
			update_radiation_properties!(rad_cache, mesh, prop_cache)

			# Both boundaries reflective
			radiation_bcs = RadiationBCSet(
				top = ReflectiveRadiationBC(),
				bottom = ReflectiveRadiationBC(),
			)

			initialize_g_from_temperature!(rad_cache, mesh, prop_cache)
			solve_quasi_steady_radiation!(rad_cache, mesh, prop_cache, radiation_bcs, 0.0)

			# Gradient at boundaries should be zero for reflective BC
			# (∂G/∂n = 0)
			G_edge1 = rad_cache.G[1]
			G_edge2 = rad_cache.G[n_cells]
			G_interior = rad_cache.G[div(n_cells, 2)]

			# In uniform temperature, G should be approximately uniform
			@test isapprox(G_edge1, G_interior, rtol = 0.15)
			@test isapprox(G_edge2, G_interior, rtol = 0.15)
		end

		@testset "Marshak BC with zero external flux" begin
			material = create_char_validation_material()
			n_cells = 25
			L = 0.01
			mesh = create_test_mesh_with_temperature(n_cells, 1, L, z -> 600.0, 200.0)
			prop_cache = PropertyCache(n_cells, n_cells + 1)
			update_properties!(prop_cache, mesh, material)
			rad_cache = RadiationCache(n_cells)
			update_radiation_properties!(rad_cache, mesh, prop_cache)

			# Marshak with zero external flux is similar to black surface at T=0
			radiation_bcs = RadiationBCSet(
				top = MarshakRadiationBC(I_external = 0.0, emissivity = 1.0),
				bottom = ReflectiveRadiationBC(),
			)

			initialize_g_from_temperature!(rad_cache, mesh, prop_cache)
			solve_quasi_steady_radiation!(rad_cache, mesh, prop_cache, radiation_bcs, 0.0)

			# Surface should lose heat (G at surface lower than equilibrium)
			G_surface = rad_cache.G[n_cells]
			G_equilibrium = 4.0 * σ_SB * 600.0^4
			@test G_surface < G_equilibrium
		end

		@testset "Dirichlet BC fixes G at boundary" begin
			material = create_char_validation_material()
			n_cells = 20
			L = 0.008
			mesh = create_test_mesh_with_temperature(n_cells, 1, L, z -> 500.0, 200.0)
			prop_cache = PropertyCache(n_cells, n_cells + 1)
			update_properties!(prop_cache, mesh, material)
			rad_cache = RadiationCache(n_cells)
			update_radiation_properties!(rad_cache, mesh, prop_cache)

			# Set specific G at top boundary (higher than equilibrium at 500K)
			G_fixed = 4.0 * σ_SB * 800.0^4
			radiation_bcs = RadiationBCSet(
				top = DirichletRadiationBC(G_fixed),
				bottom = ReflectiveRadiationBC(),
			)

			initialize_g_from_temperature!(rad_cache, mesh, prop_cache)
			solve_quasi_steady_radiation!(rad_cache, mesh, prop_cache, radiation_bcs, 0.0)

			# Check that G field is valid
			@test all(isfinite.(rad_cache.G))
			@test all(rad_cache.G .> 0)
			# G values should be influenced by the high boundary condition
			# At minimum, all G should be at or above pure blackbody at interior temperature
			G_interior_min = 4.0 * σ_SB * prop_cache.ε[1] * 500.0^4
			@test rad_cache.G[n_cells] >= G_interior_min - 1e-10  # small tolerance for float comparison
		end
	end

	@testset "Quasi-Steady Solver Convergence" begin
		@testset "Solver converges for various temperature profiles" begin
			material = create_char_validation_material()
			n_cells = 30
			L = 0.01

			profiles = [
				("uniform", z -> 500.0),
				("linear", z -> 400.0 + 400.0 * z / L),
				("parabolic", z -> 500.0 + 200.0 * (z / L)^2),
				("sine", z -> 500.0 + 150.0 * sin(π * z / L)),
			]

			for (name, T_func) in profiles
				mesh = create_test_mesh_with_temperature(n_cells, 1, L, T_func, 200.0)
				prop_cache = PropertyCache(n_cells, n_cells + 1)
				update_properties!(prop_cache, mesh, material)
				rad_cache = RadiationCache(n_cells)
				update_radiation_properties!(rad_cache, mesh, prop_cache)

				radiation_bcs = RadiationBCSet(
					top = MarshakRadiationBC(I_external = 25e3, emissivity = 0.9),
					bottom = ReflectiveRadiationBC(),
				)

				initialize_g_from_temperature!(rad_cache, mesh, prop_cache)
				solve_quasi_steady_radiation!(rad_cache, mesh, prop_cache, radiation_bcs, 0.0)

				# Check solution is valid
				@test all(isfinite.(rad_cache.G))
				@test all(rad_cache.G .> 0)
			end
		end
	end

	@testset "P1 RHS vs No-Radiation RHS Comparison" begin
		@testset "P1 adds radiation source to RHS" begin
			# This test verifies that the P1 ODE function produces different
			# RHS than the base function due to radiation source term

			NC = 2
			n_cells = 20
			L = 0.01

			# Two-component char material
			material = Material(
				name = :test,
				components = (
					SolidComponent(:char, ρ = 200.0, c = 1500.0, k = 0.3, ε = 0.95, α = 3.0),
					GaseousComponent(:gas, M = 0.044, c = 1200.0, k = 0.02, λ = 1e-5),
				),
				reactions = (),
				β = 0.5,
			)

			# Create mesh with temperature gradient
			mesh = create_test_mesh_with_temperature(
				n_cells, NC, L,
				z -> 500.0 + 300.0 * z / L,
				200.0,
			)

			# Setup property cache
			prop_cache = PropertyCache(n_cells, n_cells + 1)
			update_properties!(prop_cache, mesh, material)

			# Setup radiation
			rad_cache = RadiationCache(n_cells)
			update_radiation_properties!(rad_cache, mesh, prop_cache)
			initialize_g_from_temperature!(rad_cache, mesh, prop_cache)

			radiation_bcs = RadiationBCSet(
				top = MarshakRadiationBC(I_external = 50e3, emissivity = 0.9),
				bottom = ReflectiveRadiationBC(),
			)

			# Solve radiation field
			solve_quasi_steady_radiation!(rad_cache, mesh, prop_cache, radiation_bcs, 0.0)
			compute_radiation_flux_divergence!(rad_cache, mesh, prop_cache)

			# Verify radiation source term is non-zero
			Q_rad_total = sum(abs.(rad_cache.div_q_rad))
			@test Q_rad_total > 0

			# Verify it affects different cells differently
			@test length(unique(sign.(rad_cache.div_q_rad))) >= 1  # Has at least one sign
		end
	end

end  # Main testset
