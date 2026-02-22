"""
Tests for Discretization Module

Tests cover:
- Divergence operator
- Gradient reconstruction  
- Flux limiters
- RHS assembly
- MMS verification (Method of Manufactured Solutions)
- Conservation checks
"""

using Test
using Pyrolysis
using Pyrolysis.Core
using Pyrolysis.Geometry
using Pyrolysis.Properties
using Pyrolysis.Physics
using Pyrolysis.Discretization

@testset "Discretization Module" begin

	@testset "Flux Limiters" begin
		@testset "Limiter Functions" begin
			# Minmod limiter
			@test minmod(0.5) ≈ 0.5
			@test minmod(1.5) ≈ 1.0
			@test minmod(-0.5) ≈ 0.0
			@test minmod(0.0) ≈ 0.0
			@test @inferred minmod(0.5) isa Float64

			# Van Leer limiter
			@test vanleer(1.0) ≈ 1.0
			@test vanleer(0.0) ≈ 0.0
			@test vanleer(-1.0) ≈ 0.0
			@test @inferred vanleer(0.5) isa Float64

			# Superbee limiter
			@test superbee(1.0) ≈ 1.0
			@test superbee(0.5) ≈ 1.0
			@test superbee(2.0) ≈ 2.0
			@test superbee(0.0) ≈ 0.0
			@test @inferred superbee(0.5) isa Float64
		end

		@testset "Limiter Enum Dispatch" begin
			@test apply_limiter(0.5, MINMOD) ≈ minmod(0.5)
			@test apply_limiter(0.5, VANLEER) ≈ vanleer(0.5)
			@test apply_limiter(0.5, SUPERBEE) ≈ superbee(0.5)
			@test apply_limiter(0.5, NONE) ≈ 0.0
			@test apply_limiter(0.5, CENTRAL) ≈ 1.0
		end

		@testset "Limited Face Value" begin
			# Uniform field: should return φ_L
			@test limited_face_value(1.0, 1.0, 1.0, MINMOD) ≈ 1.0

			# Linear field: φ = [0, 1, 2], r = 1
			# limited_face_value(0, 1, 2, MINMOD) = 1 + 0.5*minmod(1)*(2-1) = 1.5
			@test limited_face_value(0.0, 1.0, 2.0, MINMOD) ≈ 1.5

			# Non-monotone: φ = [2, 1, 2], r = -1 (negative slope ratio)
			# minmod(-1) = 0, so result should be 1.0 (first order)
			@test limited_face_value(2.0, 1.0, 2.0, MINMOD) ≈ 1.0
		end
	end

	@testset "Gradient Reconstruction" begin
		# Create a simple mesh
		mesh = create_uniform_mesh(0.01, 10, Val(2))
		n = n_cells(mesh)

		# Set up linear temperature profile: T = 300 + 20000*z
		# dT/dz = 20000 K/m = 200 K/cm
		for i in 1:n
			z = cell_center(mesh, i)
			T = 300.0 + 20000.0 * z
			mesh.cell_states[i].T = T
		end

		# Compute gradients
		∇T = zeros(n)
		compute_temperature_gradients!(∇T, mesh)

		# Interior cells should have gradient ≈ 20000
		for i in 2:(n-1)
			@test ∇T[i] ≈ 20000.0 rtol=0.01
		end

		@testset "Face Value Interpolation" begin
			# Central interpolation
			@test face_value_central(1.0, 3.0) ≈ 2.0
			@test @inferred face_value_central(1.0, 3.0) isa Float64

			# Upwind interpolation
			@test face_value_upwind(1.0, 3.0, 1.0) ≈ 1.0   # v > 0: use left
			@test face_value_upwind(1.0, 3.0, -1.0) ≈ 3.0  # v < 0: use right
			@test face_value_upwind(1.0, 3.0, 0.0) ≈ 1.0   # v = 0: use left
			@test @inferred face_value_upwind(1.0, 3.0, 1.0) isa Float64

			# Weighted interpolation
			@test face_value_weighted(1.0, 3.0, 0.5) ≈ 2.0
			@test face_value_weighted(1.0, 3.0, 1.0) ≈ 1.0
			@test face_value_weighted(1.0, 3.0, 0.0) ≈ 3.0
		end
	end

	@testset "Divergence Operator" begin
		# Create mesh and test divergence
		mesh = create_uniform_mesh(0.01, 10, Val(2))
		n_cells_val = n_cells(mesh)
		n_faces_val = n_faces(mesh)

		@testset "Constant Flux" begin
			# Constant flux field: divergence should be zero
			fluxes = fill(1.0, n_faces_val)
			div = zeros(n_cells_val)

			divergence_1d!(div, fluxes, mesh)

			for i in 1:n_cells_val
				@test div[i] ≈ 0.0 atol=1e-12
			end
		end

		@testset "Linear Flux" begin
			# Linear flux field: F = z
			# dF/dz = 1.0
			fluxes = copy(mesh.face_positions)

			div = zeros(n_cells_val)
			divergence_1d!(div, fluxes, mesh)

			# Each cell should have div ≈ 1.0
			for i in 1:n_cells_val
				@test div[i] ≈ 1.0 rtol=0.01
			end
		end

		@testset "Conservation Check" begin
			# For any flux field, sum of (divergence * volume) = boundary flux difference
			# ∫ ∇·F dV = ∮ F·n dA
			fluxes = [sin(2π * f / n_faces_val) for f in 1:n_faces_val]
			div = zeros(n_cells_val)
			divergence_1d!(div, fluxes, mesh)

			# Integral of divergence
			int_div = sum(div[i] * cell_volume(mesh, i) for i in 1:n_cells_val)

			# Boundary flux difference (right boundary - left boundary) * area
			A = mesh.cross_section_area
			boundary_flux_diff = (fluxes[end] - fluxes[1]) * A

			@test int_div ≈ boundary_flux_diff rtol=1e-10
		end
	end

	@testset "RHS Cache" begin
		# Create a simple material
		virgin = SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86)
		gas = GaseousComponent(:MMA, M = 0.100, c = 1500, k = 0.02, λ = 1e-5)
		rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
		mat = Material(name = :PMMA, components = (virgin, gas), reactions = (rxn,))

		# Create mesh
		mesh = create_uniform_mesh(0.01, 10, Val(2))

		# Create caches
		rhs_cache = RHSCache(mesh, mat)

		@test length(rhs_cache.q_cond) == n_faces(mesh)
		@test length(rhs_cache.Q_rxn) == n_cells(mesh)
		@test size(rhs_cache.J_gas) == (n_components(mat), n_faces(mesh))
		@test size(rhs_cache.dξ_rxn) == (n_components(mat), n_cells(mesh))
	end

	@testset "Steady State Heat Conduction" begin
		# Test case: steady state conduction with fixed temperatures
		# T(0) = T_hot, T(L) = T_cold
		# Exact solution: T(z) = T_hot + (T_cold - T_hot) * z / L

		L = 0.01  # 1 cm
		n = 20
		T_hot = 600.0
		T_cold = 300.0

		# Create mesh with simple material (no reactions)
		mesh = create_uniform_mesh(L, n, Val(1))

		# Single inert component
		inert = SolidComponent(:inert, ρ = 1000, c = 1000, k = 0.2, ε = 0.9)
		mat = Material(name = :test, components = (inert,))

		# Set initial linear profile (exact steady state solution)
		for i in 1:n
			z = mesh.cell_centers[i]
			T_exact = T_hot + (T_cold - T_hot) * z / L
			mesh.cell_states[i].T = T_exact
			mesh.cell_states[i].ξ = (1000.0,)
		end

		# Create property cache and update
		prop_cache = PropertyCache(mesh)
		update_properties!(prop_cache, mesh, mat)

		# Compute conductive fluxes
		q_cond = zeros(n_faces(mesh))
		compute_conductive_fluxes!(q_cond, mesh, prop_cache)

		# For linear T profile, INTERIOR fluxes should be constant
		# q = -k * dT/dz = -0.2 * (T_cold - T_hot) / L = -0.2 * (-300) / 0.01 = 6000 W/m²
		expected_flux = -0.2 * (T_cold - T_hot) / L

		# Check interior faces only (faces 2 to n_faces-1)
		for f in 2:(n_faces(mesh)-1)
			@test q_cond[f] ≈ expected_flux rtol=0.01
		end

		# Divergence of constant flux should be zero for interior cells
		div_q = zeros(n)
		divergence_1d!(div_q, q_cond, mesh)

		# Interior cells (not affected by boundary zeros) should have near-zero div
		for i in 2:(n-1)
			@test div_q[i] ≈ 0.0 atol=1e-6
		end
	end

	@testset "MMS Verification - Conduction" begin
		# Method of Manufactured Solutions
		# For a sinusoidal temperature profile, check that flux computation
		# gives expected gradient

		L = 0.01
		T0 = 300.0
		A = 100.0  # Amplitude
		k = 0.2

		# Create mesh
		n = 40
		mesh = create_uniform_mesh(L, n, Val(1))

		# Set manufactured solution: T(z) = T0 + A * sin(π * z / L)
		for i in 1:n
			z = mesh.cell_centers[i]
			T_exact = T0 + A * sin(π * z / L)
			mesh.cell_states[i].T = T_exact
			mesh.cell_states[i].ξ = (1000.0,)
		end

		# Create material and caches
		inert = SolidComponent(:inert, ρ = 1000, c = 1000, k = k, ε = 0.9)
		mat = Material(name = :test, components = (inert,))
		prop_cache = PropertyCache(mesh)
		update_properties!(prop_cache, mesh, mat)

		# Compute conductive fluxes
		q_cond = zeros(n_faces(mesh))
		compute_conductive_fluxes!(q_cond, mesh, prop_cache)

		# Check interior fluxes match expected gradient
		# q = -k * dT/dz = -k * A * (π/L) * cos(π * z / L)
		max_flux = k * A * (π / L)  # Maximum expected flux magnitude
		for f in 2:(n_faces(mesh)-1)
			z_f = mesh.face_positions[f]
			expected_q = -k * A * (π / L) * cos(π * z_f / L)
			# Use both relative and absolute tolerance to handle near-zero values
			@test isapprox(q_cond[f], expected_q, rtol = 0.05, atol = 1e-10 * max_flux)
		end
	end

	@testset "ALE Flux Functions" begin
		@testset "ALE Convective Flux" begin
			φ_L, φ_R = 1.0, 2.0

			# Case 1: v > w (material moving faster than mesh to the right)
			# Relative velocity positive → upwind is left
			@test ale_convective_flux(φ_L, φ_R, 2.0, 1.0) ≈ 1.0 * (2.0 - 1.0)

			# Case 2: v < w (mesh moving faster to the right)
			# Relative velocity negative → upwind is right
			@test ale_convective_flux(φ_L, φ_R, 1.0, 2.0) ≈ 2.0 * (1.0 - 2.0)

			# Case 3: v = w (Lagrangian)
			# No relative flux
			@test ale_convective_flux(φ_L, φ_R, 1.0, 1.0) ≈ 0.0

			@test @inferred ale_convective_flux(φ_L, φ_R, 1.0, 0.5) isa Float64
		end
	end

	@testset "Darcy-Fick RHS Functions" begin
		@testset "update_darcy_properties!" begin
			# Create a mesh with 10 cells for a charring material
			L = 0.01  # 1 cm
			n = 10
			mesh = create_uniform_mesh(L, n, Val(2))

			# Set up a temperature profile and component distribution
			T_hot = 600.0
			T_cold = 300.0
			for i in 1:n
				z = mesh.cell_centers[i]
				T = T_cold + (T_hot - T_cold) * (1 - z / L)  # Hot at bottom, cold at top
				# Solid component dominant, small amount of gas
				vol = mesh.cell_states[i].volume
				mesh.cell_states[i] = CellState{2}(T, (950.0, 50.0), vol)  # solid: 950, gas: 50 kg/m³
			end

			# Create material with solid and gas component
			solid = SolidComponent(:char, ρ = 1000.0, c = 1000.0, k = 0.2, ε = 0.9, κ = 1e-12)
			gas = GaseousComponent(:volatile, M = 0.030, c = 1500.0, k = 0.03, λ = 1e-5)
			mat = Material(name = :charring, components = (solid, gas))

			# Create property cache and update basic properties
			prop_cache = PropertyCache(mesh)
			update_properties!(prop_cache, mesh, mat)


			# Run update_darcy_properties! using type-stable API
			bc_set = BoundaryConditionSet()
			top_face_ids = mesh.boundary_faces[1]
			bottom_face_ids = mesh.boundary_faces[2]
			update_darcy_properties!(prop_cache, mesh, mat, bc_set, top_face_ids, bottom_face_ids)

			# Test 1: Pressure should be positive everywhere
			@test all(prop_cache.P .> 0)

			# Test 2: Porosity should be between 0 and 1
			@test all(0.0 .<= prop_cache.φ .<= 1.0)

			# Test 3: Permeability should be non-negative
			@test all(prop_cache.κ .>= 0.0)

			# Test 4: Viscosity should be positive
			@test all(prop_cache.μ .> 0.0)

			# Test 5: Pressure gradient should be defined at all faces
			@test length(prop_cache.∇P) == n_faces(mesh)

			# Test 6: Darcy velocity should be defined at all faces
			@test length(prop_cache.v_gas) == n_faces(mesh)

			# Test 7: Interior faces with temperature gradient should have pressure gradient
			# (higher T → higher P for ideal gas)
			# Since T decreases with z, P should decrease with z
			# So ∇P should be negative for interior faces (P_right < P_left in z direction)
			for f in 2:(n_faces(mesh)-1)
				# Interior faces should have non-zero gradient with temp gradient
				@test prop_cache.∇P[f] != 0.0
			end
		end

		@testset "compute_rhs_darcy_fick!" begin
			# Create a mesh for testing full RHS
			L = 0.01  # 1 cm
			n = 20
			mesh = create_uniform_mesh(L, n, Val(2))

			# Set up uniform initial state (should give zero RHS)
			T_uniform = 400.0
			for i in 1:n
				vol = mesh.cell_states[i].volume
				mesh.cell_states[i] = CellState{2}(T_uniform, (900.0, 100.0), vol)
			end

			# Create material
			solid = SolidComponent(:char, ρ = 1000.0, c = 1000.0, k = 0.2, ε = 0.9, κ = 1e-12)
			gas = GaseousComponent(:volatile, M = 0.030, c = 1500.0, k = 0.03, λ = 1e-5)
			mat = Material(name = :test, components = (solid, gas))

			# Create caches
			prop_cache = PropertyCache(mesh)
			update_properties!(prop_cache, mesh, mat)
			rhs_cache = RHSCache(n, n_faces(mesh), 2, 0)

			# Boundary conditions
			bc_set = BoundaryConditionSet(
				thermal_top = HeatFluxBC(0.0),
				thermal_bottom = AdiabaticBC(),
				mass_top = ImpermeableBC(),
				mass_bottom = ImpermeableBC(),
			)

			# Pre-cache face IDs
			top_face_ids = mesh.boundary_faces[1]
			bottom_face_ids = mesh.boundary_faces[2]

			# Output arrays
			dT = zeros(n)
			dξ = zeros(2, n)

			# Run compute_rhs_darcy_fick!
			compute_rhs_darcy_fick!(dT, dξ, mesh, prop_cache, rhs_cache,
				mat, bc_set, 0.0, top_face_ids, bottom_face_ids)

			# Test 1: Output arrays should have correct size
			@test length(dT) == n
			@test size(dξ) == (2, n)

			# Test 2: For uniform state with no reactions and no heat flux,
			# RHS should be small (only numerical artifacts)
			@test all(abs.(dT) .< 1e-6)

			# Test 3: For uniform concentrations with impermeable BCs,
			# species RHS should be small
			# Solid component (j=1) should have zero RHS always (no reactions)
			@test all(abs.(dξ[1, :]) .< 1e-10)
		end

		@testset "Darcy-Fick vs Diffusion Comparison" begin
			# Test that Darcy-Fick reduces to diffusion when velocity is small
			L = 0.01
			n = 20
			mesh = create_uniform_mesh(L, n, Val(2))

			# Uniform temperature but with concentration gradient for gas
			T = 400.0
			for i in 1:n
				z = mesh.cell_centers[i]
				# Linear gas concentration profile
				ξ_gas = 50.0 + 50.0 * z / L
				vol = mesh.cell_states[i].volume
				mesh.cell_states[i] = CellState{2}(T, (900.0, ξ_gas), vol)
			end

			# Create material (κ=0 disables Darcy, compare to diffusion)
			solid = SolidComponent(:char, ρ = 1000.0, c = 1000.0, k = 0.2, ε = 0.9, κ = 0.0)
			gas = GaseousComponent(:volatile, M = 0.030, c = 1500.0, k = 0.03, λ = 1e-5)
			mat = Material(name = :test, components = (solid, gas))

			# Create caches
			prop_cache = PropertyCache(mesh)
			update_properties!(prop_cache, mesh, mat)
			rhs_cache = RHSCache(n, n_faces(mesh), 2, 0)

			# BCs
			bc_set = BoundaryConditionSet(
				thermal_top = HeatFluxBC(0.0),
				thermal_bottom = AdiabaticBC(),
				mass_top = ImpermeableBC(),
				mass_bottom = ImpermeableBC(),
			)

			# Pre-cache face IDs
			top_face_ids = mesh.boundary_faces[1]
			bottom_face_ids = mesh.boundary_faces[2]

			dT_darcy = zeros(n)
			dξ_darcy = zeros(2, n)
			compute_rhs_darcy_fick!(dT_darcy, dξ_darcy, mesh, prop_cache, rhs_cache,
				mat, bc_set, 0.0, top_face_ids, bottom_face_ids)

			# Compute pure diffusion RHS
			dT_diff = zeros(n)
			dξ_diff = zeros(2, n)
			compute_rhs!(dT_diff, dξ_diff, mesh, prop_cache, rhs_cache,
				mat, bc_set, 0.0, top_face_ids, bottom_face_ids)

			# Results should be similar when permeability is very low
			# Gas species RHS should be comparable
			for i in 1:n
				# Allow small difference due to different code paths
				@test isapprox(dξ_darcy[2, i], dξ_diff[2, i], rtol = 0.1, atol = 1e-10)
			end
		end

		@testset "Pressure-Driven Flow Direction" begin
			# Test that Darcy flow goes from high to low pressure
			L = 0.01
			n = 20
			mesh = create_uniform_mesh(L, n, Val(2))

			# Hot gas at top (high pressure) should flow toward bottom (low pressure)
			# In our convention: z=0 is bottom, z=L is top
			for i in 1:n
				z = mesh.cell_centers[i]
				# Temperature increases with z (hot at top)
				T = 400.0 + 400.0 * z / L
				# Gas concentration uniform
				vol = mesh.cell_states[i].volume
				mesh.cell_states[i] = CellState{2}(T, (800.0, 200.0), vol)
			end

			# Create material
			solid = SolidComponent(:char, ρ = 1000.0, c = 1000.0, k = 0.2, ε = 0.9, κ = 1e-11)
			gas = GaseousComponent(:volatile, M = 0.030, c = 1500.0, k = 0.03, λ = 1e-5)
			mat = Material(name = :test, components = (solid, gas))

			# Create caches
			prop_cache = PropertyCache(mesh)
			update_properties!(prop_cache, mesh, mat)

			bc_set = BoundaryConditionSet()
			top_face_ids = mesh.boundary_faces[1]
			bottom_face_ids = mesh.boundary_faces[2]
			update_darcy_properties!(prop_cache, mesh, mat, bc_set, top_face_ids, bottom_face_ids)

			# Pressure should increase with z (due to higher temperature)
			# P = ξ_gas * R * T / M → higher T means higher P for same ξ_gas
			@test prop_cache.P[end] > prop_cache.P[1]

			# Darcy velocity should be negative (flow from high z to low z)
			# v_g = -(κ/μ) * ∇P
			# Since ∇P > 0 (P increases with z), v_g < 0
			interior_faces = 2:(n_faces(mesh)-1)
			avg_velocity = sum(prop_cache.v_gas[interior_faces]) / length(interior_faces)
			@test avg_velocity < 0.0
		end

		@testset "DarcyVelocityJacobianCoeffs" begin
			# Test Darcy velocity Jacobian coefficients
			# Use higher gas fraction to ensure non-zero porosity and permeability
			n = 10
			thickness = 0.006
			mesh = create_uniform_mesh(thickness, n, Val(2);
				T_initial = 400.0,
				ξ_initial = (500.0, 100.0),  # Moderate gas fraction for non-zero porosity
			)

			mat = Material(
				name = :test,
				components = (
					SolidComponent(:solid, ρ = ConstantProperty(1000.0),
						c = ConstantProperty(1500.0), k = ConstantProperty(0.2), ε = 0.9,
						κ = 1e-10),
					GaseousComponent(:gas, M = 0.044, c = ConstantProperty(1000.0),
						k = ConstantProperty(0.03), λ = ConstantProperty(1e-5)),
				),
				reactions = (),
			)

			prop_cache = PropertyCache(mesh)
			update_properties!(prop_cache, mesh, mat)

			bc_set = BoundaryConditionSet()
			top_face_ids = mesh.boundary_faces[1]
			bottom_face_ids = mesh.boundary_faces[2]
			update_darcy_properties!(prop_cache, mesh, mat, bc_set, top_face_ids, bottom_face_ids)

			# Verify we have non-zero porosity and permeability
			@test all(prop_cache.φ .> 0.0)
			@test all(prop_cache.κ .> 0.0)

			# Create Jacobian coefficient storage
			n_f = n_faces(mesh)
			jac_coeffs = DarcyVelocityJacobianCoeffs(n_f)

			# Compute Jacobian coefficients
			compute_darcy_velocity_jacobian!(jac_coeffs, mesh, prop_cache.κ)

			# Test 1: Boundary faces should have zero coefficients
			@test jac_coeffs.coeff_L[1] == 0.0
			@test jac_coeffs.coeff_R[1] == 0.0
			@test jac_coeffs.coeff_L[n_f] == 0.0
			@test jac_coeffs.coeff_R[n_f] == 0.0

			# Test 2: Interior faces should have non-zero coefficients
			for f in 2:(n_f-1)
				@test jac_coeffs.coeff_L[f] > 0.0  # ∂v/∂P_L > 0 (increasing P_L increases v)
				@test jac_coeffs.coeff_R[f] < 0.0  # ∂v/∂P_R < 0 (increasing P_R decreases v)
			end

			# Test 3: coeff_L = -coeff_R (equal magnitude, opposite sign)
			for f in 2:(n_f-1)
				@test isapprox(jac_coeffs.coeff_L[f], -jac_coeffs.coeff_R[f], rtol = 1e-10)
			end

			# Test 4: Cell indices should be correct
			for f in 2:(n_f-1)
				face = mesh.faces[f]
				@test jac_coeffs.left_ids[f] == face.cell_ids[1]
				@test jac_coeffs.right_ids[f] == face.cell_ids[2]
			end

			# Test 5: Verify against finite difference
			# ∂v/∂P_L ≈ (v(P_L + ε) - v(P_L)) / ε
			ε = 1.0  # 1 Pa perturbation
			f_test = n_f ÷ 2  # Middle face
			left_id = jac_coeffs.left_ids[f_test]
			right_id = jac_coeffs.right_ids[f_test]

			# Get current velocity
			v_base = prop_cache.v_gas[f_test]

			# Perturb left pressure
			P_orig = prop_cache.P[left_id]
			prop_cache.P[left_id] = P_orig + ε

			# Recompute velocity at this face
			d_LR = abs(mesh.cell_centers[right_id] - mesh.cell_centers[left_id])
			v_pert = darcy_velocity_face(
				prop_cache.κ[left_id], prop_cache.κ[right_id],
				mesh.cell_states[left_id].T, mesh.cell_states[right_id].T,
				prop_cache.P[left_id], prop_cache.P[right_id], d_LR,
			)

			# Restore
			prop_cache.P[left_id] = P_orig

			# Check FD matches analytical coefficient
			dv_dP_fd = (v_pert - v_base) / ε
			@test isapprox(jac_coeffs.coeff_L[f_test], dv_dP_fd, rtol = 1e-4)
		end
	end

end
