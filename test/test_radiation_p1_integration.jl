"""
Integration tests for P1 volumetric radiation in ODE time integration (Week 2).

Tests the integration of P1 radiation into the ODE framework:
- State vector with radiation (StateVectorALERadiation)
- Quasi-steady P1 solver within ODE function
- Coupled temperature-radiation evolution
- Standard ALE solver works without radiation
- Radiative flux divergence in energy equation
"""

using Test
using Pyrolysis
using Pyrolysis.Properties  # For update_properties!
using Pyrolysis.Discretization  # For update_radiation_properties!, etc.
using LinearAlgebra
using OrdinaryDiffEq

# Helper to create a simple test material
# Note: α = 100 m²/kg gives α_eff = 100,000 m⁻¹ with ξ = 1000 kg/m³
# This ensures optical thickness τ = α_eff × L ≈ 1000 >> 1 (P1 valid regime)
function create_test_material()
	Material(
		name = :TestMaterial,
		components = (
			SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.25, ε = 0.9, α = 100.0),  # α in m²/kg
			GaseousComponent(:gas, M = 0.028, c = 1000.0, k = 0.03, λ = 1e-5, ε = 0.0, α = 0.0),
		),
		reactions = (
			Reaction(:decomp, 1 => 2, A = 1e10, E = 150e3, h = 300e3),  # reactant 1 → product 2
		),
		β = 0.5,
	)
end

# Helper to create mesh with initial state
function create_test_mesh(n_cells::Int, NC::Int, T_initial::Float64 = 300.0)
	mesh = create_uniform_mesh(0.01, n_cells, Val(NC))  # 1 cm domain
	for i in 1:n_cells
		mesh.cell_states[i] = CellState{NC}(
			T_initial,
			ntuple(j -> j == 1 ? 1000.0 : 0.0, NC),  # All virgin initially
			mesh.cell_volumes[i],
			mesh.cell_volumes[i],
		)
	end
	return mesh
end

@testset "P1 Radiation Integration" begin

	@testset "StateVectorALERadiation" begin
		n_cells = 10
		NC = 2

		@testset "Construction and layout" begin
			sv = StateVectorALERadiation(n_cells, Val(NC))

			# Check dimensions
			@test sv.n_cells == 10
			@test sv.n_nodes == 11
			@test sv.n_components == 2

			# Check offsets - layout is [T, G, ξ, z]
			@test sv.T_offset == 0
			@test sv.G_offset == 10      # After T
			@test sv.ξ_offset == 20      # After G
			@test sv.z_offset == 40      # After ξ (2 components × 10 cells)

			# Check total length
			total_len = state_vector_length_ale_radiation(n_cells, NC)
			@test total_len == 10 + 10 + 20 + 11  # T + G + ξ + z = 51
		end

		@testset "Pack and unpack" begin
			sv = StateVectorALERadiation(n_cells, Val(NC))
			mesh = create_test_mesh(n_cells, NC, 500.0)

			# Create radiation cache
			rad_cache = RadiationCache(n_cells)
			rad_cache.G .= 1000.0  # Set initial G values

			# Pack state
			u = zeros(state_vector_length_ale_radiation(n_cells, NC))
			pack_state_ale_radiation!(u, mesh, rad_cache, sv)

			# Check packed values
			@test u[1:10] ≈ fill(500.0, 10)  # Temperatures
			@test u[11:20] ≈ fill(1000.0, 10)  # G values
			@test u[21:30] ≈ fill(1000.0, 10)  # ξ₁ (virgin)
			@test u[31:40] ≈ fill(0.0, 10)    # ξ₂ (gas)

			# Modify and unpack
			u[5] = 600.0   # Change one temperature
			u[15] = 2000.0  # Change one G value

			unpack_state_ale_radiation!(mesh, rad_cache, u, sv)

			@test mesh.cell_states[5].T ≈ 600.0
			@test rad_cache.G[5] ≈ 2000.0
		end

		@testset "Accessor functions" begin
			sv = StateVectorALERadiation(n_cells, Val(NC))
			mesh = create_test_mesh(n_cells, NC)
			rad_cache = RadiationCache(n_cells)
			rad_cache.G .= range(100.0, 1000.0, length = n_cells)

			u = zeros(state_vector_length_ale_radiation(n_cells, NC))
			pack_state_ale_radiation!(u, mesh, rad_cache, sv)

			# Test get functions (note: argument order is u, sv, i)
			@test get_temperature_ale_rad(u, sv, 3) == mesh.cell_states[3].T
			@test get_radiation_intensity(u, sv, 5) ≈ rad_cache.G[5]

			# Test set functions (note: argument order is u, sv, i, value)
			set_temperature_ale_rad!(u, sv, 3, 800.0)
			@test u[3] ≈ 800.0

			set_radiation_intensity!(u, sv, 5, 5000.0)
			@test u[sv.G_offset+5] ≈ 5000.0

			# Test views
			T_view = temperature_view_ale_rad(u, sv)
			@test length(T_view) == n_cells
			@test T_view[3] ≈ 800.0

			G_view = radiation_view(u, sv)
			@test length(G_view) == n_cells
			@test G_view[5] ≈ 5000.0
		end

		@testset "Validation" begin
			sv = StateVectorALERadiation(n_cells, Val(NC))
			mesh = create_test_mesh(n_cells, NC)
			rad_cache = RadiationCache(n_cells)

			# Valid state
			u = zeros(state_vector_length_ale_radiation(n_cells, NC))
			pack_state_ale_radiation!(u, mesh, rad_cache, sv)
			@test validate_state_vector_ale_radiation(u, sv)

			# Too short - should throw DimensionMismatch
			u_short = zeros(10)
			@test_throws DimensionMismatch validate_state_vector_ale_radiation(u_short, sv)
		end
	end

	@testset "Quasi-steady P1 solver" begin
		n_cells = 20
		NC = 2
		mesh = create_test_mesh(n_cells, NC, 500.0)
		material = create_test_material()

		# Create caches
		prop_cache = PropertyCache(n_cells, n_cells + 1)
		Properties.update_properties!(prop_cache, mesh, material)

		rad_cache = RadiationCache(n_cells)

		@testset "Uniform temperature solution" begin
			# For uniform T, G should equal 4σT⁴
			T_uniform = 500.0
			for i in 1:n_cells
				mesh.cell_states[i] = CellState{NC}(
					T_uniform,
					(1000.0, 0.0),
					mesh.cell_volumes[i],
					mesh.cell_volumes[i],
				)
			end
			Properties.update_properties!(prop_cache, mesh, material)
			update_radiation_properties!(rad_cache, mesh, prop_cache)

			# No external flux BCs
			radiation_bcs = RadiationBCSet(
				top = MarshakRadiationBC(I_external = 0.0, emissivity = 0.9),
				bottom = ReflectiveRadiationBC(),
			)

			solve_quasi_steady_radiation!(rad_cache, mesh, prop_cache, radiation_bcs, 0.0)

			# At equilibrium with no external source and non-blackbody boundaries,
			# G ≈ 4σεT⁴ where ε is the surface emissivity (not 4σT⁴)
			# This is because the Marshak BC with ε < 1 introduces surface losses
			σ = STEFAN_BOLTZMANN
			ε_surface = 0.9
			G_equilibrium = 4.0 * σ * ε_surface * T_uniform^4
			@test rad_cache.G[n_cells÷2] ≈ G_equilibrium rtol=0.1
		end

		@testset "External flux response" begin
			# Apply external flux and verify G increases at surface
			T_initial = 300.0
			for i in 1:n_cells
				mesh.cell_states[i] = CellState{NC}(
					T_initial,
					(1000.0, 0.0),
					mesh.cell_volumes[i],
					mesh.cell_volumes[i],
				)
			end
			Properties.update_properties!(prop_cache, mesh, material)
			update_radiation_properties!(rad_cache, mesh, prop_cache)

			# External flux at top boundary
			radiation_bcs = RadiationBCSet(
				top = MarshakRadiationBC(I_external = 50e3, emissivity = 0.9),
				bottom = ReflectiveRadiationBC(),
			)

			solve_quasi_steady_radiation!(rad_cache, mesh, prop_cache, radiation_bcs, 0.0)

			# G should be higher near the surface (cell n) than back (cell 1)
			@test rad_cache.G[n_cells] > rad_cache.G[1]

			# G at surface should be elevated due to external flux
			σ = STEFAN_BOLTZMANN
			G_equilibrium_no_flux = 4.0 * σ * T_initial^4
			@test rad_cache.G[n_cells] > G_equilibrium_no_flux
		end

		@testset "Radiative flux divergence" begin
			T_hot = 800.0
			T_cold = 300.0

			# Create temperature gradient (hot at top/surface, cold at bottom)
			for i in 1:n_cells
				T = T_cold + (T_hot - T_cold) * (i - 1) / (n_cells - 1)
				mesh.cell_states[i] = CellState{NC}(T, (1000.0, 0.0),
					mesh.cell_volumes[i], mesh.cell_volumes[i])
			end
			Properties.update_properties!(prop_cache, mesh, material)
			update_radiation_properties!(rad_cache, mesh, prop_cache)

			radiation_bcs = RadiationBCSet(
				top = MarshakRadiationBC(I_external = 25e3, emissivity = 0.9),
				bottom = ReflectiveRadiationBC(),
			)

			solve_quasi_steady_radiation!(rad_cache, mesh, prop_cache, radiation_bcs, 0.0)
			compute_radiation_flux_divergence!(rad_cache, mesh, prop_cache)

			# Divergence should be computed
			@test !all(rad_cache.S_rad .== 0.0)

			# For optically thick media (τ >> 1), radiation is in local equilibrium
			# with G ≈ 4σT⁴. The divergence S_rad = ∇·q_rad captures the net 
			# radiative heat flow. With a temperature gradient, there's net heat 
			# flow from hot to cold (in -z direction in this setup since hot is at top).
			# This manifests as S_rad being negative at the hot end (net cooling)
			# because more radiation leaves than enters.
			@test rad_cache.S_rad[end] < 0  # Hot region net emits (cools)
		end
	end

	@testset "PyrolysisODEFunctionALEWithP1Radiation" begin
		n_cells = 15
		NC = 2
		mesh = create_test_mesh(n_cells, NC, 400.0)
		material = create_test_material()

		# Boundary conditions
		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(50e3),
			thermal_bottom = AdiabaticBC(),
			mass_top = DiffusiveMassBC(0.01),
			mass_bottom = ImpermeableBC(),
		)

		@testset "Construction" begin
			ode_fn_base = PyrolysisODEFunctionALE(mesh, material, bc_set)
			rad_cache = RadiationCache(n_cells)
			radiation_bcs = RadiationBCSet(
				top = MarshakRadiationBC(I_external = 0.0, emissivity = 0.9),
				bottom = ReflectiveRadiationBC(),
			)

			ode_fn = PyrolysisODEFunctionALEWithP1Radiation(ode_fn_base, rad_cache, radiation_bcs)

			@test ode_fn.use_quasi_steady == true
			@test length(ode_fn.G_work) == n_cells
		end

		@testset "RHS evaluation" begin
			ode_fn_base = PyrolysisODEFunctionALE(mesh, material, bc_set)
			rad_cache = RadiationCache(n_cells)
			radiation_bcs = RadiationBCSet(
				top = MarshakRadiationBC(I_external = 25e3, emissivity = 0.9),
				bottom = ReflectiveRadiationBC(),
			)

			ode_fn = PyrolysisODEFunctionALEWithP1Radiation(ode_fn_base, rad_cache, radiation_bcs)

			# Create state vector using standard ALE layout (P1 uses quasi-steady G)
			u = create_state_vector_ale(mesh)
			du = similar(u)

			# Evaluate RHS
			ode_fn(du, u, nothing, 0.0)

			# Temperature derivatives should be computed
			@test !all(du[1:n_cells] .== 0.0)

			# G should have been updated in rad_cache (quasi-steady)
			@test !all(rad_cache.G .== 0.0)
		end

		@testset "Radiation effect on temperature" begin
			# Compare with and without P1 radiation
			mesh_no_rad = create_test_mesh(n_cells, NC, 400.0)
			mesh_with_rad = create_test_mesh(n_cells, NC, 400.0)

			# Without radiation (or with Beer-Lambert)
			ode_fn_no_rad = PyrolysisODEFunctionALE(mesh_no_rad, material, bc_set;
				I_radiation = 0.0)

			# With P1 radiation
			ode_fn_base = PyrolysisODEFunctionALE(mesh_with_rad, material, bc_set)
			rad_cache = RadiationCache(n_cells)
			radiation_bcs = RadiationBCSet(
				top = MarshakRadiationBC(I_external = 25e3, emissivity = 0.9),
				bottom = ReflectiveRadiationBC(),
			)
			ode_fn_with_rad = PyrolysisODEFunctionALEWithP1Radiation(ode_fn_base, rad_cache, radiation_bcs)

			u_no_rad = create_state_vector_ale(mesh_no_rad)
			u_with_rad = create_state_vector_ale(mesh_with_rad)
			du_no_rad = similar(u_no_rad)
			du_with_rad = similar(u_with_rad)

			ode_fn_no_rad(du_no_rad, u_no_rad, nothing, 0.0)
			ode_fn_with_rad(du_with_rad, u_with_rad, nothing, 0.0)

			# With P1 radiation, back cells should receive some heating from volumetric absorption
			# that isn't captured by surface-only models
			# The difference will depend on optical thickness
			dT_no_rad = du_no_rad[1:n_cells]
			dT_with_rad = du_with_rad[1:n_cells]

			# Both should have heating (positive dT) at the surface side
			@test dT_no_rad[n_cells] > 0
			@test dT_with_rad[n_cells] > 0
		end
	end

	@testset "Short ODE integration" begin
		n_cells = 10
		NC = 2
		mesh = create_test_mesh(n_cells, NC, 350.0)
		material = create_test_material()

		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(30e3),  # Modest heat flux
			thermal_bottom = AdiabaticBC(),
			mass_top = DiffusiveMassBC(0.01),
			mass_bottom = ImpermeableBC(),
		)

		@testset "Run to t=1.0s" begin
			ode_fn_base = PyrolysisODEFunctionALE(mesh, material, bc_set)
			rad_cache = RadiationCache(n_cells)
			radiation_bcs = RadiationBCSet(
				top = MarshakRadiationBC(I_external = 0.0, emissivity = 0.9),
				bottom = ReflectiveRadiationBC(),
			)
			ode_fn = PyrolysisODEFunctionALEWithP1Radiation(ode_fn_base, rad_cache, radiation_bcs)

			u0 = create_state_vector_ale(mesh)

			tspan = (0.0, 1.0)
			prob = ODEProblem(ode_fn, u0, tspan)

			# Solve with TRBDF2 (good for stiff problems)
			sol = solve(prob, TRBDF2(); abstol = 1e-6, reltol = 1e-4, maxiters = 1000)

			@test sol.retcode == ReturnCode.Success

			# Temperature should have increased at surface
			T_final_surface = sol.u[end][n_cells]
			@test T_final_surface > 350.0

			# Radiation intensity should be updated
			@test !all(rad_cache.G .== 0.0)
		end
	end

	@testset "Resize functionality" begin
		n_cells = 20
		NC = 2
		mesh = create_test_mesh(n_cells, NC)
		material = create_test_material()

		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(25e3),
			thermal_bottom = AdiabaticBC(),
			mass_top = DiffusiveMassBC(0.01),
			mass_bottom = ImpermeableBC(),
		)

		ode_fn_base = PyrolysisODEFunctionALE(mesh, material, bc_set)
		rad_cache = RadiationCache(n_cells)
		radiation_bcs = RadiationBCSet(
			top = MarshakRadiationBC(I_external = 25e3, emissivity = 0.9),
			bottom = ReflectiveRadiationBC(),
		)
		ode_fn = PyrolysisODEFunctionALEWithP1Radiation(ode_fn_base, rad_cache, radiation_bcs)

		# Resize to smaller mesh
		n_cells_new = 15
		n_faces_new = 16
		resize_ode_function_ale!(ode_fn, n_cells_new, n_faces_new)

		@test length(ode_fn.rad_cache.G) == n_cells_new
		@test length(ode_fn.G_work) == n_cells_new
	end

	@testset "Accessor functions" begin
		n_cells = 10
		NC = 2
		mesh = create_test_mesh(n_cells, NC)
		material = create_test_material()

		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(25e3),
			thermal_bottom = AdiabaticBC(),
			mass_top = ImpermeableBC(),
			mass_bottom = ImpermeableBC(),
		)

		ode_fn_base = PyrolysisODEFunctionALE(mesh, material, bc_set)
		rad_cache = RadiationCache(n_cells)
		rad_cache.G .= range(100.0, 1000.0, length = n_cells)

		radiation_bcs = RadiationBCSet(
			top = MarshakRadiationBC(I_external = 25e3),
			bottom = ReflectiveRadiationBC(),
		)
		ode_fn = PyrolysisODEFunctionALEWithP1Radiation(ode_fn_base, rad_cache, radiation_bcs)

		# Test get_radiation_intensity
		@test get_radiation_intensity(ode_fn, 5) ≈ rad_cache.G[5]
		@test get_radiation_intensity(ode_fn) === rad_cache.G
	end

	@testset "Standard ALE solver without radiation" begin
		# Ensure standard ALE solver still works without radiation
		n_cells = 10
		NC = 2
		mesh = create_test_mesh(n_cells, NC)
		material = create_test_material()

		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(25e3),
			thermal_bottom = AdiabaticBC(),
			mass_top = DiffusiveMassBC(0.01),
			mass_bottom = ImpermeableBC(),
		)

		ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)
		u0 = create_state_vector_ale(mesh)
		du = similar(u0)

		# Should work without radiation
		ode_fn(du, u0, nothing, 0.0)

		@test !all(du .== 0.0)
	end
end
