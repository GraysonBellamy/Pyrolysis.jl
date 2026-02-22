"""
Tests for Problem and Solver Interface

Tests:
1. PyrolysisProblem creation and validation
2. PyrolysisSolution structure
3. Initial condition application
4. Basic solve functionality
"""

using Test
using Pyrolysis

@testset "Problem and Solver Interface" begin

	# ========================================================================
	# Helper: Create a simple test material
	# ========================================================================
	function create_test_material()
		virgin = SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86)
		gas = GaseousComponent(:MMA, M = 0.100, c = 1500, k = 0.02, λ = 1e-5)
		rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)

		Material(
			name = :PMMA,
			components = (virgin, gas),
			reactions = (rxn,),
		)
	end

	# ========================================================================
	@testset "PyrolysisProblem Creation" begin
		material = create_test_material()
		mesh = create_uniform_mesh(0.01, 20, Val(2))

		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(50e3),
			thermal_bottom = AdiabaticBC(),
		)

		# Basic problem creation
		@testset "Keyword constructor" begin
			problem = PyrolysisProblem(
				mesh = mesh,
				material = material,
				bc_set = bc_set,
				tspan = (0.0, 10.0),
			)

			@test problem.mesh === mesh
			@test problem.material === material
			@test problem.tspan == (0.0, 10.0)
			@test n_components(problem.material) == 2
		end

		@testset "Positional constructor" begin
			problem = PyrolysisProblem(mesh, material, bc_set, (0.0, 10.0))
			@test problem.tspan == (0.0, 10.0)
		end

		@testset "Default mass BCs" begin
			problem = PyrolysisProblem(
				mesh = mesh,
				material = material,
				bc_set = bc_set,
				tspan = (0.0, 10.0),
			)

			# BoundaryConditionSet defaults mass BCs to ImpermeableBC
			@test problem.bc_set.mass_top isa ImpermeableBC
			@test problem.bc_set.mass_bottom isa ImpermeableBC
		end

		@testset "Custom initial conditions" begin
			problem = PyrolysisProblem(
				mesh = mesh,
				material = material,
				bc_set = bc_set,
				T_initial = z -> 400.0 + 100.0 * z,
				ξ_initial = [z -> 1000.0, z -> 0.0],
				tspan = (0.0, 10.0),
			)

			@test problem.T_initial(0.0) ≈ 400.0
			@test problem.T_initial(0.01) ≈ 401.0
			@test problem.ξ_initial[1](0.0) ≈ 1000.0
		end

		@testset "Auto-generated output times" begin
			problem = PyrolysisProblem(
				mesh = mesh,
				material = material,
				bc_set = bc_set,
				tspan = (0.0, 100.0),
			)

			@test !isempty(problem.output_times)
			@test problem.output_times[1] ≈ 0.0
			@test problem.output_times[end] ≈ 100.0
		end

		@testset "Display methods" begin
			problem = PyrolysisProblem(mesh, material, bc_set, (0.0, 10.0))

			# Should not throw
			str = sprint(show, problem)
			@test occursin("PyrolysisProblem", str)
			@test occursin("PMMA", str)

			str_full = sprint(show, MIME("text/plain"), problem)
			@test occursin("Material:", str_full)
			@test occursin("Components:", str_full)
		end
	end

	# ========================================================================
	@testset "Problem Validation" begin
		material = create_test_material()

		@testset "Valid problem passes" begin
			mesh = create_uniform_mesh(0.01, 20, Val(2))
			bc_set = BoundaryConditionSet(thermal_top = HeatFluxBC(50e3), thermal_bottom = AdiabaticBC())
			problem = PyrolysisProblem(mesh, material, bc_set, (0.0, 10.0))

			@test validate_problem(problem) == true
			@test is_valid_problem(problem) == true
		end

		@testset "BoundaryConditionSet defaults" begin
			mesh = create_uniform_mesh(0.01, 20, Val(2))
			# Only provide thermal_top BC - others should default
			bc_set = BoundaryConditionSet(thermal_top = HeatFluxBC(50e3))

			problem = PyrolysisProblem(
				mesh = mesh,
				material = material,
				bc_set = bc_set,
				tspan = (0.0, 10.0),
			)

			# BoundaryConditionSet always has all BCs - missing ones filled with defaults
			@test problem.bc_set.thermal_top isa HeatFluxBC
			@test problem.bc_set.thermal_bottom isa AdiabaticBC  # Default
			@test problem.bc_set.mass_top isa ImpermeableBC  # Default
			@test problem.bc_set.mass_bottom isa ImpermeableBC  # Default
			@test validate_problem(problem) == true
			@test is_valid_problem(problem) == true
		end

		@testset "Invalid tspan detected" begin
			mesh = create_uniform_mesh(0.01, 20, Val(2))
			bc_set = BoundaryConditionSet(
				thermal_top = HeatFluxBC(50e3),
				thermal_bottom = AdiabaticBC(),
				mass_top = ImpermeableBC(),
				mass_bottom = ImpermeableBC(),
			)

			# t_start >= t_end - construct with invalid tspan
			problem = PyrolysisProblem(
				mesh = mesh,
				material = material,
				bc_set = bc_set,
				T_initial = z -> 300.0,
				ξ_initial = [z -> 1190.0, z -> 0.0],
				external_flux = t -> 0.0,
				ambient_temperature = t -> 300.0,
				tspan = (10.0, 5.0),
				dt_initial = 1e-4,
				output_times = Float64[],
			)

			@test_throws ErrorException validate_problem(problem)
		end
	end

	# ========================================================================
	@testset "Initial Conditions" begin
		material = create_test_material()
		mesh = create_uniform_mesh(0.01, 10, Val(2))

		@testset "Apply from problem" begin
			bc_set = BoundaryConditionSet(thermal_top = HeatFluxBC(50e3), thermal_bottom = AdiabaticBC())
			problem = PyrolysisProblem(
				mesh = mesh,
				material = material,
				bc_set = bc_set,
				T_initial = z -> 350.0,
				ξ_initial = [z -> 1000.0, z -> 100.0],
				tspan = (0.0, 10.0),
			)

			set_initial_conditions!(problem)

			# Check first cell
			@test mesh.cell_states[1].T ≈ 350.0
			@test mesh.cell_states[1].ξ[1] ≈ 1000.0
			@test mesh.cell_states[1].ξ[2] ≈ 100.0
		end

		@testset "Uniform initial conditions" begin
			mesh2 = create_uniform_mesh(0.01, 10, Val(2))
			apply_initial_conditions!(mesh2, material; T = 400.0, ξ = (500.0, 50.0))

			for state in mesh2.cell_states
				@test state.T ≈ 400.0
				@test state.ξ[1] ≈ 500.0
				@test state.ξ[2] ≈ 50.0
			end
		end

		@testset "Profile helpers" begin
			# Step profile
			@test step_profile(0.003; z_interface = 0.005, value_left = 100.0, value_right = 200.0) ≈ 100.0
			@test step_profile(0.007; z_interface = 0.005, value_left = 100.0, value_right = 200.0) ≈ 200.0

			# Linear profile
			@test linear_profile(0.0; z_min = 0.0, z_max = 0.01, value_min = 300.0, value_max = 400.0) ≈ 300.0
			@test linear_profile(0.005; z_min = 0.0, z_max = 0.01, value_min = 300.0, value_max = 400.0) ≈ 350.0
			@test linear_profile(0.01; z_min = 0.0, z_max = 0.01, value_min = 300.0, value_max = 400.0) ≈ 400.0

			# Gaussian profile
			val = gaussian_profile(0.005; z_center = 0.005, width = 0.001, amplitude = 100.0, baseline = 300.0)
			@test val ≈ 400.0  # At center, should be baseline + amplitude

			# Tanh profile
			val_left = tanh_profile(0.0; z_center = 0.005, width = 0.001, value_left = 100.0, value_right = 200.0)
			val_right = tanh_profile(0.01; z_center = 0.005, width = 0.001, value_left = 100.0, value_right = 200.0)
			@test val_left < 110.0  # Close to left value
			@test val_right > 190.0  # Close to right value
		end
	end

	# ========================================================================
	@testset "PyrolysisSolution" begin
		material = create_test_material()
		mesh = create_uniform_mesh(0.01, 10, Val(2))
		bc_set = BoundaryConditionSet(thermal_top = HeatFluxBC(50e3), thermal_bottom = AdiabaticBC())
		problem = PyrolysisProblem(mesh, material, bc_set, (0.0, 10.0))

		@testset "Empty solution creation" begin
			sol = PyrolysisSolution(problem)
			@test sol.retcode == :Default
			@test isempty(sol.times)
			@test sol.solve_time == 0.0
		end

		@testset "Display methods" begin
			sol = PyrolysisSolution(problem)
			str = sprint(show, sol)
			@test occursin("PyrolysisSolution", str)
		end
	end

	# ========================================================================
	@testset "Solve Interface" begin
		# Create a simple heat conduction problem (no reactions)
		virgin = SolidComponent(:virgin, ρ = 1000, c = 1000, k = 1.0, ε = 0.9)
		simple_material = Material(
			name = :Simple,
			components = (virgin,),
			reactions = (),
		)

		mesh = create_uniform_mesh(0.01, 10, Val(1); T_initial = 300.0)
		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(10000.0),  # 10 kW/m² at top
			thermal_bottom = AdiabaticBC(),
		)

		problem = PyrolysisProblem(
			mesh = mesh,
			material = simple_material,
			bc_set = bc_set,
			T_initial = z -> 300.0,
			ξ_initial = [z -> 1000.0],
			tspan = (0.0, 1.0),
			dt_initial = 1e-3,
			output_times = [0.0, 0.5, 1.0],
		)

		@testset "Basic solve" begin
			sol = solve(problem; progress = false)

			@test sol.retcode == :Success
			@test length(sol.times) >= 3
			@test sol.solve_time > 0.0

			# Temperature should increase at surface
			@test sol.surface_temperature[end] > sol.surface_temperature[1]

			# Check solution structure
			@test !isempty(sol.temperatures)
			@test !isempty(sol.concentrations)
			@test length(sol.temperatures) == length(sol.times)
		end

		@testset "Solve with different tolerance" begin
			sol = solve(problem; abstol = 1e-6, reltol = 1e-4, progress = false)
			@test sol.retcode == :Success
		end

		@testset "Get profile utility" begin
			sol = solve(problem; progress = false)

			T_profile = get_profile(sol, :T, 1)
			@test length(T_profile) == 10

			ξ_profile = get_profile(sol, :ξ1, 1)
			@test length(ξ_profile) == 10
		end

		@testset "Get cell positions" begin
			sol = solve(problem; progress = false)
			positions = get_cell_positions(sol)
			@test length(positions) == 10
			@test positions[1] < positions[end]
		end
	end

	# ========================================================================
	@testset "Full Integration Test" begin
		# Test complete workflow with reactions
		virgin = SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86)
		gas = GaseousComponent(:gas, M = 0.050, c = 1500, k = 0.02, λ = 1e-5)

		# Very fast reaction for testing (low activation energy)
		rxn = Reaction(:decomp, 1 => 2, A = 1e6, E = 50e3, h = -500e3)

		material = Material(
			name = :TestMaterial,
			components = (virgin, gas),
			reactions = (rxn,),
		)

		mesh = create_uniform_mesh(0.005, 20, Val(2))
		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(100e3),  # High flux for fast heating
			thermal_bottom = AdiabaticBC(),
			mass_top = ConvectiveMassBC(h_m = 0.1, ξ_inf = 0.0),  # Gas escapes at surface
			mass_bottom = ImpermeableBC(),
		)

		problem = PyrolysisProblem(
			mesh = mesh,
			material = material,
			bc_set = bc_set,
			T_initial = z -> 300.0,
			ξ_initial = [z -> 1190.0, z -> 0.0],
			tspan = (0.0, 0.5),
			dt_initial = 1e-5,
			output_times = collect(range(0.0, 0.5, length = 11)),
		)

		sol = solve(problem; progress = false, maxiters = 100000)

		@test sol.retcode == :Success

		# Temperature should increase (heat flux applied)
		@test sol.surface_temperature[end] > sol.surface_temperature[1]

		# Check reaction is occurring: solid concentration should decrease
		initial_solid = sum(sol.concentrations[1][1, :])
		final_solid = sum(sol.concentrations[end][1, :])
		@test final_solid < initial_solid  # Solid consumed by reaction

		# Gas should be generated (even if some escapes)
		# Check that gas has been produced at some point
		max_gas = maximum(sum(sol.concentrations[i][2, :]) for i in 1:length(sol.times))
		@test max_gas > 0.0  # Gas was produced
	end

end
