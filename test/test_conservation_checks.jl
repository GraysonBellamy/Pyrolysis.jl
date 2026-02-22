"""
Conservation checks validation tests.

Tests mass and energy conservation tracking and verification.

Reference: Physics Audit Recommendations (2026)
"""

using Test
using Pyrolysis
using Pyrolysis.PostProcess: MassConservationTracker, EnergyConservationTracker,
	compute_domain_mass, compute_sensible_energy,
	update_mass_tracker!, update_energy_tracker!,
	check_mass_conservation, check_energy_conservation,
	mass_conservation_report, energy_conservation_report

@testset "Conservation Checks" begin

	@testset "Mass Conservation Tracker" begin
		@testset "Initialization and Domain Mass" begin
			# Create simple test problem
			virgin = SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2, ε = 0.9)
			gas = GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.025, λ = 1e-5)
			rxn = Reaction(:decomp, 1 => 2, A = 1e10, E = 100e3, h = 100e3)
			material = Material(name = :test, components = (virgin, gas), reactions = (rxn,))

			n_cells = 10
			L = 0.01
			mesh = create_uniform_mesh(L, n_cells, Val(2); T_initial = 300.0, ξ_initial = (1000.0, 0.0))

			# Initialize tracker
			tracker = MassConservationTracker(mesh, material)

			# Check initial mass calculation
			expected_mass_virgin = 1000.0 * L * mesh.cross_section_area  # ρ × V × A
			expected_mass_gas = 0.0

			@test tracker.initial_mass[1] ≈ expected_mass_virgin
			@test tracker.initial_mass[2] ≈ expected_mass_gas
			@test tracker.cumulative_efflux == (0.0, 0.0)
			@test tracker.last_time == 0.0

			# Verify compute_domain_mass
			current_mass = compute_domain_mass(mesh, material)
			@test current_mass[1] ≈ expected_mass_virgin
			@test current_mass[2] ≈ expected_mass_gas
		end

		@testset "Mass Tracking with Boundary Fluxes" begin
			virgin = SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2)
			gas = GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.025, λ = 1e-5)
			material = Material(name = :test, components = (virgin, gas), reactions = tuple())

			mesh = create_uniform_mesh(0.01, 10, Val(2); T_initial = 300.0, ξ_initial = (1000.0, 0.0))
			tracker = MassConservationTracker(mesh, material)

			# Simulate boundary fluxes (gas escaping at top, nothing at bottom)
			A = mesh.cross_section_area
			dt = 0.1
			J_boundary = [0.0 0.0;    # virgin (solid, no flux)
			              -0.001 0.0]  # gas (escaping at top)

			update_mass_tracker!(tracker, J_boundary, A, dt)

			# Expected efflux: -(-0.001) * A * dt for gas
			expected_efflux_gas = 0.001 * A * dt
			@test tracker.cumulative_efflux[1] ≈ 0.0
			@test tracker.cumulative_efflux[2] ≈ expected_efflux_gas
			@test tracker.last_time ≈ dt

			# Simulate more time steps
			for _ in 1:10
				update_mass_tracker!(tracker, J_boundary, A, dt)
			end

			@test tracker.cumulative_efflux[2] ≈ 0.001 * A * dt * 11  # 11 total steps
			@test tracker.last_time ≈ 1.1
		end

		@testset "Mass Conservation Check - Perfect Conservation" begin
			virgin = SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2)
			gas = GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.025, λ = 1e-5)
			material = Material(name = :test, components = (virgin, gas), reactions = tuple())

			mesh = create_uniform_mesh(0.01, 10, Val(2); T_initial = 300.0, ξ_initial = (1000.0, 10.0))
			tracker = MassConservationTracker(mesh, material)

			# Simulate gas leaving
			A = mesh.cross_section_area
			dt = 0.1
			J_boundary = [0.0 0.0; -0.001 0.0]

			# Manually remove gas from domain to match flux
			for _ in 1:10
				update_mass_tracker!(tracker, J_boundary, A, dt)
				mass_removed_per_cell = 0.001 * A * dt / length(mesh.cells)
				for i in eachindex(mesh.cells)
					old_ξ = mesh.cell_states[i].ξ
					new_ξ2 = max(0.0, old_ξ[2] - mass_removed_per_cell / mesh.cell_volumes[i])
					mesh.cell_states[i] = Pyrolysis.Geometry.CellState(
						mesh.cell_states[i].T,
						(old_ξ[1], new_ξ2),
						mesh.cell_states[i].volume,
						mesh.cell_states[i].volume_ref
					)
				end
			end

			# Check conservation
			result = check_mass_conservation(tracker, mesh, material)
			@test result.conserved == true
			@test result.max_error < 1e-10

			# Gas component should have error near zero
			@test result.relative_error[2] < 1e-10
		end

		@testset "Mass Conservation Report" begin
			virgin = SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2)
			gas = GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.025, λ = 1e-5)
			material = Material(name = :test, components = (virgin, gas), reactions = tuple())

			mesh = create_uniform_mesh(0.01, 10, Val(2); T_initial = 300.0, ξ_initial = (1000.0, 10.0))
			tracker = MassConservationTracker(mesh, material)

			# Generate report
			report = mass_conservation_report(tracker, mesh, material)

			# Check report contains key information
			@test occursin("MASS CONSERVATION REPORT", report)
			@test occursin("virgin", report)
			@test occursin("gas", report)
			@test occursin("Initial mass:", report)
			@test occursin("PASSED", report) || occursin("FAILED", report)
		end
	end

	@testset "Energy Conservation Tracker" begin
		@testset "Initialization and Sensible Energy" begin
			virgin = SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2, ε = 0.9)
			gas = GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.025, λ = 1e-5)
			material = Material(name = :test, components = (virgin, gas), reactions = tuple())

			n_cells = 10
			L = 0.01
			T_init = 400.0
			mesh = create_uniform_mesh(L, n_cells, Val(2); T_initial = T_init, ξ_initial = (1000.0, 0.0))

			# Initialize tracker
			tracker = EnergyConservationTracker(mesh, material)

			# Check initial energy calculation
			# E = ρcₚT × V × A (for single component, uniform T)
			expected_energy = 1000.0 * 1500.0 * T_init * L * mesh.cross_section_area

			@test tracker.initial_energy ≈ expected_energy rtol=0.01
			@test tracker.cumulative_boundary_flux == 0.0
			@test tracker.cumulative_reaction_heat == 0.0
			@test tracker.last_time == 0.0

			# Verify compute_sensible_energy
			current_energy = compute_sensible_energy(mesh, material)
			@test current_energy ≈ expected_energy rtol=0.01
		end

		@testset "Energy Tracking with Sources" begin
			virgin = SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2)
			gas = GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.025, λ = 1e-5)
			material = Material(name = :test, components = (virgin, gas), reactions = tuple())

			mesh = create_uniform_mesh(0.01, 10, Val(2); T_initial = 300.0, ξ_initial = (1000.0, 0.0))
			tracker = EnergyConservationTracker(mesh, material)

			# Simulate energy sources
			A = mesh.cross_section_area
			dt = 0.1
			q_boundary = 5000.0  # 5 kW/m² into domain
			Q_rxn_total = -1000.0  # -1 kW reaction heat (endothermic)

			update_energy_tracker!(tracker, q_boundary, Q_rxn_total, A, dt)

			# Expected cumulative values
			@test tracker.cumulative_boundary_flux ≈ q_boundary * A * dt
			@test tracker.cumulative_reaction_heat ≈ Q_rxn_total * dt
			@test tracker.last_time ≈ dt

			# Multiple time steps
			for _ in 1:10
				update_energy_tracker!(tracker, q_boundary, Q_rxn_total, A, dt)
			end

			@test tracker.cumulative_boundary_flux ≈ q_boundary * A * dt * 11
			@test tracker.cumulative_reaction_heat ≈ Q_rxn_total * dt * 11
			@test tracker.last_time ≈ 1.1
		end

		@testset "Energy Conservation Check - Perfect Conservation" begin
			virgin = SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2)
			gas = GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.025, λ = 1e-5)
			material = Material(name = :test, components = (virgin, gas), reactions = tuple())

			mesh = create_uniform_mesh(0.01, 10, Val(2); T_initial = 300.0, ξ_initial = (1000.0, 0.0))
			tracker = EnergyConservationTracker(mesh, material)

			# Add energy to domain via boundary flux
			A = mesh.cross_section_area
			dt = 0.1
			q_boundary = 5000.0
			Q_rxn_total = 0.0

			for _ in 1:10
				update_energy_tracker!(tracker, q_boundary, Q_rxn_total, A, dt)

				# Manually increase temperature to match energy input
				# ΔE = q × A × dt
				# ΔE = ρ × c × ΔT × V (for entire domain)
				total_volume = sum(mesh.cell_volumes)
				ρcp = 1000.0 * 1500.0
				ΔT = (q_boundary * A * dt) / (ρcp * total_volume * mesh.cross_section_area)

				for i in eachindex(mesh.cells)
					old_T = mesh.cell_states[i].T
					mesh.cell_states[i] = Pyrolysis.Geometry.CellState(
						old_T + ΔT,
						mesh.cell_states[i].ξ,
						mesh.cell_states[i].volume,
						mesh.cell_states[i].volume_ref
					)
				end
			end

			# Check conservation
			result = check_energy_conservation(tracker, mesh, material)
			@test result.conserved == true
			@test result.relative_error < 1e-6
		end

		@testset "Energy Conservation Report" begin
			virgin = SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2)
			gas = GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.025, λ = 1e-5)
			material = Material(name = :test, components = (virgin, gas), reactions = tuple())

			mesh = create_uniform_mesh(0.01, 10, Val(2); T_initial = 300.0, ξ_initial = (1000.0, 0.0))
			tracker = EnergyConservationTracker(mesh, material)

			# Generate report
			report = energy_conservation_report(tracker, mesh, material)

			# Check report contains key information
			@test occursin("ENERGY CONSERVATION REPORT", report)
			@test occursin("Initial energy:", report)
			@test occursin("Boundary flux", report)
			@test occursin("Reaction heat", report)
			@test occursin("PASSED", report) || occursin("FAILED", report)
		end
	end

	@testset "Conservation with Realistic Simulation" begin
		# Test conservation tracking in a small realistic pyrolysis simulation
		@testset "Simple Pyrolysis with Gas Generation" begin
			# Create material
			virgin = SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86)
			gas = GaseousComponent(:volatiles, M = 0.100, c = 1500.0, k = 0.02, λ = 2e-5)
			rxn = Reaction(:decomp, 1 => 2, A = 1e10, E = 150e3, h = 500e3)
			material = Material(name = :PMMA, components = (virgin, gas), reactions = (rxn,))

			# Create mesh
			mesh = create_uniform_mesh(0.005, 10, Val(2); T_initial = 300.0, ξ_initial = (1190.0, 0.0))

			# Initialize trackers
			mass_tracker = MassConservationTracker(mesh, material)
			energy_tracker = EnergyConservationTracker(mesh, material)

			# Store initial values
			initial_mass_total = sum(mass_tracker.initial_mass)
			initial_energy = energy_tracker.initial_energy

			# The simulation would update these trackers at each step
			# For this test, just verify they're initialized correctly
			@test mass_tracker.initial_mass[1] > 0.0  # Virgin material present
			@test mass_tracker.initial_mass[2] ≈ 0.0  # No initial gas
			@test energy_tracker.initial_energy > 0.0

			# Verify conservation checks work
			mass_result = check_mass_conservation(mass_tracker, mesh, material)
			energy_result = check_energy_conservation(energy_tracker, mesh, material)

			# Initially perfect conservation (no time evolution yet)
			@test mass_result.conserved == true
			@test energy_result.conserved == true
		end
	end
end
