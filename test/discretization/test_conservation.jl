"""
Conservation checks validation tests.

Tests mass and energy conservation tracking and verification.
"""

using Test
using Pyrolysis
using Pyrolysis.Diagnostics:
    MassConservationTracker,
    EnergyConservationTracker,
    compute_domain_mass,
    compute_sensible_energy,
    update_mass_tracker!,
    update_energy_tracker!,
    check_mass_conservation,
    check_energy_conservation,
    mass_conservation_report,
    energy_conservation_report,
    DiagnosticsState,
    record_substep!

@testset "Conservation Checks" begin

    @testset "Mass Conservation Tracker" begin
        @testset "Initialization and Domain Mass" begin
            # Create simple test problem
            virgin = SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2, ε = 0.9)
            gas = GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.025, λ = 1e-5)
            rxn = Reaction(:decomp, 1 => 2, A = 1e10, E = 100e3, h = 100e3)
            material =
                Material(name = :test, components = (virgin, gas), reactions = (rxn,))

            n_cells = 10
            L = 0.01
            mesh = create_uniform_mesh(
                L,
                n_cells,
                Val(2);
                T_initial = 300.0,
                ξ_initial = (1000.0, 0.0),
            )

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
            material =
                Material(name = :test, components = (virgin, gas), reactions = tuple())

            mesh = create_uniform_mesh(
                0.01,
                10,
                Val(2);
                T_initial = 300.0,
                ξ_initial = (1000.0, 0.0),
            )
            tracker = MassConservationTracker(mesh, material)

            # Simulate boundary fluxes (gas escaping at top, nothing at bottom)
            A = mesh.cross_section_area
            dt = 0.1
            J_boundary = [
                0.0 0.0;    # virgin (solid, no flux)
                -0.001 0.0
            ]  # gas (escaping at top)

            update_mass_tracker!(tracker, J_boundary, A, dt)

            # Expected efflux: -(-0.001) * A * dt for gas
            expected_efflux_gas = 0.001 * A * dt
            @test tracker.cumulative_efflux[1] ≈ 0.0
            @test tracker.cumulative_efflux[2] ≈ expected_efflux_gas
            @test tracker.last_time ≈ dt

            # Simulate more time steps
            for _ = 1:10
                update_mass_tracker!(tracker, J_boundary, A, dt)
            end

            @test tracker.cumulative_efflux[2] ≈ 0.001 * A * dt * 11  # 11 total steps
            @test tracker.last_time ≈ 1.1
        end

        @testset "Mass Conservation Check - Perfect Conservation" begin
            virgin = SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2)
            gas = GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.025, λ = 1e-5)
            material =
                Material(name = :test, components = (virgin, gas), reactions = tuple())

            mesh = create_uniform_mesh(
                0.01,
                10,
                Val(2);
                T_initial = 300.0,
                ξ_initial = (1000.0, 10.0),
            )
            tracker = MassConservationTracker(mesh, material)

            # Simulate gas leaving
            A = mesh.cross_section_area
            dt = 0.1
            J_boundary = [0.0 0.0; -0.001 0.0]

            # Manually remove gas from domain to match flux
            for _ = 1:10
                update_mass_tracker!(tracker, J_boundary, A, dt)
                mass_removed_per_cell = 0.001 * A * dt / length(mesh.cells)
                for i in eachindex(mesh.cells)
                    old_ξ = mesh.cell_states[i].ξ
                    new_ξ2 =
                        max(0.0, old_ξ[2] - mass_removed_per_cell / mesh.cell_volumes[i])
                    mesh.cell_states[i] = Pyrolysis.Geometry.CellState(
                        mesh.cell_states[i].T,
                        (old_ξ[1], new_ξ2),
                        mesh.cell_states[i].volume,
                        mesh.cell_states[i].volume_ref,
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
            material =
                Material(name = :test, components = (virgin, gas), reactions = tuple())

            mesh = create_uniform_mesh(
                0.01,
                10,
                Val(2);
                T_initial = 300.0,
                ξ_initial = (1000.0, 10.0),
            )
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
            material =
                Material(name = :test, components = (virgin, gas), reactions = tuple())

            n_cells = 10
            L = 0.01
            T_init = 400.0
            mesh = create_uniform_mesh(
                L,
                n_cells,
                Val(2);
                T_initial = T_init,
                ξ_initial = (1000.0, 0.0),
            )

            # Initialize tracker
            tracker = EnergyConservationTracker(mesh, material)

            # ρcp sums every component; with no gas present the initial
            # energy is the solid contribution.
            expected = 1000.0 * 1500.0 * T_init * L * mesh.cross_section_area

            @test tracker.initial_energy ≈ expected rtol=0.01
            @test tracker.cumulative_boundary_conduction == 0.0
            @test tracker.cumulative_reaction_heat == 0.0
            @test tracker.cumulative_advection_sink == 0.0
            @test tracker.last_time == 0.0

            @test compute_sensible_energy(mesh, material) ≈ expected rtol=0.01

            # Pore gas stores sensible heat like every other component.
            mesh_gas = create_uniform_mesh(
                L,
                n_cells,
                Val(2);
                T_initial = T_init,
                ξ_initial = (1000.0, 2.0),
            )
            expected_gas =
                (1000.0 * 1500.0 + 2.0 * 1000.0) * T_init * L * mesh.cross_section_area
            @test compute_sensible_energy(mesh_gas, material) ≈ expected_gas rtol=0.01
        end

        @testset "Energy Tracking with Sources" begin
            virgin = SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2)
            gas = GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.025, λ = 1e-5)
            material =
                Material(name = :test, components = (virgin, gas), reactions = tuple())

            mesh = create_uniform_mesh(
                0.01,
                10,
                Val(2);
                T_initial = 300.0,
                ξ_initial = (1000.0, 0.0),
            )
            tracker = EnergyConservationTracker(mesh, material)

            # Simulate one substep's contributions across the ledger fields.
            A = mesh.cross_section_area
            dt = 0.1
            q_cond_boundary = 5000.0   # conducted into slab
            Q_rxn_total = -1000.0  # endothermic reaction heat
            S_conv_total = -200.0   # gas advection sink

            update_energy_tracker!(
                tracker;
                q_cond_boundary = q_cond_boundary,
                Q_rxn_total = Q_rxn_total,
                S_conv_total = S_conv_total,
                A = A,
                dt = dt,
            )

            @test tracker.cumulative_boundary_conduction ≈ q_cond_boundary * A * dt
            @test tracker.cumulative_reaction_heat ≈ Q_rxn_total * dt
            @test tracker.cumulative_advection_sink ≈ S_conv_total * dt
            @test tracker.last_time ≈ dt

            # Multiple substeps accumulate linearly.
            for _ = 1:10
                update_energy_tracker!(
                    tracker;
                    q_cond_boundary = q_cond_boundary,
                    Q_rxn_total = Q_rxn_total,
                    S_conv_total = S_conv_total,
                    A = A,
                    dt = dt,
                )
            end

            @test tracker.cumulative_boundary_conduction ≈ q_cond_boundary * A * dt * 11
            @test tracker.cumulative_reaction_heat ≈ Q_rxn_total * dt * 11
            @test tracker.last_time ≈ 1.1
        end

        @testset "Diagnostics duplicate time does not advance stride" begin
            solid = SolidComponent(:solid, ρ = 1000.0, c = 1500.0, k = 0.2)
            material = Material(name = :test, components = (solid,), reactions = tuple())
            mesh = create_uniform_mesh(
                0.01,
                2,
                Val(1);
                T_initial = 300.0,
                ξ_initial = (1000.0,),
            )
            state = DiagnosticsState(mesh, material; stride = 2)
            J_boundary = zeros(1, 2)

            function record_at!(t)
                record_substep!(
                    state,
                    mesh,
                    material;
                    q_cond_boundary = 0.0,
                    Q_rxn_total = 0.0,
                    S_conv_total = 0.0,
                    J_boundary = J_boundary,
                    A = mesh.cross_section_area,
                    t = t,
                )
            end

            record_at!(0.0)
            @test state.step_counter == 0
            @test isempty(state.log)

            record_at!(0.1)
            @test state.step_counter == 1
            @test isempty(state.log)

            record_at!(0.2)
            @test state.step_counter == 2
            @test length(state.log) == 1
        end

        @testset "Energy Conservation Check - Perfect Closure" begin
            virgin = SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2)
            gas = GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.025, λ = 1e-5)
            material =
                Material(name = :test, components = (virgin, gas), reactions = tuple())

            mesh = create_uniform_mesh(
                0.01,
                10,
                Val(2);
                T_initial = 300.0,
                ξ_initial = (1000.0, 0.0),
                cross_section_area = 0.25,
            )
            tracker = EnergyConservationTracker(mesh, material)

            # Add energy to domain matching what the tracker accumulates as
            # boundary conduction + advection sink.
            A = mesh.cross_section_area
            dt = 0.1
            q_cond_boundary = 5000.0

            for _ = 1:10
                update_energy_tracker!(
                    tracker;
                    q_cond_boundary = q_cond_boundary,
                    Q_rxn_total = 0.0,
                    S_conv_total = 0.0,
                    A = A,
                    dt = dt,
                )

                # Match the energy delivered to the tracker: ΔE = q · A · dt
                total_volume = sum(mesh.cell_volumes)
                ρcp = 1000.0 * 1500.0   # no gas present in this mesh
                ΔT = (q_cond_boundary * A * dt) / (ρcp * total_volume)

                for i in eachindex(mesh.cells)
                    old_T = mesh.cell_states[i].T
                    mesh.cell_states[i] = Pyrolysis.Geometry.CellState(
                        old_T + ΔT,
                        mesh.cell_states[i].ξ,
                        mesh.cell_states[i].volume,
                        mesh.cell_states[i].volume_ref,
                    )
                end
            end

            # The ledger residual must close to ~zero (constant ρcp and V).
            result = check_energy_conservation(tracker, mesh, material)
            @test result.conserved
            @test result.relative < 1e-6
        end

        @testset "Energy Conservation Report" begin
            virgin = SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2)
            gas = GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.025, λ = 1e-5)
            material =
                Material(name = :test, components = (virgin, gas), reactions = tuple())

            mesh = create_uniform_mesh(
                0.01,
                10,
                Val(2);
                T_initial = 300.0,
                ξ_initial = (1000.0, 0.0),
            )
            tracker = EnergyConservationTracker(mesh, material)

            report = energy_conservation_report(tracker, mesh, material)

            @test occursin("ENERGY CONSERVATION REPORT", report)
            @test occursin("E(0):", report)
            @test occursin("Residual", report)
        end
    end

    @testset "Conservation with Realistic Simulation" begin
        # Test conservation tracking in a small realistic pyrolysis simulation
        @testset "Simple Pyrolysis with Gas Generation" begin
            # Create material
            virgin = SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86)
            gas = GaseousComponent(:volatiles, M = 0.100, c = 1500.0, k = 0.02, λ = 2e-5)
            rxn = Reaction(:decomp, 1 => 2, A = 1e10, E = 150e3, h = 500e3)
            material =
                Material(name = :PMMA, components = (virgin, gas), reactions = (rxn,))

            # Create mesh
            mesh = create_uniform_mesh(
                0.005,
                10,
                Val(2);
                T_initial = 300.0,
                ξ_initial = (1190.0, 0.0),
            )

            # Initialize trackers
            mass_tracker = MassConservationTracker(mesh, material)
            energy_tracker = EnergyConservationTracker(mesh, material)

            # Store initial values
            initial_mass_total = sum(mass_tracker.initial_mass)

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

    # Under variable cross-section, trackers receive the
    # trapezoidal-in-A average over the substep. The trackers themselves take
    # `A::Float64` per call — these tests pin that a *time-varying* A passed
    # substep-by-substep still closes the books exactly.
    @testset "Variable cross-section — time-varying A passes through trackers" begin
        virgin = SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2)
        gas = GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.025, λ = 1e-5)
        material = Material(name = :test, components = (virgin, gas), reactions = tuple())

        mesh = create_uniform_mesh(
            0.01,
            10,
            Val(2);
            T_initial = 300.0,
            ξ_initial = (1000.0, 10.0),
        )
        tracker = MassConservationTracker(mesh, material)

        # Hand-driven schedule: A shrinks linearly from A_0 to 0.85·A_0 over
        # 10 substeps. Each substep uses the trapezoidal midpoint A_mid.
        A_0 = mesh.cross_section_area_initial
        A_schedule = collect(range(A_0, 0.85 * A_0, length = 11))
        dt = 0.1
        J_boundary = [0.0 0.0; -0.001 0.0]

        expected_efflux = 0.0
        for k = 1:10
            A_mid = 0.5 * (A_schedule[k] + A_schedule[k+1])
            update_mass_tracker!(tracker, J_boundary, A_mid, dt)
            expected_efflux += 0.001 * A_mid * dt
        end

        # Closure under varying A: cumulative efflux equals the analytic
        # trapezoidal sum, independent of the choice of A schedule.
        @test tracker.cumulative_efflux[2] ≈ expected_efflux rtol = 1e-12
        @test tracker.last_time ≈ 1.0
    end
end
