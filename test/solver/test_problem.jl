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
using Pyrolysis.Materials: n_components
using Pyrolysis.Problem:
    apply_initial_conditions!,
    validate_problem,
    is_valid_problem,
    step_profile,
    linear_profile,
    gaussian_profile,
    tanh_profile
using Pyrolysis.Solver: get_profile, get_cell_positions

@testset "Problem and Solver Interface" begin

    # ========================================================================
    # Helper: Create a simple test material
    # ========================================================================
    function create_test_material()
        virgin = SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86)
        gas = GaseousComponent(:MMA, M = 0.100, c = 1500, k = 0.02, λ = 1e-5)
        rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)

        Material(name = :PMMA, components = (virgin, gas), reactions = (rxn,))
    end

    # ========================================================================
    @testset "PyrolysisProblem Creation" begin
        material = create_test_material()
        mesh = create_uniform_mesh(0.01, 20, Val(2))

        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(50000.0), bottom = AdiabaticBC()),
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
            @test problem.bc_set.mass.top isa ImpermeableBC
            @test problem.bc_set.mass.bottom isa ImpermeableBC
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

        @testset "Generation-sink datum warning" begin
            problem = PyrolysisProblem(mesh, material, bc_set, (0.0, 10.0))

            # :with_generation_sink is consistent only with T_ref-datum
            # reaction heats; selecting it must warn about the datum
            # requirement (DSC-style heats double-count gas sensible enthalpy)
            @test_logs (:warn, r"T_ref.*datum"s) match_mode = :any Pyrolysis.Problem.to_problem_def(
                problem;
                energy_form = :with_generation_sink,
            )

            # The default :standard form stays silent
            @test_logs Pyrolysis.Problem.to_problem_def(problem)
        end
    end

    # ========================================================================
    @testset "Problem Validation" begin
        material = create_test_material()

        @testset "Valid problem passes" begin
            mesh = create_uniform_mesh(0.01, 20, Val(2))
            bc_set = BoundaryConditionSet(
                thermal = (top = HeatFluxBC(50000.0), bottom = AdiabaticBC()),
            )
            problem = PyrolysisProblem(mesh, material, bc_set, (0.0, 10.0))

            @test validate_problem(problem) == true
            @test is_valid_problem(problem) == true
        end

        @testset "BoundaryConditionSet defaults" begin
            mesh = create_uniform_mesh(0.01, 20, Val(2))
            # Only provide :top thermal BC - others should default
            bc_set = BoundaryConditionSet(thermal = (top = HeatFluxBC(50000.0),))

            problem = PyrolysisProblem(
                mesh = mesh,
                material = material,
                bc_set = bc_set,
                tspan = (0.0, 10.0),
            )

            # BoundaryConditionSet always has all BCs - missing ones filled with defaults
            @test problem.bc_set.thermal.top isa HeatFluxBC
            @test problem.bc_set.thermal.bottom isa AdiabaticBC  # Default
            @test problem.bc_set.mass.top isa ImpermeableBC  # Default
            @test problem.bc_set.mass.bottom isa ImpermeableBC  # Default
            @test validate_problem(problem) == true
            @test is_valid_problem(problem) == true
        end

        @testset "Invalid tspan detected" begin
            mesh = create_uniform_mesh(0.01, 20, Val(2))
            bc_set = BoundaryConditionSet(
                thermal = (top = HeatFluxBC(50000.0), bottom = AdiabaticBC()),
                mass = (top = ImpermeableBC(), bottom = ImpermeableBC()),
            )

            # t_start >= t_end - construct with invalid tspan
            problem = PyrolysisProblem(
                mesh = mesh,
                material = material,
                bc_set = bc_set,
                T_initial = z -> 300.0,
                ξ_initial = [z -> 1190.0, z -> 0.0],
                tspan = (10.0, 5.0),
                dt_initial = 1e-4,
                output_times = Float64[],
            )

            @test_throws ErrorException validate_problem(problem)
        end

        @testset "Darcy zero-porosity rejection" begin
            # Permeable wood with NO ρ_intrinsic: intrinsic density falls back
            # to bulk, so the solid fills the cell exactly → φ₀ = 0 while
            # κ_eff = 1e-14 > 0. Physically inconsistent and numerically fatal
            # (the φ floor in the ideal-gas pressure inversion amplifies P
            # unboundedly once any trial gas mass appears).
            function darcy_material(; ρ_intrinsic = nothing)
                kwargs = (; ρ = 500.0, c = 1500.0, k = 0.15, ε = 0.9, κ = 1e-14)
                virgin =
                    ρ_intrinsic === nothing ? SolidComponent(:wood; kwargs...) :
                    SolidComponent(:wood; kwargs..., ρ_intrinsic = ρ_intrinsic)
                gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5)
                rxn = Reaction(:pyrolysis, 1 => 2, A = 1e10, E = 150e3, h = -300e3)
                return Material(
                    name = :wood,
                    components = (virgin, gas),
                    reactions = (rxn,),
                )
            end
            mesh20() = create_uniform_mesh(0.01, 20, Val(2))
            bc_set = BoundaryConditionSet(
                thermal = (top = HeatFluxBC(50000.0), bottom = AdiabaticBC()),
            )
            make_problem(mat; kw...) = PyrolysisProblem(;
                mesh = mesh20(),
                material = mat,
                bc_set = bc_set,
                tspan = (0.0, 10.0),
                kw...,
            )

            @testset "φ₀ = 0 with κ > 0 rejected" begin
                problem = make_problem(darcy_material())
                @test_throws ErrorException validate_problem(problem)
                @test is_valid_problem(problem) == false
                # The error must teach the fix
                err = try
                    validate_problem(problem)
                    nothing
                catch e
                    sprint(showerror, e)
                end
                @test occursin("ρ_intrinsic", err)
                @test occursin("background_gas_concentration", err)
                @test occursin("porosity", err)
            end

            @testset "ρ_intrinsic restores validity" begin
                # φ₀ = 1 - 500/1500 = 2/3 — same material otherwise
                problem = make_problem(darcy_material(ρ_intrinsic = 1500.0))
                @test validate_problem(problem) == true
            end

            @testset "partial zero-porosity region rejected" begin
                # Only the lower half starts fully dense; still an error
                mat = darcy_material(ρ_intrinsic = 1500.0)
                problem = make_problem(
                    mat;
                    ξ_initial = [z -> z < 0.005 ? 1500.0 : 500.0, z -> 0.0],
                )
                @test_throws ErrorException validate_problem(problem)
            end

            @testset "impermeable nonporous material passes (Fickian)" begin
                # No κ anywhere → Fickian transport; φ₀ = 0 is irrelevant
                virgin = SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86)
                gas = GaseousComponent(:MMA, M = 0.100, c = 1500, k = 0.02, λ = 1e-5)
                rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
                mat = Material(name = :PMMA, components = (virgin, gas), reactions = (rxn,))
                @test validate_problem(make_problem(mat)) == true
            end

            @testset "κ_eff = 0 at the initial state passes" begin
                # Virgin is impermeable (κ = 0, fully dense); only the char
                # product carries κ. has_darcy_flow is true, but the initial
                # composition has zero effective permeability, so the
                # zero-porosity cells are inert and must not be flagged.
                virgin = SolidComponent(:virgin, ρ = 500.0, c = 1500.0, k = 0.15, ε = 0.9)
                char = SolidComponent(
                    :char,
                    ρ = 100.0,
                    c = 1100.0,
                    k = 0.1,
                    ε = 0.95,
                    κ = 1e-12,
                    ρ_intrinsic = 1950.0,
                )
                gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5)
                rxn = Reaction(
                    :pyrolysis,
                    1 => (2, 3),
                    A = 1e10,
                    E = 150e3,
                    h = -300e3,
                    yields = (0.3, 0.7),
                )
                mat = Material(
                    name = :charring,
                    components = (virgin, char, gas),
                    reactions = (rxn,),
                )
                mesh3 = create_uniform_mesh(0.01, 20, Val(3))
                problem = PyrolysisProblem(;
                    mesh = mesh3,
                    material = mat,
                    bc_set = bc_set,
                    tspan = (0.0, 10.0),
                )
                @test validate_problem(problem) == true
            end
        end

        @testset "degenerate Fickian gas transport warns" begin
            # The Landau ablation deck: a gas-producing reaction, Fickian
            # transport, and a solid without gas-transfer coefficient.
            function ablation_material(; λ_solid = 0.0, κ = 0.0, reactive = true)
                solid = SolidComponent(
                    :solid,
                    ρ = 1000.0,
                    c = 1000.0,
                    k = 1.0,
                    λ = λ_solid,
                    κ = κ,
                    ρ_intrinsic = 2000.0,
                )
                gas = GaseousComponent(:gas, M = 0.05, c = 1000.0, k = 0.02, λ = 1e-4)
                rxns =
                    reactive ? (Reaction(:ablate, 1 => 2; A = 1e10, E = 1e5, h = 1e6),) : ()
                return Material(
                    name = :ablation,
                    components = (solid, gas),
                    reactions = rxns,
                )
            end
            problem_for(mat) = PyrolysisProblem(
                mesh = create_uniform_mesh(0.01, 10, Val(2)),
                material = mat,
                bc_set = BoundaryConditionSet(
                    thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()),
                ),
                tspan = (0.0, 1.0),
            )
            @test_logs (:warn, r"gas-transfer coefficient") validate_problem(
                problem_for(ablation_material()),
            )
            # A condensed gas-transfer coefficient, a permeable (Darcy) matrix,
            # or no gas production keeps validation silent.
            @test_logs validate_problem(problem_for(ablation_material(λ_solid = 1e-6)))
            @test_logs validate_problem(problem_for(ablation_material(κ = 1e-14)))
            @test_logs validate_problem(problem_for(ablation_material(reactive = false)))
        end
    end

    # ========================================================================
    @testset "Initial Conditions" begin
        material = create_test_material()
        mesh = create_uniform_mesh(0.01, 10, Val(2))

        @testset "Apply from problem" begin
            bc_set = BoundaryConditionSet(
                thermal = (top = HeatFluxBC(50000.0), bottom = AdiabaticBC()),
            )
            problem = PyrolysisProblem(
                mesh = mesh,
                material = material,
                bc_set = bc_set,
                T_initial = z -> 350.0,
                ξ_initial = [z -> 1000.0, z -> 100.0],
                tspan = (0.0, 10.0),
            )

            apply_initial_conditions!(problem)

            # Check first cell
            @test mesh.cell_states[1].T ≈ 350.0
            @test mesh.cell_states[1].ξ[1] ≈ 1000.0
            @test mesh.cell_states[1].ξ[2] ≈ 100.0

            # The problem-level method also populates the χ progress
            # denominator from the freshly written cell states (solid ξ only).
            @test length(mesh.m_dry_solid_initial) == length(mesh.cells)
            @test mesh.m_dry_solid_initial[1] ≈ 1000.0 * mesh.cell_volumes[1]
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
            @test step_profile(
                0.003;
                z_interface = 0.005,
                value_left = 100.0,
                value_right = 200.0,
            ) ≈ 100.0
            @test step_profile(
                0.007;
                z_interface = 0.005,
                value_left = 100.0,
                value_right = 200.0,
            ) ≈ 200.0

            # Linear profile
            @test linear_profile(
                0.0;
                z_min = 0.0,
                z_max = 0.01,
                value_min = 300.0,
                value_max = 400.0,
            ) ≈ 300.0
            @test linear_profile(
                0.005;
                z_min = 0.0,
                z_max = 0.01,
                value_min = 300.0,
                value_max = 400.0,
            ) ≈ 350.0
            @test linear_profile(
                0.01;
                z_min = 0.0,
                z_max = 0.01,
                value_min = 300.0,
                value_max = 400.0,
            ) ≈ 400.0

            # Gaussian profile
            val = gaussian_profile(
                0.005;
                z_center = 0.005,
                width = 0.001,
                amplitude = 100.0,
                baseline = 300.0,
            )
            @test val ≈ 400.0  # At center, should be baseline + amplitude

            # Tanh profile
            val_left = tanh_profile(
                0.0;
                z_center = 0.005,
                width = 0.001,
                value_left = 100.0,
                value_right = 200.0,
            )
            val_right = tanh_profile(
                0.01;
                z_center = 0.005,
                width = 0.001,
                value_left = 100.0,
                value_right = 200.0,
            )
            @test val_left < 110.0  # Close to left value
            @test val_right > 190.0  # Close to right value
        end
    end

    # ========================================================================
    @testset "PyrolysisSolution" begin
        material = create_test_material()
        mesh = create_uniform_mesh(0.01, 10, Val(2))
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(50000.0), bottom = AdiabaticBC()),
        )
        problem = PyrolysisProblem(mesh, material, bc_set, (0.0, 10.0))

        @testset "Empty solution creation" begin
            sol = PyrolysisSolution(problem)
            @test sol.retcode == :Default
            @test isempty(sol.t)
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
        simple_material = Material(name = :Simple, components = (virgin,), reactions = ())

        mesh = create_uniform_mesh(0.01, 10, Val(1); T_initial = 300.0)
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(10000.0), bottom = AdiabaticBC()),
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
            @test length(sol.t) >= 3
            @test sol.solve_time > 0.0

            # Temperature should increase at surface
            @test sol.surface_T[end] > sol.surface_T[1]

            # Check solution structure
            @test size(sol.T, 1) > 0
            @test size(sol.ξ, 1) > 0
            @test size(sol.T, 2) == length(sol.t)
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

        material =
            Material(name = :TestMaterial, components = (virgin, gas), reactions = (rxn,))

        mesh = create_uniform_mesh(0.005, 20, Val(2))
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(100000.0), bottom = AdiabaticBC()),
            mass = (
                top = ConvectiveMassBC(h_m = 0.1, ξ_inf = 0.0),
                bottom = ImpermeableBC(),
            ),
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
        @test sol.surface_T[end] > sol.surface_T[1]

        # Check reaction is occurring: solid concentration should decrease
        initial_solid = sum(sol.ξ[:, :, 1][1, :])
        final_solid = sum(sol.ξ[:, :, end][1, :])
        @test final_solid < initial_solid  # Solid consumed by reaction

        # Gas should be generated (even if some escapes)
        # Check that gas has been produced at some point
        max_gas = maximum(sum(sol.ξ[:, :, i][2, :]) for i = 1:length(sol.t))
        @test max_gas > 0.0  # Gas was produced
    end

end
