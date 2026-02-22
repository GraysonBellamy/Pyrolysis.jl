"""
Tests for radiation model selection and implementation.

This test suite validates:
- Radiation model enum and helper functions
- Surface absorption physics and energy conservation
- Beer-Lambert physics and energy conservation
- Warning system for model misconfigurations
- Jacobian sparsity patterns for each model
- End-to-end solver integration

Reference: docs/Radiation_Model_Refactor_Plan.md
"""

using Test
using Pyrolysis
using SparseArrays

@testset "Radiation Model Enum and Helpers" begin

    @testset "Enum Values" begin
        @test Int(NO_RADIATION) == 0
        @test Int(SURFACE_ABSORPTION) == 1
        @test Int(BEER_LAMBERT) == 2
        @test Int(P1_QUASI_STEADY) == 3
    end

    @testset "requires_absorption_coefficients" begin
        @test requires_absorption_coefficients(NO_RADIATION) == false
        @test requires_absorption_coefficients(SURFACE_ABSORPTION) == false
        @test requires_absorption_coefficients(BEER_LAMBERT) == true
        @test requires_absorption_coefficients(P1_QUASI_STEADY) == true
    end

    @testset "has_nonlocal_jacobian_coupling" begin
        @test has_nonlocal_jacobian_coupling(NO_RADIATION) == false
        @test has_nonlocal_jacobian_coupling(SURFACE_ABSORPTION) == false
        @test has_nonlocal_jacobian_coupling(BEER_LAMBERT) == true
        @test has_nonlocal_jacobian_coupling(P1_QUASI_STEADY) == false
    end

    @testset "requires_radiation_cache" begin
        @test requires_radiation_cache(NO_RADIATION) == false
        @test requires_radiation_cache(SURFACE_ABSORPTION) == false
        @test requires_radiation_cache(BEER_LAMBERT) == false
        @test requires_radiation_cache(P1_QUASI_STEADY) == true
    end
end

@testset "Surface Absorption Model" begin

    @testset "Energy Conservation" begin
        # Create test mesh
        mesh = create_uniform_mesh(0.01, 10, Val(2); T_initial=300.0, ξ_initial=(1200.0, 0.0))
        q_rad = zeros(10)
        I_surface = 50e3  # W/m²

        # Compute radiation source
        Pyrolysis.Physics.compute_surface_radiation_source!(q_rad, mesh, I_surface)

        # Only surface cell (cell 10, at z=L) should have non-zero source
        @test all(q_rad[1:9] .== 0.0)
        @test q_rad[10] > 0.0

        # Total absorbed should equal incident flux
        Δz_surface = Pyrolysis.Geometry.cell_size(mesh, 10)
        total_absorbed = q_rad[10] * Δz_surface
        @test total_absorbed ≈ I_surface rtol=1e-10
    end

    @testset "Warning for Ignored α" begin
        # Create material with absorption coefficients
        mat = Material(
            name = :TestMat,
            components = (
                SolidComponent(:solid, ρ=1000.0, c=1000.0, k=0.2, ε=0.9, α=1.5),
                GaseousComponent(:gas, M=0.044, c=1000.0, k=0.02, λ=1e-5)
            ),
            reactions = (
                Reaction(:decomp, 1 => 2, A=1e10, E=150e3, h=500e3, n=1.0),
            )
        )

        @test has_absorption_coefficients(mat.components) == true

        # Warning should be issued when using SURFACE_ABSORPTION with α > 0
        @test_logs (:warn, r"SURFACE_ABSORPTION"s) begin
            Pyrolysis.Solvers._validate_radiation_model(SURFACE_ABSORPTION, mat, nothing)
        end
    end
end

@testset "Beer-Lambert Model" begin

    @testset "In-Depth Penetration" begin
        mesh = create_uniform_mesh(0.01, 10, Val(2); T_initial=300.0, ξ_initial=(1200.0, 0.0))
        prop_cache = Pyrolysis.Properties.PropertyCache(10, 11)

        # Set uniform absorption coefficient
        α_eff = 100.0  # m⁻¹
        fill!(prop_cache.α, α_eff)

        q_rad = zeros(10)
        I_surface = 50e3  # W/m²

        # Compute Beer-Lambert radiation source
        Pyrolysis.Physics.compute_beer_lambert_source!(q_rad, mesh, prop_cache, I_surface)

        # All cells should have some absorption
        @test all(q_rad .> 0.0)

        # Surface cell (cell 10) should have highest absorption
        @test q_rad[10] == maximum(q_rad)

        # Absorption should decrease monotonically with depth
        @test issorted(q_rad)  # Cells 1-10 are ordered bottom to top

        # Energy conservation: absorbed + transmitted = incident
        total_absorbed = sum(q_rad[i] * Pyrolysis.Geometry.cell_size(mesh, i) for i in 1:10)
        transmitted = Pyrolysis.Physics.transmitted_radiation(I_surface, mesh, prop_cache)
        @test total_absorbed + transmitted ≈ I_surface rtol=1e-10
    end

    @testset "Warning for Missing α" begin
        # Create material without absorption coefficients
        mat = Material(
            name = :TestMat,
            components = (
                SolidComponent(:solid, ρ=1000.0, c=1000.0, k=0.2, ε=0.9, α=0.0),
                GaseousComponent(:gas, M=0.044, c=1000.0, k=0.02, λ=1e-5)
            ),
            reactions = (
                Reaction(:decomp, 1 => 2, A=1e10, E=150e3, h=500e3, n=1.0),
            )
        )

        @test has_absorption_coefficients(mat.components) == false

        # Warning should be issued when using BEER_LAMBERT without α
        @test_logs (:warn, r"BEER_LAMBERT"s) begin
            Pyrolysis.Solvers._validate_radiation_model(BEER_LAMBERT, mat, nothing)
        end
    end
end

@testset "Jacobian Sparsity Patterns" begin

    n_cells = 10
    NC = 2

    @testset "Tridiagonal Models" begin
        # NO_RADIATION, SURFACE_ABSORPTION, P1_QUASI_STEADY should all have same sparsity
        J_none = Pyrolysis.TimeIntegration.jacobian_sparsity(n_cells, NC; radiation_model=NO_RADIATION)
        J_surface = Pyrolysis.TimeIntegration.jacobian_sparsity(n_cells, NC; radiation_model=SURFACE_ABSORPTION)
        J_p1 = Pyrolysis.TimeIntegration.jacobian_sparsity(n_cells, NC; radiation_model=P1_QUASI_STEADY)

        nnz_none = nnz(J_none)
        nnz_surface = nnz(J_surface)
        nnz_p1 = nnz(J_p1)

        @test nnz_surface == nnz_none
        @test nnz_p1 == nnz_none

        # All entries should match exactly
        @test J_surface.rowval == J_none.rowval
        @test J_surface.colptr == J_none.colptr
        @test J_p1.rowval == J_none.rowval
        @test J_p1.colptr == J_none.colptr
    end

    @testset "BEER_LAMBERT Non-Local Coupling" begin
        J_none = Pyrolysis.TimeIntegration.jacobian_sparsity(n_cells, NC; radiation_model=NO_RADIATION)
        J_beer = Pyrolysis.TimeIntegration.jacobian_sparsity(n_cells, NC; radiation_model=BEER_LAMBERT)

        nnz_none = nnz(J_none)
        nnz_beer = nnz(J_beer)

        # BEER_LAMBERT should have MORE entries due to non-local coupling
        @test nnz_beer > nnz_none

        # For 10 cells, 2 components, expect roughly 1.5x more entries
        ratio = nnz_beer / nnz_none
        @test 1.4 < ratio < 1.6
    end

    @testset "ALE Sparsity Patterns" begin
        # ALE versions should also follow same pattern
        J_ale_none = Pyrolysis.TimeIntegration.jacobian_sparsity_ale(n_cells, NC; radiation_model=NO_RADIATION)
        J_ale_surface = Pyrolysis.TimeIntegration.jacobian_sparsity_ale(n_cells, NC; radiation_model=SURFACE_ABSORPTION)
        J_ale_beer = Pyrolysis.TimeIntegration.jacobian_sparsity_ale(n_cells, NC; radiation_model=BEER_LAMBERT)
        J_ale_p1 = Pyrolysis.TimeIntegration.jacobian_sparsity_ale(n_cells, NC; radiation_model=P1_QUASI_STEADY)

        nnz_ale_none = nnz(J_ale_none)

        # Tridiagonal models should match
        @test nnz(J_ale_surface) == nnz_ale_none
        @test nnz(J_ale_p1) == nnz_ale_none

        # BEER_LAMBERT should have more entries
        @test nnz(J_ale_beer) > nnz_ale_none
    end
end

@testset "End-to-End Solver Integration" begin

    # Create simple test material
    material = Material(
        name = :test_material,
        components = (
            SolidComponent(:solid, ρ=1200.0, c=1500.0, k=0.2, ε=0.9, α=10.0),
            GaseousComponent(:gas, M=0.044, c=1000.0, k=0.02, λ=1e-5)
        ),
        reactions = (
            Reaction(:decomposition, 1 => 2, A=1e10, E=150e3, h=500e3, n=1.0),
        )
    )

    bc_set = BoundaryConditionSet(
        thermal_top = HeatFluxBC(50e3),
        thermal_bottom = AdiabaticBC(),
        mass_top = DiffusiveMassBC(Inf),
        mass_bottom = ImpermeableBC()
    )

    @testset "Fixed Mesh - All Models" begin
        for model in [NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT]
            mesh = create_uniform_mesh(0.01, 5, Val(2); T_initial=300.0, ξ_initial=(1200.0, 0.0))
            problem = PyrolysisProblem(mesh, material, bc_set, (0.0, 0.5))

            sol = solve(problem;
                radiation_model=model,
                use_ale=false,
                saveat=[0.0, 0.5],
                progress=false)

            @test sol.retcode == :Success
        end
    end

    @testset "ALE - All Models" begin
        for model in [NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT]
            mesh = create_uniform_mesh(0.01, 5, Val(2); T_initial=300.0, ξ_initial=(1200.0, 0.0))
            problem = PyrolysisProblem(mesh, material, bc_set, (0.0, 0.5))

            sol = solve(problem;
                radiation_model=model,
                use_ale=true,
                saveat=[0.0, 0.5],
                progress=false)

            @test sol.retcode == :Success
        end
    end

    @testset "With Sparse Jacobian" begin
        for model in [SURFACE_ABSORPTION, BEER_LAMBERT]
            mesh = create_uniform_mesh(0.01, 5, Val(2); T_initial=300.0, ξ_initial=(1200.0, 0.0))
            problem = PyrolysisProblem(mesh, material, bc_set, (0.0, 0.5))

            sol = solve(problem;
                radiation_model=model,
                use_ale=false,
                use_sparse_jacobian=true,
                saveat=[0.0, 0.5],
                progress=false)

            @test sol.retcode == :Success
        end
    end
end

@testset "Model Dispatch Correctness" begin

    mesh = create_uniform_mesh(0.01, 10, Val(2); T_initial=300.0, ξ_initial=(1200.0, 0.0))
    prop_cache = Pyrolysis.Properties.PropertyCache(10, 11)
    fill!(prop_cache.α, 100.0)

    q_rad = zeros(10)
    I_surface = 50e3

    @testset "NO_RADIATION dispatch" begin
        Pyrolysis.Physics.compute_radiation_source!(q_rad, mesh, prop_cache, I_surface, NO_RADIATION)
        @test all(q_rad .== 0.0)
    end

    @testset "SURFACE_ABSORPTION dispatch" begin
        Pyrolysis.Physics.compute_radiation_source!(q_rad, mesh, prop_cache, I_surface, SURFACE_ABSORPTION)
        @test sum(q_rad[1:9]) == 0.0  # Only surface cell
        @test q_rad[10] > 0.0
    end

    @testset "BEER_LAMBERT dispatch" begin
        Pyrolysis.Physics.compute_radiation_source!(q_rad, mesh, prop_cache, I_surface, BEER_LAMBERT)
        @test all(q_rad .> 0.0)  # All cells have absorption
        @test issorted(q_rad)  # Monotonically increasing toward surface
    end

    @testset "P1_QUASI_STEADY dispatch" begin
        # P1 is handled separately via radiation cache, so this should zero q_rad
        Pyrolysis.Physics.compute_radiation_source!(q_rad, mesh, prop_cache, I_surface, P1_QUASI_STEADY)
        @test all(q_rad .== 0.0)
    end
end
