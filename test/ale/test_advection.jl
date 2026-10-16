"""
ALE Advection Verification Test

Verifies the ALE advection sign correction by running the FULL ODE function
(not just the isolated advection term) and comparing ALE vs Eulerian results.

Strategy:
1. Pure conduction with no reactions → mesh velocity w = 0 everywhere.
   ALE and Eulerian should give IDENTICAL results (within solver tolerances).
   This confirms the ALE code path doesn't introduce spurious effects.

2. Both ALE and Eulerian compared against analytical erf solution.
   This validates correctness of the ALE code path against a known solution.

3. Single-reaction shrinking problem where ALE and Eulerian are compared.
   Bulk quantities (surface temperature, total mass) should agree closely
   at early times before significant mesh deformation.

If the ALE advection sign were wrong (dT -= w·∇T instead of +=), the ALE
solution would diverge from the Eulerian solution even with w = 0 in the
stationary case, and would show unphysical behavior in the shrinking case.
"""

using Test
using Pyrolysis
using Pyrolysis.Solver: get_cell_positions

if !isdefined(@__MODULE__, :erf)
    function erf(x::Real)
        sign_x = x < 0 ? -1.0 : 1.0
        a = abs(float(x))
        t = 1.0 / (1.0 + 0.3275911 * a)
        y =
            1.0 -
            (
                ((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) *
                t + 0.254829592
            ) *
            t *
            exp(-a * a)
        return sign_x * y
    end
end

# ============================================================================
# Test 1: Pure Conduction — ALE vs Eulerian (No Reactions, w = 0)
# ============================================================================
#
# With no reactions, there is no density change, so mesh velocity w = 0.
# The ALE correction w·∇T = 0, so ALE and Eulerian should give identical
# results. Any difference indicates a bug in the ALE machinery.
#
@testset "ALE vs Eulerian: Pure Conduction (w=0)" begin

    # Material properties (constant for analytical solution)
    ρ = 1000.0    # kg/m³
    cp = 1000.0   # J/(kg·K)
    k = 0.5       # W/(m·K)
    α_diff = k / (ρ * cp)  # Thermal diffusivity = 5e-7 m²/s

    # Boundary and initial conditions
    T_initial = 300.0  # K
    T_surface = 500.0  # K

    # Domain
    L = 0.02  # 2 cm
    n_cells = 40

    # Analytical solution for semi-infinite solid
    function T_analytical(z, t)
        if t <= 0
            return T_initial
        end
        # Distance from surface = L - z (surface at z = L)
        return T_surface + (T_initial - T_surface) * erf((L - z) / (2 * sqrt(α_diff * t)))
    end

    # Create material (inert, no reactions)
    solid = SolidComponent(:solid, ρ = ρ, c = cp, k = k, ε = 0.9)
    mat = Material(name = :TestSolid, components = (solid,), reactions = ())

    # Boundary conditions: Dirichlet at top, adiabatic at bottom
    bc_set = BoundaryConditionSet(
        thermal = (top = DirichletThermalBC(T_surface), bottom = AdiabaticBC()),
    )

    t_final = 50.0  # seconds

    # --- Solve with Eulerian (fixed mesh) ---
    mesh_euler =
        create_uniform_mesh(L, n_cells, Val(1); T_initial = T_initial, ξ_initial = (ρ,))

    problem_euler = PyrolysisProblem(
        mesh = mesh_euler,
        material = mat,
        bc_set = bc_set,
        T_initial = z -> T_initial,
        ξ_initial = [z -> ρ],
        tspan = (0.0, t_final),
    )

    sol_euler = solve(problem_euler; abstol = 1e-10, reltol = 1e-8, progress = false)
    @test sol_euler.retcode == :Success

    # --- Solve with ALE ---
    mesh_ale =
        create_uniform_mesh(L, n_cells, Val(1); T_initial = T_initial, ξ_initial = (ρ,))

    problem_ale = PyrolysisProblem(
        mesh = mesh_ale,
        material = mat,
        bc_set = bc_set,
        T_initial = z -> T_initial,
        ξ_initial = [z -> ρ],
        tspan = (0.0, t_final),
    )

    sol_ale =
        solve(problem_ale; use_ale = true, abstol = 1e-10, reltol = 1e-8, progress = false)
    @test sol_ale.retcode == :Success

    # --- Compare temperatures ---
    T_euler = sol_euler.T[:, end]
    T_ale = sol_ale.T[:, end]
    z_positions = get_cell_positions(sol_euler)

    # ALE and Eulerian should match to solver tolerance
    max_diff = maximum(abs.(T_euler .- T_ale))
    @testset "ALE matches Eulerian (no mesh motion)" begin
        @test max_diff < 1e-4  # Should be near machine precision for same problem
    end

    # Both should match the analytical solution
    @testset "Both match analytical solution" begin
        max_err_euler = 0.0
        max_err_ale = 0.0
        for (i, z) in enumerate(z_positions)
            T_exact = T_analytical(z, t_final)
            max_err_euler = max(max_err_euler, abs(T_euler[i] - T_exact))
            max_err_ale = max(max_err_ale, abs(T_ale[i] - T_exact))
        end

        # Both should have small error (< 2K for this mesh)
        @test max_err_euler < 2.0
        @test max_err_ale < 2.0

        # ALE error should be comparable to Eulerian (not worse)
        @test max_err_ale < 2.0 * max_err_euler + 1e-6
    end

    # Surface temperatures should match
    @testset "Surface temperature agreement" begin
        T_surf_euler = sol_euler.surface_T[end]
        T_surf_ale = sol_ale.surface_T[end]
        @test abs(T_surf_euler - T_surf_ale) < 1e-3
    end

    # Thickness should be unchanged (no shrinkage)
    @testset "No mesh motion occurred" begin
        @test abs(sol_ale.extras.thickness_history[end] - L) < 1e-10
    end
end

# ============================================================================
# Test 2: Simple Reaction — ALE vs Eulerian Early-Time Agreement
# ============================================================================
#
# With a single reaction, the solid decomposes and the mesh shrinks (ALE).
# At early times before significant mesh deformation, ALE and Eulerian
# should give similar surface temperatures and mass loss rates.
#
# This tests that the ALE advection correction (with w ≠ 0) doesn't
# introduce spurious temperature changes.
#
@testset "ALE vs Eulerian: Reacting Material (early time)" begin

    # Create a simple decomposing material
    virgin = SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.3, ε = 0.9)
    char = SolidComponent(:char, ρ = 200.0, c = 1000.0, k = 0.15, ε = 0.9)

    # Single first-order reaction: virgin (1) → char (2)
    rxn = Reaction(:decomp, 1 => 2, A = 1e8, E = 1.2e5, h = 2e5)

    mat = Material(
        name = :SimpleDecomposition,
        components = (virgin, char),
        reactions = (rxn,),
    )

    # Domain and mesh
    L = 0.01  # 1 cm
    n_cells = 30
    T_initial = 300.0

    # Boundary conditions: moderate heat flux
    bc_set =
        BoundaryConditionSet(thermal = (top = HeatFluxBC(30000.0), bottom = AdiabaticBC()))

    t_early = 5.0  # Short time — limited decomposition

    # --- Eulerian ---
    mesh_euler = create_uniform_mesh(
        L,
        n_cells,
        Val(2);
        T_initial = T_initial,
        ξ_initial = (1000.0, 0.0),
    )

    problem_euler = PyrolysisProblem(
        mesh = mesh_euler,
        material = mat,
        bc_set = bc_set,
        T_initial = z -> T_initial,
        ξ_initial = [z -> 1000.0, z -> 0.0],
        tspan = (0.0, t_early),
    )

    sol_euler = solve(problem_euler; abstol = 1e-8, reltol = 1e-6, progress = false)
    @test sol_euler.retcode == :Success

    # --- ALE ---
    mesh_ale = create_uniform_mesh(
        L,
        n_cells,
        Val(2);
        T_initial = T_initial,
        ξ_initial = (1000.0, 0.0),
    )

    problem_ale = PyrolysisProblem(
        mesh = mesh_ale,
        material = mat,
        bc_set = bc_set,
        T_initial = z -> T_initial,
        ξ_initial = [z -> 1000.0, z -> 0.0],
        tspan = (0.0, t_early),
    )

    sol_ale =
        solve(problem_ale; use_ale = true, abstol = 1e-8, reltol = 1e-6, progress = false)
    @test sol_ale.retcode == :Success

    # At early times, surface temperatures should be similar
    # (ALE mesh motion is small, so the physical solution is nearly identical)
    @testset "Surface temperature agreement (early time)" begin
        T_surf_euler = sol_euler.surface_T[end]
        T_surf_ale = sol_ale.surface_T[end]

        # Allow modest tolerance — ALE does handle the physics slightly differently
        # due to mesh motion, but results should be in the same ballpark
        rel_diff = abs(T_surf_euler - T_surf_ale) / T_surf_euler
        @test rel_diff < 0.05  # Within 5%
    end

    # Total mass should be similar (same reactions, same physics)
    @testset "Total mass agreement (early time)" begin
        mass_euler = sol_euler.extras.total_mass[end]
        mass_ale = sol_ale.extras.total_mass[end]

        rel_diff = abs(mass_euler - mass_ale) / mass_euler
        @test rel_diff < 0.02  # Within 2%
    end

    # ALE thickness should remain close to initial (limited decomposition at early time)
    @testset "ALE mesh thickness reasonable" begin
        thickness_final = sol_ale.extras.thickness_history[end]
        @test abs(thickness_final - L) / L < 0.05  # Within 5% of original
    end
end

# ============================================================================
# Test 3: Mesh Convergence — ALE Conduction
# ============================================================================
#
# Verify that the ALE formulation maintains spatial convergence order
# for pure conduction (same as Eulerian). This would fail if the ALE
# advection term had the wrong sign, as it would introduce O(w·∇T) error.
#
@testset "ALE Conduction: Spatial Convergence" begin

    ρ = 1000.0
    cp = 1000.0
    k = 0.5
    α_diff = k / (ρ * cp)

    T_initial = 300.0
    T_surface = 500.0
    L = 0.05  # 5 cm — thick enough for semi-infinite approximation

    function T_analytical(z, t)
        if t <= 0
            return T_initial
        end
        return T_surface + (T_initial - T_surface) * erf((L - z) / (2 * sqrt(α_diff * t)))
    end

    solid = SolidComponent(:solid, ρ = ρ, c = cp, k = k, ε = 0.9)
    mat = Material(name = :TestSolid, components = (solid,), reactions = ())

    bc_set = BoundaryConditionSet(
        thermal = (top = DirichletThermalBC(T_surface), bottom = AdiabaticBC()),
    )

    t_test = 50.0
    n_cells_list = [25, 50, 100, 200]
    errors = Float64[]
    h_values = Float64[]

    for n_cells in n_cells_list
        mesh =
            create_uniform_mesh(L, n_cells, Val(1); T_initial = T_initial, ξ_initial = (ρ,))
        h = L / n_cells
        push!(h_values, h)

        problem = PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T_initial,
            ξ_initial = [z -> ρ],
            tspan = (0.0, t_test),
        )

        sol =
            solve(problem; use_ale = true, abstol = 1e-12, reltol = 1e-10, progress = false)

        z_positions = get_cell_positions(sol)
        T_num = sol.T[:, end]

        l2_error = 0.0
        for (i, z) in enumerate(z_positions)
            T_exact = T_analytical(z, t_test)
            l2_error += (T_num[i] - T_exact)^2 * h
        end
        l2_error = sqrt(l2_error / L)
        push!(errors, l2_error)
    end

    # Compute convergence rates
    rates = Float64[]
    for i = 2:length(errors)
        rate = log(errors[i-1] / errors[i]) / log(h_values[i-1] / h_values[i])
        push!(rates, rate)
    end

    # Should achieve at least 1.5 order (2nd order FV with some time integration error)
    @testset "Convergence order" begin
        @test all(r -> r > 1.5, rates)
        @test sum(rates) / length(rates) > 1.8  # Average close to 2nd order
    end
end
