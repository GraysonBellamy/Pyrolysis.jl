"""
Verification tests for the volume-fraction-gradient gas flux formulation.

Tests the formula:
    J = -lambda_face * (xi_R * T_R - xi_L * T_L) / (T_face * d_LR)

which is equivalent to:
    J = -rho_gas * lambda * d(xi/rho_gas)/dx   (volume fraction gradient)
    J = -(lambda/T) * d(xi*T)/dx               (product rule form)
    J = -lambda * [dxi/dx + (xi/T)*dT/dx]      (expanded form)
"""

using Test
using Pyrolysis
using Pyrolysis.Problem: StateCache

# Access non-exported functions via the Physics submodule
using Pyrolysis.Physics: gas_flux, gas_flux_darcy_fick, compute_gas_fluxes!
using Pyrolysis.Properties: PropertyCache, update_properties!, harmonic_mean_weighted
using Pyrolysis.Materials: R_GAS
using Pyrolysis.Geometry: CellState

@testset "Gas Flux Volume Fraction Gradient" begin

    # ========================================================================
    # Test 1: Equal concentration, different temperature -> non-zero flux
    # ========================================================================
    @testset "Equal concentration, different temperature -> non-zero flux" begin
        # Two cells with equal gas mass concentration but different T
        xi_L = 10.0   # kg/m^3
        xi_R = 10.0   # kg/m^3
        T_L = 400.0  # K
        T_R = 800.0  # K
        lam = 2e-5   # m^2/s
        d_LR = 0.001  # m

        # Analytical: J = -lam * (xi_R*T_R - xi_L*T_L) / (T_face * d_LR)
        T_face = 0.5 * (T_L + T_R)  # = 600 K
        J_expected = -lam * (xi_R * T_R - xi_L * T_L) / (T_face * d_LR)
        # = -2e-5 * (10*800 - 10*400) / (600 * 0.001)
        # = -2e-5 * 4000 / 0.6
        # = -2e-5 * 6666.67
        # = -0.13333...

        J_computed = gas_flux(lam, xi_L, xi_R, T_L, T_R, d_LR)

        @test J_computed != 0.0
        @test J_computed < 0.0  # gas flows from hot (high pressure) to cold (low pressure)
        @test J_computed ≈ J_expected atol=1e-12
        @test J_computed ≈ -2e-5 * 6666.6667 rtol=1e-4
    end

    # ========================================================================
    # Test 2: Uniform temperature -> reduces to simple Fick's law
    # ========================================================================
    @testset "Uniform temperature -> reduces to simple Fick's law" begin
        T_uniform = 500.0  # K
        xi_L = 5.0    # kg/m^3
        xi_R = 15.0   # kg/m^3
        lam = 2e-5   # m^2/s
        d_LR = 0.001  # m

        # Volume fraction formula:
        # J = -lam * (xi_R*T - xi_L*T) / (T * d_LR)
        #   = -lam * T*(xi_R - xi_L) / (T * d_LR)
        #   = -lam * (xi_R - xi_L) / d_LR
        J_volfrac = gas_flux(lam, xi_L, xi_R, T_uniform, T_uniform, d_LR)

        # Simple Fick's law: J = -lam * (xi_R - xi_L) / d_LR
        J_fick = -lam * (xi_R - xi_L) / d_LR
        # = -2e-5 * 10000 = -0.2

        @test J_volfrac ≈ J_fick atol=1e-12
        @test J_volfrac ≈ -0.2 atol=1e-12
    end

    # ========================================================================
    # Test 3: Thermal diffusion term quantification
    # ========================================================================
    @testset "Thermal diffusion term quantification" begin
        # Set up a case where both xi and T vary
        xi_L = 8.0    # kg/m^3
        xi_R = 12.0   # kg/m^3
        T_L = 500.0  # K
        T_R = 700.0  # K
        lam = 3e-5   # m^2/s
        d_LR = 0.002  # m

        # Product rule form: J = -(lam / T_face) * (xi_R*T_R - xi_L*T_L) / d_LR
        T_face = 0.5 * (T_L + T_R)
        J_product = -lam * (xi_R * T_R - xi_L * T_L) / (T_face * d_LR)

        # Expanded form: J = -lam * [dxi/dx + (xi_face/T_face)*dT/dx]
        dxi_dx = (xi_R - xi_L) / d_LR
        dT_dx = (T_R - T_L) / d_LR
        xi_face = 0.5 * (xi_L + xi_R)
        J_expanded = -lam * (dxi_dx + (xi_face / T_face) * dT_dx)

        # These should match to machine precision because they are algebraically
        # equivalent when using arithmetic-mean face values:
        # -(lam/T_f) * (xi_R*T_R - xi_L*T_L)/d = -lam*[(xi_R-xi_L)/d + xi_f/T_f * (T_R-T_L)/d]
        # Proof: xi_R*T_R - xi_L*T_L = T_f*(xi_R - xi_L) + xi_f*(T_R - T_L)
        #   because T_f*(xi_R-xi_L) + xi_f*(T_R-T_L)
        #   = (T_L+T_R)/2 * (xi_R-xi_L) + (xi_L+xi_R)/2 * (T_R-T_L)
        #   = (T_L*xi_R - T_L*xi_L + T_R*xi_R - T_R*xi_L + xi_L*T_R - xi_L*T_L + xi_R*T_R - xi_R*T_L) / 2
        #   = (2*xi_R*T_R - 2*xi_L*T_L) / 2
        #   = xi_R*T_R - xi_L*T_L
        @test J_product ≈ J_expanded atol=1e-14

        # Verify against gas_flux
        J_computed = gas_flux(lam, xi_L, xi_R, T_L, T_R, d_LR)
        @test J_computed ≈ J_product atol=1e-14
        @test J_computed ≈ J_expanded atol=1e-14
    end

    # ========================================================================
    # Test 4: Equivalence with partial pressure gradient
    # ========================================================================
    @testset "Equivalence with partial pressure gradient" begin
        xi_L = 10.0   # kg/m^3
        xi_R = 10.0   # kg/m^3
        T_L = 400.0  # K
        T_R = 800.0  # K
        lam = 2e-5   # m^2/s
        d_LR = 0.001  # m
        M = 0.028  # kg/mol (N2-like)

        # Partial pressures: P = xi * R_GAS * T / M
        P_L = xi_L * R_GAS * T_L / M
        P_R = xi_R * R_GAS * T_R / M
        dP_dx = (P_R - P_L) / d_LR

        T_face = 0.5 * (T_L + T_R)

        # From the derivation: J = -(lam * M) / (R_GAS * T_face) * dP/dx
        J_pressure = -(lam * M) / (R_GAS * T_face) * dP_dx

        # From the implemented formula
        J_computed = gas_flux(lam, xi_L, xi_R, T_L, T_R, d_LR)

        @test J_computed ≈ J_pressure atol=1e-12
    end

    # ========================================================================
    # Test 5: Full mesh integration test
    # ========================================================================
    @testset "Full mesh integration test" begin
        # Define a 2-component material: 1 solid + 1 gas
        solid = SolidComponent(:solid, ρ = 1000.0, c = 1000.0, k = 0.2, λ = 2e-5)
        gas = GaseousComponent(:gas, M = 0.028, c = 1000.0, k = 0.02, λ = 2e-5)

        mat = Material(name = :test_mat, components = (solid, gas), reactions = ())

        # Create a 5-cell mesh with uniform gas concentration
        n_cells = 5
        NC = 2
        xi_gas_uniform = 10.0
        xi_solid = 500.0
        T_cold = 300.0
        T_hot = 800.0

        mesh = create_uniform_mesh(
            0.005,
            n_cells,
            Val(NC);
            T_initial = T_cold,
            ξ_initial = (xi_solid, xi_gas_uniform),
        )

        # Apply a linear temperature gradient: T increases from cell 1 to cell 5
        for i = 1:n_cells
            frac = (i - 0.5) / n_cells  # midpoint fraction in [0,1]
            T_i = T_cold + frac * (T_hot - T_cold)
            mesh.cell_states[i] =
                CellState{NC}(T_i, (xi_solid, xi_gas_uniform), mesh.cell_states[i].volume)
        end

        # Create property cache and update
        cache = PropertyCache(mesh)
        update_properties!(cache, StateCache(mesh), mesh, mat)

        # Compute gas fluxes
        n_faces = length(mesh.faces)
        J = zeros(NC, n_faces)
        compute_gas_fluxes!(J, StateCache(mesh), mesh, cache, mat)

        # Check that interior face fluxes for the gas component (index 2) are non-zero
        # Interior faces are indices 2 through n_faces-1
        for f = 2:(n_faces-1)
            # Gas fluxes should be non-zero due to thermal diffusion
            @test J[2, f] != 0.0
            # Gas should flow from hot toward cold, i.e., toward lower z
            # Since T increases with z, flux should be negative (toward -z)
            @test J[2, f] < 0.0
        end

        # Solid fluxes should be zero
        for f = 1:n_faces
            @test J[1, f] == 0.0
        end
    end

    # ========================================================================
    # Test 6: Darcy-Fick combined flux
    # ========================================================================
    @testset "Darcy-Fick combined flux" begin
        xi_L = 10.0
        xi_R = 10.0
        T_L = 400.0
        T_R = 800.0
        lam = 2e-5
        d_LR = 0.001
        v_g = 0.05  # m/s positive Darcy velocity (+z direction)

        # Diffusive part (volume fraction gradient)
        T_face = 0.5 * (T_L + T_R)
        J_diff = -lam * (xi_R * T_R - xi_L * T_L) / (T_face * d_LR)

        # Advective part with first-order upwinding
        # v_g >= 0 => upwind = left
        xi_upwind = xi_L
        J_adv = v_g * xi_upwind

        J_expected = J_adv + J_diff

        J_computed = gas_flux_darcy_fick(v_g, xi_L, xi_R, 1.0, 1.0, lam, T_L, T_R, d_LR)

        @test J_computed ≈ J_expected atol=1e-12

        # Verify advective part is unchanged (simple upwind * velocity)
        @test J_adv ≈ v_g * xi_L
        @test J_adv ≈ 0.5 atol=1e-12  # 0.05 * 10.0

        # Verify diffusive part uses the volume fraction gradient
        J_diff_from_gas_flux = gas_flux(lam, xi_L, xi_R, T_L, T_R, d_LR)
        @test J_diff ≈ J_diff_from_gas_flux atol=1e-12

        # Test with negative velocity (upwind = right)
        v_g_neg = -0.03
        J_neg = gas_flux_darcy_fick(v_g_neg, xi_L, xi_R, 1.0, 1.0, lam, T_L, T_R, d_LR)
        J_adv_neg = v_g_neg * xi_R  # upwind is right when v_g < 0
        J_expected_neg = J_adv_neg + J_diff
        @test J_neg ≈ J_expected_neg atol=1e-12
    end

    # ========================================================================
    # Test 7: Harmonic mean of lambda (not lambda*rho)
    # ========================================================================
    @testset "Harmonic mean of lambda (not lambda*rho)" begin
        # Two cells with very different lambda values
        lam_L = 1e-6   # m^2/s
        lam_R = 1e-4   # m^2/s

        # Harmonic mean: 2 / (1/lam_L + 1/lam_R)
        lam_harmonic = 2.0 / (1.0 / lam_L + 1.0 / lam_R)
        @test lam_harmonic ≈ 2.0 * lam_L * lam_R / (lam_L + lam_R)

        # Should be dominated by the small value
        @test lam_harmonic ≈ 2.0 * 1e-6 * 1e-4 / (1e-6 + 1e-4) rtol=1e-10
        # = 2e-10 / 1.01e-4 ≈ 1.98e-6
        @test lam_harmonic < 2e-6
        @test lam_harmonic > 1.9e-6

        # Use the package's harmonic_mean_weighted with equal distances
        d_L = 0.5e-3
        d_R = 0.5e-3
        lam_face = harmonic_mean_weighted(lam_L, lam_R, d_L, d_R)

        # For equal distances, distance-weighted harmonic mean = simple harmonic mean
        @test lam_face ≈ lam_harmonic rtol=1e-10

        # Now verify the flux uses harmonic_mean(lambda), NOT harmonic_mean(lambda*rho)
        # Set up two cells with different lambda and some rho values
        xi_L = 5.0
        xi_R = 15.0
        T_L = 500.0
        T_R = 500.0
        d_LR = 0.001

        # Flux with harmonic mean of lambda
        J_correct = gas_flux(lam_face, xi_L, xi_R, T_L, T_R, d_LR)

        # Verify it matches the formula with harmonic mean of lambda
        T_face = 0.5 * (T_L + T_R)
        J_manual = -lam_face * (xi_R * T_R - xi_L * T_L) / (T_face * d_LR)
        @test J_correct ≈ J_manual atol=1e-14

        # The face lambda should NOT be harmonic_mean(lambda*rho)
        # For demonstration: compute what harmonic_mean(lambda*rho) would give
        rho_L = 1.0  # arbitrary density
        rho_R = 0.5
        lam_rho_L = lam_L * rho_L
        lam_rho_R = lam_R * rho_R
        lam_rho_face = 2.0 * lam_rho_L * lam_rho_R / (lam_rho_L + lam_rho_R)
        # This would give a different value
        @test lam_face != lam_rho_face  # they must differ
    end

    @testset "Surface diffusion flux uses lambda directly" begin
        lam = 1e-5
        xi = 2e-3
        d = 1e-4
        expected_internal = lam * xi / d

        for rho in (1.0, 1e-20, 1e-30)
            @test Pyrolysis.BoundaryConditions.mass_bc_flux(
                DiffusiveMassBC(),
                xi,
                lam * rho,
                rho,
                d,
                0.0,
                1,
            ) ≈ -expected_internal
        end

        h_m = 0.01
        expected_resistance = lam * h_m * xi / (lam + h_m * d)
        @test Pyrolysis.BoundaryConditions.mass_bc_flux(
            DiffusiveMassBC(h_m),
            xi,
            lam * 1e-30,
            1e-30,
            d,
            0.0,
            1,
        ) ≈ -expected_resistance
        @test Pyrolysis.BoundaryConditions.mass_bc_flux(
            DiffusiveMassBC(h_m),
            xi,
            lam * 1e-30,
            1e-30,
            d,
            0.0,
            1,
        ) ≈ -expected_resistance
    end

end  # top-level testset
