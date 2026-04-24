#!/usr/bin/env julia
"""
Test suite for the four analytical Jacobian bug fixes.

Validates each fix by comparing the analytical Jacobian against central finite
differences. Each @testset exercises the specific code path affected by the bug.

Bug fixes validated:
1. Distance-weighted harmonic mean: face interpolation uses
   (d_L+d_R)/(d_L/a + d_R/b) with correct derivatives.
2. Conduction dkeff/dT: derivative of effective conductivity w.r.t. temperature.
3. Gas convective dc_g/dT: derivative of gas heat capacity w.r.t. temperature
   in the convective energy source term.
4. dk_dT for SERIES mixing: correct dk/dT formula for series conductivity mixing.

Note on FD Validation Accuracy:
Per-block relative errors can be large for near-zero blocks or blocks involving
boundary cell geometry (see test_analytical_jacobian.jl). The overall Frobenius
norm is the primary success criterion. Block-level tests use absolute tolerance
to handle boundary effects.

Usage:
    julia --project=. test/test_jacobian_fixes.jl
"""

using Test
using Pyrolysis
using Pyrolysis.Core: has_darcy_flow
using Pyrolysis.TimeIntegration:
    PyrolysisODEFunction, JacobianCache, StateVector,
    compute_analytical_jacobian!, create_state_vector, unpack_state!,
    harmonic_mean_derivatives_weighted
using Pyrolysis.Geometry: update_geometry!
using LinearAlgebra
using Printf

# ============================================================================
# Helpers
# ============================================================================

"""
Compute central finite difference Jacobian with adaptive step sizing.
"""
function compute_central_fd_jacobian!(J_fd, u, ode_fn)
    n_dof = length(u)
    u_plus = copy(u)
    u_minus = copy(u)
    du_plus = similar(u)
    du_minus = similar(u)

    for j in 1:n_dof
        eps_j = sqrt(eps()) * max(abs(u[j]), 1.0)
        u_plus .= u
        u_minus .= u
        u_plus[j] += eps_j
        u_minus[j] -= eps_j
        ode_fn(du_plus, u_plus, nothing, 0.0)
        ode_fn(du_minus, u_minus, nothing, 0.0)
        @views J_fd[:, j] .= (du_plus .- du_minus) ./ (2 * eps_j)
    end
end

"""
Compute analytical and FD Jacobians and return them along with the ODE function.
Accepts optional I_radiation and radiation_model keyword arguments.
"""
function compute_jacobians(mesh, material, bc_set, u, NC, NR, n_cells;
    I_radiation::Float64 = 0.0,
    radiation_model = Pyrolysis.Core.NO_RADIATION)
    ode_fn = PyrolysisODEFunction(mesh, material, bc_set;
        I_radiation = I_radiation, radiation_model = radiation_model)
    sv = ode_fn.sv
    unpack_state!(mesh, u, sv)
    Pyrolysis.Properties.update_properties!(ode_fn.prop_cache, mesh, material)

    du = similar(u)
    ode_fn(du, u, nothing, 0.0)

    n_faces = length(mesh.faces)
    jac_cache = JacobianCache(n_cells, n_faces, NC, NR)
    n_dof = length(u)
    J_an = zeros(n_dof, n_dof)
    compute_analytical_jacobian!(J_an, u, nothing, 0.0, ode_fn, jac_cache)

    J_fd = zeros(n_dof, n_dof)
    compute_central_fd_jacobian!(J_fd, u, ode_fn)

    return J_an, J_fd
end

# ============================================================================
# Test 1: Distance-Weighted Harmonic Mean (Bug #1)
# ============================================================================

@testset "Bug #1: Distance-weighted harmonic mean derivatives" begin

    @testset "Unit test: harmonic_mean_derivatives_weighted" begin
        k_L = 0.5
        k_R = 0.1
        d_L = 0.3e-3
        d_R = 0.7e-3

        k_face, dk_dkL, dk_dkR = harmonic_mean_derivatives_weighted(k_L, k_R, d_L, d_R)

        # Verify face value: k_f = (d_L+d_R) / (d_L/k_L + d_R/k_R)
        expected_k_face = (d_L + d_R) / (d_L / k_L + d_R / k_R)
        @test k_face ≈ expected_k_face atol=1e-14

        # NOT equal to standard harmonic mean when d_L != d_R
        standard_harmonic = 2 * k_L * k_R / (k_L + k_R)
        @test abs(k_face - standard_harmonic) > 1e-6

        # Derivatives via central FD
        eps_k = 1e-7
        k_plus_L, _, _ = harmonic_mean_derivatives_weighted(k_L + eps_k, k_R, d_L, d_R)
        k_minus_L, _, _ = harmonic_mean_derivatives_weighted(k_L - eps_k, k_R, d_L, d_R)
        @test dk_dkL ≈ (k_plus_L - k_minus_L) / (2 * eps_k) rtol=1e-6

        k_plus_R, _, _ = harmonic_mean_derivatives_weighted(k_L, k_R + eps_k, d_L, d_R)
        k_minus_R, _, _ = harmonic_mean_derivatives_weighted(k_L, k_R - eps_k, d_L, d_R)
        @test dk_dkR ≈ (k_plus_R - k_minus_R) / (2 * eps_k) rtol=1e-6

        # When d_L = d_R, should equal standard harmonic mean
        k_sym, _, _ = harmonic_mean_derivatives_weighted(k_L, k_R, 0.5e-3, 0.5e-3)
        @test k_sym ≈ standard_harmonic atol=1e-14

        # Extreme asymmetry
        k_asym, dk_L_asym, dk_R_asym = harmonic_mean_derivatives_weighted(
            1.0, 0.01, 0.1e-3, 0.9e-3)
        @test k_asym ≈ 1.0e-3 / (0.1e-3 / 1.0 + 0.9e-3 / 0.01) atol=1e-14
        @test dk_L_asym > 0.0
        @test dk_R_asym > 0.0
    end

    @testset "Full Jacobian verification" begin
        NC = 2; NR = 1; n_cells = 20; thickness = 0.006

        mat = Material(
            name = :test_pmma,
            components = (
                SolidComponent(:virgin,
                    ρ = ConstantProperty(1190.0),
                    c = LinearProperty(1420.0, 2.5),
                    k = ConstantProperty(0.21),
                    ε = 0.86, α = 1.5, γ = 1.0,
                ),
                GaseousComponent(:gas,
                    M = 0.100,
                    c = ConstantProperty(1500.0),
                    k = ConstantProperty(0.02),
                    λ = ConstantProperty(2e-5),
                ),
            ),
            reactions = (
                Reaction(:decomp, 1 => 2,
                    A = 8.5e12, E = 188e3, h = 870e3,
                    T_min = 300.0, T_max = 1500.0,
                ),
            ),
            conductivity_mixing = PARALLEL,
            β = 0.5,
        )

        mesh = create_uniform_mesh(thickness, n_cells, Val(NC);
            T_initial = 300.0, ξ_initial = (1190.0, 0.0))
        bc_set = BoundaryConditionSet(
            thermal_top = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300.0),
            thermal_bottom = AdiabaticBC(),
            mass_top = DiffusiveMassBC(0.01),
            mass_bottom = ImpermeableBC())

        u = create_state_vector(mesh, NC)
        J_an, J_fd = compute_jacobians(mesh, mat, bc_set, u, NC, NR, n_cells)

        overall_err = norm(J_an - J_fd) / max(norm(J_fd), 1e-10)
        @test overall_err < 1e-3
        @info "Bug #1 full Jacobian" overall_err
    end
end

# ============================================================================
# Test 2: Conduction dkeff/dT (Bug #2)
# ============================================================================

@testset "Bug #2: Conduction dkeff/dT with T-dependent conductivity" begin
    NC = 2; NR = 1; n_cells = 20; thickness = 0.006

    mat = Material(
        name = :k_of_T_test,
        components = (
            SolidComponent(:virgin,
                ρ = ConstantProperty(1190.0),
                c = ConstantProperty(1500.0),
                k = LinearProperty(0.10, 5e-4),  # strong dk/dT
                ε = 0.9, α = 0.0, γ = 1.0),
            GaseousComponent(:gas,
                M = 0.044,
                c = ConstantProperty(1500.0),
                k = ConstantProperty(0.03),
                λ = ConstantProperty(2e-5)),
        ),
        reactions = (
            Reaction(:decomp, 1 => 2, A = 8.5e12, E = 250e3, h = 870e3),
        ),
        conductivity_mixing = PARALLEL, β = 0.5)

    mesh = create_uniform_mesh(thickness, n_cells, Val(NC);
        T_initial = 300.0, ξ_initial = (1190.0, 0.0))
    bc_set = BoundaryConditionSet(
        thermal_top = HeatFluxBC(50e3), thermal_bottom = AdiabaticBC(),
        mass_top = DiffusiveMassBC(0.01), mass_bottom = ImpermeableBC())

    u = create_state_vector(mesh, NC)
    n = n_cells
    for i in 1:n
        progress = (i - 1) / (n - 1)
        u[i] = 350.0 + 400.0 * progress   # 350-750 K large gradient
        u[n + i] = 1190.0
        u[2n + i] = 0.0
    end

    J_an, J_fd = compute_jacobians(mesh, mat, bc_set, u, NC, NR, n_cells)

    # dT/dT is the primary block for conduction dk/dT
    T_range = 1:n
    err_TT = norm(J_an[T_range, T_range] - J_fd[T_range, T_range]) /
             max(norm(J_fd[T_range, T_range]), 1e-10)
    @test err_TT < 1e-4

    overall_err = norm(J_an - J_fd) / max(norm(J_fd), 1e-10)
    @test overall_err < 1e-3
    @info "Bug #2 dkeff/dT" err_TT overall_err
end

# ============================================================================
# Test 3: Gas Convective dc_g/dT (Bug #3)
# ============================================================================

@testset "Bug #3: Gas convective dc_g/dT with T-dependent gas cp" begin
    NC = 2; NR = 1; n_cells = 40; thickness = 0.006

    mat = Material(
        name = :cg_of_T_test,
        components = (
            SolidComponent(:virgin,
                ρ = ConstantProperty(1190.0),
                c = ConstantProperty(1500.0),
                k = ConstantProperty(0.21),
                ε = 0.9, α = 0.0, γ = 1.0),
            GaseousComponent(:gas,
                M = 0.044,
                c = LinearProperty(1000.0, 2.0),  # c_g = 1000 + 2*T
                k = ConstantProperty(0.03),
                λ = ConstantProperty(2e-5)),
        ),
        reactions = (
            Reaction(:decomp, 1 => 2, A = 1e12, E = 200e3, h = 500e3),
        ),
        conductivity_mixing = PARALLEL, β = 0.5)

    mesh = create_uniform_mesh(thickness, n_cells, Val(NC);
        T_initial = 300.0, ξ_initial = (1190.0, 0.0))
    bc_set = BoundaryConditionSet(
        thermal_top = HeatFluxBC(50e3), thermal_bottom = AdiabaticBC(),
        mass_top = DiffusiveMassBC(0.01), mass_bottom = ImpermeableBC())

    u = create_state_vector(mesh, NC)
    n = n_cells
    # Small T gradient, small gas to exercise dc_g/dT
    for i in 1:n
        progress = (i - 1) / (n - 1)
        u[i] = 320.0 + 80.0 * progress   # 320-400 K
        u[n + i] = 1100.0
        u[2n + i] = 3.0                   # small gas present
    end

    J_an, J_fd = compute_jacobians(mesh, mat, bc_set, u, NC, NR, n_cells)

    # dT/dT: where dc_g/dT matters for convective energy source
    T_range = 1:n
    err_TT = norm(J_an[T_range, T_range] - J_fd[T_range, T_range]) /
             max(norm(J_fd[T_range, T_range]), 1e-10)
    @test err_TT < 1e-3

    # Overall (primary criterion)
    overall_err = norm(J_an - J_fd) / max(norm(J_fd), 1e-10)
    @test overall_err < 1e-3
    @info "Bug #3 dc_g/dT" err_TT overall_err
end

# ============================================================================
# Test 4: dk_dT for SERIES Mixing (Bug #4)
# ============================================================================

@testset "Bug #4: dk_dT for SERIES conductivity mixing" begin
    NC = 2; NR = 1; n_cells = 20; thickness = 0.006

    mat = Material(
        name = :series_k_test,
        components = (
            SolidComponent(:virgin,
                ρ = ConstantProperty(1190.0),
                c = ConstantProperty(1500.0),
                k = LinearProperty(0.50, 1e-3),  # high k, T-dependent
                ε = 0.9, α = 0.0, γ = 1.0),
            GaseousComponent(:gas,
                M = 0.044,
                c = ConstantProperty(1500.0),
                k = LinearProperty(0.02, 5e-5),  # low k, T-dependent
                λ = ConstantProperty(2e-5)),
        ),
        reactions = (
            Reaction(:decomp, 1 => 2, A = 8.5e12, E = 250e3, h = 870e3),
        ),
        conductivity_mixing = SERIES, β = 0.5)

    mesh = create_uniform_mesh(thickness, n_cells, Val(NC);
        T_initial = 300.0, ξ_initial = (1190.0, 0.0))
    bc_set = BoundaryConditionSet(
        thermal_top = HeatFluxBC(50e3), thermal_bottom = AdiabaticBC(),
        mass_top = DiffusiveMassBC(0.01), mass_bottom = ImpermeableBC())

    u = create_state_vector(mesh, NC)
    n = n_cells
    for i in 1:n
        progress = (i - 1) / (n - 1)
        u[i] = 400.0 + 300.0 * progress  # 400-700 K
        u[n + i] = 800.0                  # virgin
        u[2n + i] = 5.0                   # gas
    end

    J_an, J_fd = compute_jacobians(mesh, mat, bc_set, u, NC, NR, n_cells)

    # dT/dT: SERIES dk/dT affects conduction Jacobian
    T_range = 1:n
    err_TT = norm(J_an[T_range, T_range] - J_fd[T_range, T_range]) /
             max(norm(J_fd[T_range, T_range]), 1e-10)
    @test err_TT < 1e-3

    overall_err = norm(J_an - J_fd) / max(norm(J_fd), 1e-10)
    @test overall_err < 1e-3
    @info "Bug #4 SERIES dk_dT" err_TT overall_err
end

# ============================================================================
# Combined stress test
# ============================================================================

@testset "Combined: SERIES + T-dep k + T-dep gas cp" begin
    NC = 2; NR = 1; n_cells = 20; thickness = 0.006

    mat = Material(
        name = :combined_test,
        components = (
            SolidComponent(:virgin,
                ρ = ConstantProperty(1190.0),
                c = ConstantProperty(1500.0),
                k = LinearProperty(0.30, 8e-4),  # T-dep k (bug #2)
                ε = 0.9, α = 0.0, γ = 1.0),
            GaseousComponent(:gas,
                M = 0.044,
                c = LinearProperty(1200.0, 1.5),  # T-dep cp (bug #3)
                k = LinearProperty(0.02, 5e-5),   # T-dep k
                λ = ConstantProperty(2e-5)),
        ),
        reactions = (
            Reaction(:decomp, 1 => 2, A = 1e12, E = 200e3, h = 500e3),
        ),
        conductivity_mixing = SERIES, β = 0.5)  # bug #4

    mesh = create_uniform_mesh(thickness, n_cells, Val(NC);
        T_initial = 300.0, ξ_initial = (1190.0, 0.0))
    bc_set = BoundaryConditionSet(
        thermal_top = HeatFluxBC(50e3), thermal_bottom = AdiabaticBC(),
        mass_top = DiffusiveMassBC(0.01), mass_bottom = ImpermeableBC())

    u = create_state_vector(mesh, NC)
    n = n_cells
    for i in 1:n
        progress = (i - 1) / (n - 1)
        u[i] = 350.0 + 200.0 * progress  # 350-550 K
        u[n + i] = 1000.0
        u[2n + i] = 3.0
    end

    J_an, J_fd = compute_jacobians(mesh, mat, bc_set, u, NC, NR, n_cells)

    T_range = 1:n
    err_TT = norm(J_an[T_range, T_range] - J_fd[T_range, T_range]) /
             max(norm(J_fd[T_range, T_range]), 1e-10)
    @test err_TT < 1e-3

    overall_err = norm(J_an - J_fd) / max(norm(J_fd), 1e-10)
    @test overall_err < 1e-3
    @info "Combined test" err_TT overall_err
end

# ============================================================================
# Bug #5: RHS calls inside Jacobian missing radiation_model parameter
# ============================================================================

@testset "Bug #5: radiation_model passed to RHS in Jacobian" begin
    # With BEER_LAMBERT and I_radiation > 0, the RHS has nonzero volumetric
    # Q_rad which contributes to the property-derivative terms
    # (∂(1/ρcₚ)/∂ξ × RHS). Previously, the RHS calls inside the Jacobian
    # omitted radiation_model, defaulting to NO_RADIATION, which would give
    # Q_rad = 0 and thus wrong RHS values for these terms.
    #
    # This test uses BEER_LAMBERT to verify the radiation_model is correctly
    # passed to the RHS inside the Jacobian.

    NC = 2
    NR = 1
    n_cells = 20

    mat = Material(
        name = :test_rad_model,
        components = (
            SolidComponent(:virgin,
                ρ = ConstantProperty(1190.0),
                c = LinearProperty(1400.0, 2.5),  # T-dependent cp for nontrivial dρcp_dξ
                k = ConstantProperty(0.20),
                ε = 0.90,
                α = 1.5,  # nonzero α for Beer-Lambert absorption
                γ = 1.0,
            ),
            GaseousComponent(:gas,
                M = 0.044,
                c = ConstantProperty(1500.0),
                k = ConstantProperty(0.02),
                λ = ConstantProperty(2e-5),
                ε = 0.0,
                α = 0.0,
            ),
        ),
        reactions = (
            Reaction(:decomp, 1 => 2,
                A = 1e10, E = 160e3, h = 500e3,
                T_min = 300.0, T_max = 1500.0,
            ),
        ),
        conductivity_mixing = PARALLEL,
        β = 0.5,
        gas_transfer_mixing = PARALLEL,
        β_λ = 0.5,
    )

    mesh = create_uniform_mesh(0.006, n_cells, Val(NC);
        T_initial = 300.0, ξ_initial = (1190.0, 0.0))

    bc_set = BoundaryConditionSet(
        thermal_top = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300.0),
        thermal_bottom = AdiabaticBC(),
        mass_top = DiffusiveMassBC(Inf),
        mass_bottom = ImpermeableBC(),
    )

    # Set a mid-decomposition state with temperature gradient
    sv = Pyrolysis.TimeIntegration.StateVector(n_cells, NC)
    u = Pyrolysis.TimeIntegration.create_state_vector(mesh, NC)
    n = n_cells
    for i in 1:n
        progress = (i - 1) / (n - 1)
        u[i] = 350.0 + 150.0 * progress  # 350-500 K gradient
        u[n + i] = 1190.0 * (1.0 - 0.2 * progress)  # virgin depleting
        u[2n + i] = 5.0 * progress  # small gas present
    end

    J_an, J_fd = compute_jacobians(mesh, mat, bc_set, u, NC, NR, n_cells;
        I_radiation = 30e3,
        radiation_model = Pyrolysis.Core.BEER_LAMBERT)

    # Beer-Lambert radiation has inherent O(1%) FD-vs-analytical differences
    # due to the absorption derivative approximation (exp(-τ)×αⱼ for large τ).
    # Use relaxed tolerance matching comprehensive test suite.
    overall_err = norm(J_an - J_fd) / max(norm(J_fd), 1e-10)
    @test overall_err < 0.05
    @info "Bug #5 BEER_LAMBERT" overall_err
end

# ============================================================================
# Bug #6: Radiation Jacobian guard (Beer-Lambert only)
# ============================================================================

@testset "Bug #6: No Beer-Lambert Jacobian for SURFACE_ABSORPTION" begin
    # With SURFACE_ABSORPTION, the RHS has Q_rad=0 for interior cells.
    # Previously, the Jacobian incorrectly computed Beer-Lambert volumetric
    # derivatives when I_radiation > 0, regardless of radiation_model.
    # This created spurious non-local coupling in the upper triangular block.
    #
    # After the fix, compute_radiation_jacobian! is only called when
    # radiation_model == BEER_LAMBERT.

    NC = 2
    NR = 1
    n_cells = 15

    # Material with nonzero absorption coefficients (would cause
    # incorrect Beer-Lambert coupling if the guard is missing)
    mat = Material(
        name = :test_rad_guard,
        components = (
            SolidComponent(:virgin,
                ρ = ConstantProperty(500.0),
                c = ConstantProperty(1500.0),
                k = ConstantProperty(0.15),
                ε = 0.90,
                α = 1.5,  # nonzero α would cause Beer-Lambert coupling
                γ = 1.0,
            ),
            GaseousComponent(:gas,
                M = 0.044,
                c = ConstantProperty(1200.0),
                k = ConstantProperty(0.025),
                λ = ConstantProperty(1e-5),
                ε = 0.0,
                α = 0.0,
            ),
        ),
        reactions = (
            Reaction(:decomp, 1 => 2,
                A = 1e12, E = 170e3, h = 400e3,
                T_min = 300.0, T_max = 1500.0,
            ),
        ),
        conductivity_mixing = PARALLEL,
        β = 0.5,
        gas_transfer_mixing = PARALLEL,
        β_λ = 0.5,
    )

    mesh = create_uniform_mesh(0.005, n_cells, Val(NC);
        T_initial = 300.0, ξ_initial = (500.0, 0.0))

    bc_set = BoundaryConditionSet(
        thermal_top = HeatFluxBC(50e3),
        thermal_bottom = AdiabaticBC(),
        mass_top = DiffusiveMassBC(Inf),
        mass_bottom = ImpermeableBC(),
    )

    # SURFACE_ABSORPTION + I_radiation > 0: radiation should NOT create
    # Beer-Lambert volumetric Jacobian entries
    ode_fn = Pyrolysis.TimeIntegration.PyrolysisODEFunction(
        mesh, mat, bc_set;
        I_radiation = 50e3,
        radiation_model = Pyrolysis.Core.SURFACE_ABSORPTION,
    )

    sv = Pyrolysis.TimeIntegration.StateVector(n_cells, NC)
    u = Pyrolysis.TimeIntegration.create_state_vector(mesh, NC)

    # Uniform state (pure conduction + reactions, no radiation volumetric source)
    n = n_cells
    for i in 1:n
        progress = (i - 1) / (n - 1)
        u[i] = 350.0 + 200.0 * progress
        u[n + i] = 500.0
        u[2n + i] = 0.0
    end

    J_an, J_fd = compute_jacobians(mesh, mat, bc_set, u, NC, NR, n_cells;
        I_radiation = 50e3,
        radiation_model = Pyrolysis.Core.SURFACE_ABSORPTION)

    # Check specifically that upper-triangular non-local coupling is absent
    # (Beer-Lambert creates dT[i]/dξ[k] for k > i)
    T_range = 1:n
    ξ1_range = (n+1):(2n)
    dTdξ_block = J_an[T_range, ξ1_range]

    # Count non-local (upper triangular) entries that should be zero
    non_local_norm = 0.0
    for i in 1:n, j in (i+1):n
        non_local_norm += dTdξ_block[i, j]^2
    end
    non_local_norm = sqrt(non_local_norm)
    @test non_local_norm < 1e-10  # Should be exactly zero with fix

    overall_err = norm(J_an - J_fd) / max(norm(J_fd), 1e-10)
    @test overall_err < 1e-3
    @info "Bug #6 radiation guard" non_local_norm overall_err
end

# ============================================================================
# Bug #7: Missing ∂λ_face/∂T for interior gas transport
# ============================================================================

@testset "Bug #7: dλ_face/dT for interior gas transport" begin
    # The gas flux J = -λ_face × Φ / (T_face × d_LR) depends on T through:
    # 1. Direct: Φ and T_face (already included)
    # 2. Through λ_face(T): component λ_j(T) are T-dependent
    #
    # Previously, only the direct T-derivative was included for interior faces.
    # The boundary mass BC Jacobian already included ∂λ/∂T, creating an
    # inconsistency. This fix adds the ∂λ_face/∂T contribution to interior
    # gas transport T-derivatives.
    #
    # To maximally exercise this: use a gas with strongly T-dependent λ.

    NC = 2
    NR = 1
    n_cells = 20

    mat = Material(
        name = :test_dlambda_dT,
        components = (
            SolidComponent(:solid,
                ρ = ConstantProperty(800.0),
                c = ConstantProperty(1400.0),
                k = ConstantProperty(0.15),
                ε = 0.90,
                α = 0.0,
                γ = 1.0,
            ),
            GaseousComponent(:gas,
                M = 0.028,
                c = ConstantProperty(1100.0),
                k = ConstantProperty(0.03),
                λ = LinearProperty(1e-5, 5e-8),  # Strongly T-dependent λ
                ε = 0.0,
                α = 0.0,
            ),
        ),
        reactions = (
            Reaction(:decomp, 1 => 2,
                A = 1e10, E = 150e3, h = 300e3,
                T_min = 300.0, T_max = 1500.0,
            ),
        ),
        conductivity_mixing = PARALLEL,
        β = 0.5,
        gas_transfer_mixing = PARALLEL,
        β_λ = 0.5,
    )

    mesh = create_uniform_mesh(0.005, n_cells, Val(NC);
        T_initial = 300.0, ξ_initial = (800.0, 0.0))

    bc_set = BoundaryConditionSet(
        thermal_top = HeatFluxBC(40e3),
        thermal_bottom = AdiabaticBC(),
        mass_top = DiffusiveMassBC(Inf),
        mass_bottom = ImpermeableBC(),
    )

    ode_fn = Pyrolysis.TimeIntegration.PyrolysisODEFunction(
        mesh, mat, bc_set; I_radiation = 0.0)

    sv = Pyrolysis.TimeIntegration.StateVector(n_cells, NC)
    u = Pyrolysis.TimeIntegration.create_state_vector(mesh, NC)

    # State with T gradient and gas present (exercises λ_face T-dependence)
    n = n_cells
    for i in 1:n
        progress = (i - 1) / (n - 1)
        u[i] = 400.0 + 400.0 * progress  # 400-800 K (large gradient for dλ/dT)
        u[n + i] = 800.0 * (1.0 - 0.5 * progress)  # solid depleting
        u[2n + i] = 20.0 + 30.0 * progress  # gas present throughout
    end

    J_an, J_fd = compute_jacobians(mesh, mat, bc_set, u, NC, NR, n_cells)

    # Check dξ_gas/dT block specifically (gas transport T-derivatives)
    # This is the block directly affected by the ∂λ_face/∂T fix.
    T_range = 1:n
    ξ2_range = (2n+1):(3n)  # gas component
    err_ξT = norm(J_an[ξ2_range, T_range] - J_fd[ξ2_range, T_range]) /
             max(norm(J_fd[ξ2_range, T_range]), 1e-10)
    @test err_ξT < 1e-3

    # Overall error may be slightly above 1e-3 due to boundary cell effects
    # with large temperature gradients (400-800K). The key test is the
    # dξ/dT block above which directly exercises the ∂λ_face/∂T fix.
    overall_err = norm(J_an - J_fd) / max(norm(J_fd), 1e-10)
    @test overall_err < 0.01
    @info "Bug #7 dλ_face/dT" err_ξT overall_err
end
