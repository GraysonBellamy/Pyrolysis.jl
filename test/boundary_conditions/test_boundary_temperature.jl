"""
Boundary Temperature Solver Tests

Tests for the boundary temperature solve functionality that ensures energy balance
at boundaries in cell-centered finite volume discretization.

Reference: docs/boundary_temperature_enhancement.md

Tests cover:
1. Unit tests for solve_boundary_temperature() with each BC type
2. Helper function tests (get_surface_temperature, etc.)
3. Energy balance verification
4. Integration tests with solver
"""

using Test
using Pyrolysis
using Pyrolysis.Materials
using Pyrolysis.BoundaryConditions
using Pyrolysis.Geometry
using Pyrolysis.Discretization

# Import specific functions for testing
using Pyrolysis.Discretization:
    solve_boundary_temperature,
    BoundaryState,
    get_surface_temperature,
    get_back_face_temperature,
    check_boundary_energy_balance,
    update_boundary_state!,
    evaluate_bc_flux_and_derivative
using Pyrolysis.Materials: T_REF, GaseousComponent, SolidComponent, Material, Reaction
using Pyrolysis.BoundaryConditions: ambient_temperature

# Test helper: constructs a tiny inert single-solid material so that the
# transpiration loop is well-defined but contributes zero (no gas component).
# Tests below pass this material with `zero_J = ()` to recover the
# no-blowing limit — equivalent to the no-blowing boundary solve and
# convenient for closed-form analytic verification.
const _BC_TEST_SOLID = SolidComponent(:test_solid; ρ = 1000.0, c = 1000.0, k = 0.2)
const _BC_TEST_MAT = Material(name = :bc_test_mat, components = (_BC_TEST_SOLID,))
const _BC_ZERO_J = ntuple(_ -> 0.0, Val(1))

# Wrapper: no-blowing call, T_ref resolved from BC ambient (or T_REF fallback).
_solve_bt(T_cell, k, Δz, bc, t, ε; kwargs...) = solve_boundary_temperature(
    T_cell,
    k,
    Δz,
    bc,
    _BC_TEST_MAT,
    _BC_ZERO_J,
    ambient_temperature(bc, t, Float64(T_REF)),
    t,
    ε;
    kwargs...,
)

@testset "Boundary Temperature Solver" begin

    # ========================================================================
    # Unit Tests for solve_boundary_temperature()
    # ========================================================================
    @testset "solve_boundary_temperature" begin

        @testset "Pure Heat Flux BC" begin
            # With prescribed flux, T_s should satisfy: q_ext = k(T_s - T_cell)/Δz
            # Rearranging: T_s = T_cell + q_ext * Δz_half / k
            T_cell = 500.0
            k_eff = 0.2
            Δz_half = 50e-6
            q_ext = 10000.0  # W/m² into domain
            ε_eff = 0.9
            bc = HeatFluxBC(_ -> q_ext)

            state = _solve_bt(T_cell, k_eff, Δz_half, bc, 0.0, ε_eff)

            @test state.converged
            @test state.iterations <= 5  # Should converge quickly (linear problem)
            @test state.q_conducted ≈ q_ext rtol=1e-8

            # Expected: T_s = T_cell + q_ext * Δz_half / k_eff
            # (positive flux means heat into domain, so T_s > T_cell)
            T_s_expected = T_cell + q_ext * Δz_half / k_eff
            @test state.T_surface ≈ T_s_expected rtol=1e-8

            # Verify stored values
            @test state.T_cell ≈ T_cell
            @test state.k_eff ≈ k_eff
            @test state.Δz_half ≈ Δz_half
        end

        @testset "Constant Heat Flux BC" begin
            # Test with constant (non-function) heat flux
            T_cell = 400.0
            k_eff = 0.25
            Δz_half = 100e-6
            q_ext = 25000.0
            ε_eff = 0.9
            bc = HeatFluxBC(q_ext)

            state = _solve_bt(T_cell, k_eff, Δz_half, bc, 0.0, ε_eff)

            @test state.converged
            @test state.q_applied ≈ q_ext rtol=1e-8
            @test state.q_conducted ≈ q_ext rtol=1e-8

            # T_s = T_cell + q * Δz / k
            T_s_expected = T_cell + q_ext * Δz_half / k_eff
            @test state.T_surface ≈ T_s_expected rtol=1e-8
        end

        @testset "Convective BC" begin
            # Convective BC: q = h(T_inf - T_s)
            # Energy balance: h(T_inf - T_s) = k(T_s - T_cell)/Δz
            T_cell = 600.0
            k_eff = 0.2
            Δz_half = 50e-6
            h = 10.0
            T_inf = 300.0
            ε_eff = 0.9
            bc = ConvectiveBC(h = h, T_inf = T_inf)

            state = _solve_bt(T_cell, k_eff, Δz_half, bc, 0.0, ε_eff)

            @test state.converged
            @test state.iterations <= 5  # Linear problem
            @test state.q_applied ≈ state.q_conducted rtol=1e-8

            # Analytical solution: h(T_inf - T_s) = k(T_s - T_cell)/Δz
            # h*T_inf - h*T_s = k*T_s/Δz - k*T_cell/Δz
            # h*T_inf + k*T_cell/Δz = T_s*(h + k/Δz)
            # T_s = (h*T_inf + k*T_cell/Δz) / (h + k/Δz)
            T_s_expected = (h * T_inf + k_eff * T_cell / Δz_half) / (h + k_eff / Δz_half)
            @test state.T_surface ≈ T_s_expected rtol=1e-8

            # Surface should be between T_inf and T_cell (cooling case)
            @test T_inf < state.T_surface < T_cell
        end

        @testset "Radiative BC" begin
            # Radiative BC: q = εσ(T_inf⁴ - T_s⁴)
            T_cell = 800.0
            k_eff = 0.2
            Δz_half = 50e-6
            ε = 0.9
            T_inf = 300.0
            bc = RadiativeBC(ε = ε, T_inf = T_inf)

            state = _solve_bt(T_cell, k_eff, Δz_half, bc, 0.0, ε)

            @test state.converged
            @test state.q_applied ≈ state.q_conducted rtol=1e-6

            # Verify flux formula: q = εσ(T_inf⁴ - T_s⁴)
            σ = STEFAN_BOLTZMANN
            q_expected = ε * σ * (T_inf^4 - state.T_surface^4)
            @test state.q_applied ≈ q_expected rtol=1e-8

            # Surface should be between T_inf and T_cell (cooling case)
            @test T_inf < state.T_surface < T_cell
        end

        @testset "Combined BC (flux + radiative + convective)" begin
            # Typical cone calorimeter setup
            T_cell = 850.0
            k_eff = 0.18
            Δz_half = 23e-6

            q_ext = 50000.0   # 50 kW/m² external heat flux
            ε = 0.9
            h = 7.2           # Convective coefficient
            T_amb = 330.0     # Ambient temperature

            bc =
                HeatFluxBC(_ -> q_ext) +
                RadiativeBC(ε = ε, T_inf = T_amb) +
                ConvectiveBC(h = h, T_inf = T_amb)

            state = _solve_bt(T_cell, k_eff, Δz_half, bc, 0.0, ε)

            @test state.converged
            @test state.iterations <= 20  # Nonlinear but well-behaved
            @test abs(state.q_applied - state.q_conducted) < 1e-6

            # With net heating, T_surface > T_cell
            @test state.T_surface > T_cell
        end

        @testset "Adiabatic BC" begin
            # Adiabatic: q = 0, so T_s = T_cell (no gradient needed)
            T_cell = 500.0
            k_eff = 0.2
            Δz_half = 50e-6
            ε_eff = 0.9
            bc = AdiabaticBC()

            state = _solve_bt(T_cell, k_eff, Δz_half, bc, 0.0, ε_eff)

            @test state.converged
            @test state.iterations == 1  # Should converge immediately
            @test state.q_applied ≈ 0.0 atol=1e-12
            @test state.q_conducted ≈ 0.0 atol=1e-12
            @test state.T_surface ≈ T_cell rtol=1e-10
        end

        @testset "Dirichlet BC" begin
            # Dirichlet: T_s is prescribed, solver should return it directly
            T_cell = 500.0
            k_eff = 0.2
            Δz_half = 50e-6
            T_prescribed = 800.0
            ε_eff = 0.9
            bc = DirichletThermalBC(T_prescribed)

            state = _solve_bt(T_cell, k_eff, Δz_half, bc, 0.0, ε_eff)

            @test state.converged
            @test state.T_surface ≈ T_prescribed rtol=1e-4  # Allow small numerical tolerance

            # Conducted flux should be consistent with temperature difference
            q_expected = k_eff * (T_prescribed - T_cell) / Δz_half
            @test state.q_conducted ≈ q_expected rtol=1e-4  # Allow small numerical tolerance
        end

        @testset "Time-dependent BC" begin
            # Heat flux that varies with time
            T_cell = 400.0
            k_eff = 0.2
            Δz_half = 50e-6
            ε_eff = 0.9

            # Ramping flux: q(t) = 1000*t for t in [0, 50]
            bc = HeatFluxBC(t -> 1000.0 * t)

            # At t=0
            state0 = _solve_bt(T_cell, k_eff, Δz_half, bc, 0.0, ε_eff)
            @test state0.q_applied ≈ 0.0 atol=1e-10
            @test state0.T_surface ≈ T_cell rtol=1e-10

            # At t=10
            state10 = _solve_bt(T_cell, k_eff, Δz_half, bc, 10.0, ε_eff)
            @test state10.q_applied ≈ 10000.0 rtol=1e-10
            @test state10.T_surface > T_cell

            # At t=50
            state50 = _solve_bt(T_cell, k_eff, Δz_half, bc, 50.0, ε_eff)
            @test state50.q_applied ≈ 50000.0 rtol=1e-10
            @test state50.T_surface > state10.T_surface
        end

        @testset "Transpiration enthalpy sign" begin
            solid = SolidComponent(:test_solid; ρ = 1000.0, c = 1000.0, k = 0.2)
            gas = GaseousComponent(:test_gas; M = 0.044, c = 1500.0, k = 0.03, λ = 1e-5)
            material = Material(name = :transpiration_sign_test, components = (solid, gas))

            T_cell = 700.0
            T_ref = 300.0
            k_eff = 0.2
            Δz_half = 50e-6
            q_ext = 30000.0
            ε_eff = 0.9
            bc = HeatFluxBC(q_ext)

            no_blowing = solve_boundary_temperature(
                T_cell,
                k_eff,
                Δz_half,
                bc,
                material,
                (0.0, 0.0),
                T_ref,
                0.0,
                ε_eff,
            )
            outward_hot_gas = solve_boundary_temperature(
                T_cell,
                k_eff,
                Δz_half,
                bc,
                material,
                (0.0, 0.04),
                T_ref,
                0.0,
                ε_eff,
            )
            inward_hot_gas = solve_boundary_temperature(
                T_cell,
                k_eff,
                Δz_half,
                bc,
                material,
                (0.0, -0.04),
                T_ref,
                0.0,
                ε_eff,
            )

            @test no_blowing.converged
            @test outward_hot_gas.converged
            @test inward_hot_gas.converged

            # Positive outward gas carries sensible enthalpy away from the surface,
            # reducing the conducted heat needed to balance the imposed flux.
            @test outward_hot_gas.T_surface < no_blowing.T_surface
            @test outward_hot_gas.q_conducted < no_blowing.q_conducted
            @test outward_hot_gas.q_applied - outward_hot_gas.q_conducted > 0.0

            # Negative "outward" flux is gas entering from the ambient side. The
            # transpiration term is UPWINDED to the donor temperature: entering
            # gas crosses the face at T_∞ = T_ref carrying Δh_g(T_ref,T_ref) = 0,
            # so it makes no contribution to the surface balance (it is heated to
            # T_s volumetrically inside). Inflow therefore leaves the surface
            # temperature and conducted flux unchanged from the no-blowing case.
            @test inward_hot_gas.T_surface ≈ no_blowing.T_surface
            @test inward_hot_gas.q_conducted ≈ no_blowing.q_conducted
            @test inward_hot_gas.q_applied - inward_hot_gas.q_conducted ≈ 0.0 atol = 1e-9
            @test inward_hot_gas.q_transpired ≈ 0.0 atol = 1e-9
            # Outward gas, by contrast, records a positive transpiration sink.
            @test outward_hot_gas.q_transpired > 0.0
        end

    end

    # ========================================================================
    # Convergence and Robustness Tests
    # ========================================================================
    @testset "Convergence Robustness" begin

        @testset "High temperature radiative BC" begin
            # Extreme case: very hot surface with strong radiation
            T_cell = 2000.0
            k_eff = 0.5
            Δz_half = 100e-6
            ε = 0.95
            T_inf = 300.0
            bc = RadiativeBC(ε = ε, T_inf = T_inf)

            state = _solve_bt(T_cell, k_eff, Δz_half, bc, 0.0, ε)

            @test state.converged
            @test 300.0 < state.T_surface < 2000.0
        end

        @testset "Very thin boundary layer" begin
            # Small Δz_half means high conductance
            T_cell = 600.0
            k_eff = 0.2
            Δz_half = 1e-6  # 1 micron
            q_ext = 50000.0
            ε_eff = 0.9
            bc = HeatFluxBC(q_ext)

            state = _solve_bt(T_cell, k_eff, Δz_half, bc, 0.0, ε_eff)

            @test state.converged
            # T_s should be very close to T_cell (high conductance)
            @test abs(state.T_surface - T_cell) < 1.0
        end

        @testset "Very thick boundary layer" begin
            # Large Δz_half means low conductance
            T_cell = 600.0
            k_eff = 0.2
            Δz_half = 1e-3  # 1 mm
            q_ext = 50000.0
            ε_eff = 0.9
            bc = HeatFluxBC(q_ext)

            state = _solve_bt(T_cell, k_eff, Δz_half, bc, 0.0, ε_eff)

            @test state.converged
            # Large temperature jump expected
            @test state.T_surface > T_cell + 100.0
        end

    end

    # ========================================================================
    # RadiativeFluxBC with Radiation Models
    # ========================================================================
    @testset "RadiativeFluxBC with Radiation Models" begin
        # Test that RadiativeFluxBC correctly returns zero flux for volumetric
        # radiation models (SURFACE_ABSORPTION, BEER_LAMBERT, P1_QUASI_STEADY)
        # because radiation is handled via compute_radiation_source!, not as
        # a boundary flux.
        #
        # Reference: docs/boundary_temperature_enhancement.md

        T_s = 500.0
        t = 0.0
        ε_eff = 0.9
        I_incident = 50e3  # 50 kW/m² incident radiation

        bc = RadiativeFluxBC(I_incident)

        @testset "NO_RADIATION: flux applied at boundary" begin
            # With NO_RADIATION, RadiativeFluxBC flux is treated as boundary heat flux
            # The absorbed flux is I_incident * ε_eff (Kirchhoff's law: α = ε)
            q, dq = evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff, NO_RADIATION)

            @test q ≈ I_incident * ε_eff  # Absorbed flux = incident × absorptivity
            @test dq ≈ 0.0                # Flux doesn't depend on T_s
        end

        @testset "SURFACE_ABSORPTION: flux applied at boundary" begin
            # With SURFACE_ABSORPTION, radiation is absorbed at the boundary surface.
            # This ensures T_surface is physically correct (hotter than T_cell when net heating).
            # Same behavior as NO_RADIATION for boundary temperature purposes.
            q, dq = evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff, SURFACE_ABSORPTION)

            @test q ≈ I_incident * ε_eff  # Absorbed flux = incident × absorptivity
            @test dq ≈ 0.0                # Flux doesn't depend on T_s
        end

        @testset "BEER_LAMBERT: flux handled volumetrically" begin
            # With BEER_LAMBERT, radiation penetrates in-depth and is absorbed
            # according to Beer-Lambert law via compute_radiation_source!
            q, dq = evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff, BEER_LAMBERT)

            @test q ≈ 0.0   # No boundary flux contribution
            @test dq ≈ 0.0  # No temperature dependence
        end

        @testset "P1_QUASI_STEADY: flux handled by P1 solver" begin
            # With P1_QUASI_STEADY, radiation is solved via P1 approximation
            # with Marshak boundary conditions - not via boundary flux
            q, dq = evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff, P1_QUASI_STEADY)

            @test q ≈ 0.0   # No boundary flux contribution
            @test dq ≈ 0.0  # No temperature dependence
        end

        @testset "RadiativeFluxBC with explicit absorptivity" begin
            # When absorptivity is explicitly specified in RadiativeFluxBC,
            # it overrides the effective emissivity
            α_explicit = 0.7
            bc_with_α = RadiativeFluxBC(flux = I_incident, absorptivity = α_explicit)

            # NO_RADIATION uses explicit absorptivity
            q, dq = evaluate_bc_flux_and_derivative(bc_with_α, T_s, t, ε_eff, NO_RADIATION)
            @test q ≈ I_incident * α_explicit  # Uses explicit α, not ε_eff
            @test dq ≈ 0.0

            # SURFACE_ABSORPTION also uses explicit absorptivity (same as NO_RADIATION)
            q_surf, _ = evaluate_bc_flux_and_derivative(
                bc_with_α,
                T_s,
                t,
                ε_eff,
                SURFACE_ABSORPTION,
            )
            @test q_surf ≈ I_incident * α_explicit

            # Volumetric models (BEER_LAMBERT, P1) return zero
            q_bl, _ =
                evaluate_bc_flux_and_derivative(bc_with_α, T_s, t, ε_eff, BEER_LAMBERT)
            q_p1, _ =
                evaluate_bc_flux_and_derivative(bc_with_α, T_s, t, ε_eff, P1_QUASI_STEADY)

            @test q_bl ≈ 0.0
            @test q_p1 ≈ 0.0
        end

        @testset "Combined BC with RadiativeFluxBC" begin
            # When RadiativeFluxBC is combined with other BCs
            h_conv = 10.0
            T_inf = 300.0

            bc_combined =
                RadiativeFluxBC(I_incident) + ConvectiveBC(h = h_conv, T_inf = T_inf)

            # For NO_RADIATION: both components contribute
            q_no_rad, dq_no_rad =
                evaluate_bc_flux_and_derivative(bc_combined, T_s, t, ε_eff, NO_RADIATION)
            q_rad_expected = I_incident * ε_eff
            q_conv_expected = h_conv * (T_inf - T_s)
            @test q_no_rad ≈ q_rad_expected + q_conv_expected
            @test dq_no_rad ≈ -h_conv  # Only convective contributes to derivative

            # For SURFACE_ABSORPTION: both components also contribute (same as NO_RADIATION)
            q_surf, dq_surf = evaluate_bc_flux_and_derivative(
                bc_combined,
                T_s,
                t,
                ε_eff,
                SURFACE_ABSORPTION,
            )
            @test q_surf ≈ q_rad_expected + q_conv_expected
            @test dq_surf ≈ -h_conv

            # For BEER_LAMBERT: only convective contributes (radiation handled volumetrically)
            q_bl, dq_bl =
                evaluate_bc_flux_and_derivative(bc_combined, T_s, t, ε_eff, BEER_LAMBERT)
            @test q_bl ≈ q_conv_expected
            @test dq_bl ≈ -h_conv
        end

    end

    # ========================================================================
    # Helper Function Tests
    # ========================================================================
    @testset "Helper Functions" begin

        @testset "get_surface_temperature and get_back_face_temperature" begin
            # Create a mesh with boundary states initialized
            L = 0.01
            n_cells = 10
            mesh = create_uniform_mesh(L, n_cells, Val(2))

            # Manually set boundary states for testing
            T_top = 800.0
            T_bottom = 350.0
            mesh.boundary_states[1] = Pyrolysis.Geometry.BoundaryState(
                T_top,
                50000.0,
                50000.0,
                0.0,
                750.0,
                0.2,
                50e-6,
                true,
                3,
            )
            mesh.boundary_states[2] = Pyrolysis.Geometry.BoundaryState(
                T_bottom,
                0.0,
                0.0,
                0.0,
                350.0,
                0.2,
                50e-6,
                true,
                1,
            )

            @test get_surface_temperature(mesh) ≈ T_top
            @test get_back_face_temperature(mesh) ≈ T_bottom
        end

        @testset "check_boundary_energy_balance" begin
            # Create a mesh with boundary states for testing
            L = 0.01
            n_cells = 10
            mesh_test = create_uniform_mesh(L, n_cells, Val(2))

            # Set up a well-balanced boundary state
            mesh_test.boundary_states[1] = Pyrolysis.Geometry.BoundaryState(
                600.0,    # T_surface
                10000.0,  # q_applied
                10000.0,  # q_conducted (matching!)
                0.0,      # q_transpired
                550.0,    # T_cell
                0.2,      # k_eff
                50e-6,    # Δz_half
                true,     # converged
                3,         # iterations
            )
            mesh_test.boundary_states[2] = Pyrolysis.Geometry.BoundaryState(
                350.0,
                0.0,
                0.0,
                0.0,
                350.0,
                0.2,
                50e-6,
                true,
                1,
            )

            # check_boundary_energy_balance returns Dict{Int, Float64}
            balance_results = check_boundary_energy_balance(mesh_test)
            @test balance_results[1] < 1e-10
            @test balance_results[2] < 1e-10

            # Now test with imbalance
            mesh_test.boundary_states[1] = Pyrolysis.Geometry.BoundaryState(
                600.0,    # T_surface
                10000.0,  # q_applied
                8000.0,   # q_conducted (20% error!)
                0.0,      # q_transpired
                550.0,    # T_cell
                0.2,      # k_eff
                50e-6,    # Δz_half
                true,     # converged
                3,         # iterations
            )

            balance_results_bad = check_boundary_energy_balance(mesh_test)
            # 20% imbalance
            @test balance_results_bad[1] ≈ 0.2 atol=0.01
        end

    end

    # ========================================================================
    # evaluate_bc_flux_and_derivative Tests
    # ========================================================================
    @testset "BC Flux and Derivative Evaluation" begin

        @testset "HeatFluxBC derivatives" begin
            bc = HeatFluxBC(50000.0)
            q, dq = evaluate_bc_flux_and_derivative(bc, 500.0, 0.0, 0.9)
            @test q ≈ 50000.0
            @test dq ≈ 0.0  # No T_s dependence
        end

        @testset "ConvectiveBC derivatives" begin
            h = 10.0
            T_inf = 300.0
            bc = ConvectiveBC(h = h, T_inf = T_inf)
            T_s = 500.0
            q, dq = evaluate_bc_flux_and_derivative(bc, T_s, 0.0, 0.9)

            @test q ≈ h * (T_inf - T_s)
            @test dq ≈ -h  # dq/dT_s = -h
        end

        @testset "RadiativeBC derivatives" begin
            ε = 0.9
            T_inf = 300.0
            bc = RadiativeBC(ε = ε, T_inf = T_inf)
            T_s = 600.0
            q, dq = evaluate_bc_flux_and_derivative(bc, T_s, 0.0, ε)

            σ = STEFAN_BOLTZMANN
            @test q ≈ ε * σ * (T_inf^4 - T_s^4)
            @test dq ≈ -4 * ε * σ * T_s^3 rtol=1e-8
        end

        @testset "AdiabaticBC derivatives" begin
            bc = AdiabaticBC()
            q, dq = evaluate_bc_flux_and_derivative(bc, 500.0, 0.0, 0.9)
            @test q ≈ 0.0
            @test dq ≈ 0.0
        end

    end

    # ========================================================================
    # Integration Tests
    # ========================================================================
    @testset "Integration with Solver" begin

        @testset "Boundary states updated during solve" begin
            # Create a simple problem and verify boundary states are populated
            L = 0.005
            n_cells = 20

            solid = SolidComponent(:solid, ρ = 1000.0, c = 2000.0, k = 0.2, ε = 0.9)
            mat = Material(name = :TestSolid, components = (solid,), reactions = ())

            mesh = create_uniform_mesh(
                L,
                n_cells,
                Val(1);
                T_initial = 300.0,
                ξ_initial = (1000.0,),
            )

            bc_set = BoundaryConditionSet(
                thermal = (
                    top = HeatFluxBC(50000.0) + RadiativeBC(ε = 0.9, T_inf = 300.0),
                    bottom = AdiabaticBC(),
                ),
            )

            problem = PyrolysisProblem(
                mesh = mesh,
                material = mat,
                bc_set = bc_set,
                T_initial = z -> 300.0,
                ξ_initial = [z -> 1000.0],
                tspan = (0.0, 1.0),
            )

            solution = solve(problem; progress = false)

            @test solution.retcode == :Success

            # Verify boundary states were updated
            @test mesh.boundary_states[1].converged
            @test mesh.boundary_states[2].converged

            # Surface temperature should be different from cell temperature
            T_surface = mesh.boundary_states[1].T_surface
            T_cell_top = mesh.cell_states[n_cells].T
            @test T_surface != T_cell_top

            # With 50 kW/m² heating, surface should be hotter than cell center
            @test T_surface > T_cell_top

            # Energy balance should be satisfied - check_boundary_energy_balance returns Dict
            balance_results = check_boundary_energy_balance(mesh)
            @test balance_results[1] < 1e-6
            @test balance_results[2] < 1e-6
        end

        @testset "Solution output uses boundary temperature" begin
            # Verify that solution.surface_T uses T_surface, not T_cell
            L = 0.005
            n_cells = 20

            solid = SolidComponent(:solid, ρ = 1000.0, c = 2000.0, k = 0.2, ε = 0.9)
            mat = Material(name = :TestSolid, components = (solid,), reactions = ())

            mesh = create_uniform_mesh(
                L,
                n_cells,
                Val(1);
                T_initial = 300.0,
                ξ_initial = (1000.0,),
            )

            bc_set = BoundaryConditionSet(
                thermal = (top = HeatFluxBC(75000.0), bottom = AdiabaticBC()),
            )

            problem = PyrolysisProblem(
                mesh = mesh,
                material = mat,
                bc_set = bc_set,
                T_initial = z -> 300.0,
                ξ_initial = [z -> 1000.0],
                tspan = (0.0, 5.0),
            )

            solution = solve(problem; progress = false, saveat = [5.0])

            @test solution.retcode == :Success

            # Get final surface temperature from solution
            T_surf_output = solution.surface_T[end]

            # This should match the boundary state, not the top cell
            @test T_surf_output ≈ mesh.boundary_states[1].T_surface rtol=0.01

            # Surface should be noticeably hotter than top cell
            T_cell_top = mesh.cell_states[n_cells].T
            ΔT = T_surf_output - T_cell_top
            @test ΔT > 5.0  # Expect meaningful difference with 75 kW/m²
        end

    end

end  # main testset
