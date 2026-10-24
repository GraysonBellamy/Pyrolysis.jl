"""
Boundary Temperature Solver Tests

Tests for the surface energy balance that sets the face temperature in the
cell-centered finite volume discretization.

Tests cover:
1. Unit tests of the in-place Newton solve with each BC type, evaluated at the
   trial clamp of its iterate as the residual evaluates it
2. Energy balance verification
3. Flux, derivative and reported decomposition of each BC at a face
4. The face balances during a solve
"""

using Test
using Logging
using Pyrolysis
using Pyrolysis.Materials
using Pyrolysis.BoundaryConditions

using Pyrolysis.BoundaryConditions:
    solve_boundary_temperature!, BoundaryState, evaluate_bc_flux_and_derivative
using Pyrolysis.Materials: GaseousComponent, SolidComponent, Material, Reaction
using Pyrolysis.BoundaryConditions: compute_boundary_heat_fluxes
using Pyrolysis.Properties: mixture_emissivity, mixture_absorptivity
using Pyrolysis.States: trial_temperature

# Test helper: a tiny inert single-solid material (no gas component). Tests
# below pass it with a zero outward flux to recover the no-blowing limit,
# convenient for closed-form analytic verification.
const _BC_TEST_SOLID = SolidComponent(:test_solid; ρ = 1000.0, c = 1000.0, k = 0.2)
const _BC_TEST_MAT = Material(name = :bc_test_mat, components = (_BC_TEST_SOLID,))
const _BC_ZERO_J = ntuple(_ -> 0.0, Val(1))

# The face balance (no-blowing call by default), evaluated at the trial clamp
# of each iterate as the residual evaluates it; returns the diagnostic state.
function _solve_bt(
    T_cell,
    k,
    Δz,
    bc,
    t,
    ε;
    material = _BC_TEST_MAT,
    J = _BC_ZERO_J,
    kwargs...,
)
    state = BoundaryState()
    solve_boundary_temperature!(
        state,
        T_cell,
        k,
        Δz,
        bc,
        material,
        J,
        t,
        ε;
        trial = trial_temperature,
        kwargs...,
    )
    return state
end

# The face of the no-blowing call at surface temperature `T_s`.
_bc_face(T_s) = BlowingFace(_BC_ZERO_J, T_s, _BC_TEST_MAT)

@testset "Boundary Temperature Solver" begin

    # ========================================================================
    # Unit Tests for solve_boundary_temperature!()
    # ========================================================================
    @testset "solve_boundary_temperature!" begin

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

        @testset "Outward gas flux enters only through blowing" begin
            solid = SolidComponent(:test_solid; ρ = 1000.0, c = 1000.0, k = 0.2)
            gas = GaseousComponent(:test_gas; M = 0.044, c = 1500.0, k = 0.03, λ = 1e-5)
            material = Material(name = :gas_flux_test, components = (solid, gas))

            # A prescribed flux has no convective member, so the gas crossing
            # the face changes nothing: the surface balance carries no
            # gas-enthalpy term (the gas transport carries it in its face
            # enthalpy flux).
            bc = HeatFluxBC(30000.0)
            surface(J) = _solve_bt(700.0, 0.2, 50e-6, bc, 0.0, 0.9; material, J)
            no_flux = surface((0.0, 0.0))
            @test no_flux.converged
            for J in ((0.0, 0.04), (0.0, -0.04))
                st = surface(J)
                @test st.converged
                @test st.T_surface == no_flux.T_surface
                @test st.q_conducted == no_flux.q_conducted
            end
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
        # RadiativeFluxBC applies its absorbed flux at the face under the
        # surface models (NO_RADIATION, SURFACE_ABSORPTION) and none under
        # BEER_LAMBERT, whose in-depth radiation kernel absorbs it in the
        # cells.

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
            # according to Beer-Lambert law by the in-depth radiation kernel
            q, dq = evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff, BEER_LAMBERT)

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

            # The volumetric model returns zero
            q_bl, _ =
                evaluate_bc_flux_and_derivative(bc_with_α, T_s, t, ε_eff, BEER_LAMBERT)

            @test q_bl ≈ 0.0
        end

        @testset "Combined BC with RadiativeFluxBC" begin
            # When RadiativeFluxBC is combined with other BCs
            h_conv = 10.0
            T_inf = 300.0

            bc_combined =
                RadiativeFluxBC(I_incident) + ConvectiveBC(h = h_conv, T_inf = T_inf)
            face = _bc_face(T_s)

            # For NO_RADIATION: both components contribute
            q_no_rad, dq_no_rad = evaluate_bc_flux_and_derivative(
                bc_combined,
                T_s,
                t,
                ε_eff,
                NO_RADIATION,
                face,
            )
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
                face,
            )
            @test q_surf ≈ q_rad_expected + q_conv_expected
            @test dq_surf ≈ -h_conv

            # For BEER_LAMBERT: only convective contributes (radiation handled volumetrically)
            q_bl, dq_bl = evaluate_bc_flux_and_derivative(
                bc_combined,
                T_s,
                t,
                ε_eff,
                BEER_LAMBERT,
                face,
            )
            @test q_bl ≈ q_conv_expected
            @test dq_bl ≈ -h_conv
        end

        @testset "absorptivity = nothing absorbs at a_eff, emits at ε_eff" begin
            # `a_eff` is the material's surface-absorptivity mixture; every
            # overload takes it as a keyword defaulting to ε_eff (Kirchhoff).
            a_eff = 0.72
            q, dq = evaluate_bc_flux_and_derivative(
                bc,
                T_s,
                t,
                ε_eff,
                NO_RADIATION;
                a_eff = a_eff,
            )
            @test q ≈ I_incident * a_eff
            @test dq == 0.0
            q_default, _ =
                evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff, SURFACE_ABSORPTION)
            @test q_default ≈ I_incident * ε_eff

            # An explicit BC absorptivity ignores a_eff.
            bc_α = RadiativeFluxBC(flux = I_incident, absorptivity = 0.5)
            q_α, _ = evaluate_bc_flux_and_derivative(
                bc_α,
                T_s,
                t,
                ε_eff,
                NO_RADIATION;
                a_eff = a_eff,
            )
            @test q_α ≈ 0.5 * I_incident

            # Emission keeps ε_eff regardless of a_eff; combined BCs forward both.
            rad = RadiativeBC(T_inf = 300.0)
            q_emit, dq_emit = evaluate_bc_flux_and_derivative(
                rad,
                T_s,
                t,
                ε_eff,
                NO_RADIATION;
                a_eff = a_eff,
            )
            @test q_emit ≈ ε_eff * STEFAN_BOLTZMANN * (300.0^4 - T_s^4)
            mat_face = Material(
                name = :face,
                components = (
                    SolidComponent(:s; ρ = 1000.0, c = 1000.0, k = 0.2, ε = 0.9),
                    GaseousComponent(:g; M = 0.03, c = 1000.0, k = 0.03, λ = 1e-5),
                ),
                reactions = (),
            )
            q_c, dq_c = evaluate_bc_flux_and_derivative(
                bc + rad,
                T_s,
                t,
                ε_eff,
                NO_RADIATION,
                Pyrolysis.BoundaryConditions.BlowingFace((0.0, 0.0), T_s, mat_face);
                a_eff = a_eff,
            )
            @test q_c ≈ q + q_emit
            @test dq_c ≈ dq_emit

            # Volumetric models still contribute zero at the face.
            q_vol, _ = evaluate_bc_flux_and_derivative(
                bc,
                T_s,
                t,
                ε_eff,
                BEER_LAMBERT;
                a_eff = a_eff,
            )
            @test q_vol == 0.0
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
            q, dq = evaluate_bc_flux_and_derivative(
                bc,
                T_s,
                0.0,
                0.9,
                NO_RADIATION,
                _bc_face(T_s),
            )

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

        # An inert solid heated through the condition `top`, adiabatic below.
        function _heated_solid(top; tspan, output_times = nothing)
            solid = SolidComponent(:solid, ρ = 1000.0, c = 2000.0, k = 0.2, ε = 0.9)
            mat = Material(name = :TestSolid, components = (solid,), reactions = ())
            times = output_times === nothing ? (;) : (; output_times)
            return PyrolysisProblem(;
                mesh = create_uniform_mesh(0.005, 20),
                material = mat,
                bc_set = BoundaryConditionSet(
                    thermal = (top = top, bottom = AdiabaticBC()),
                ),
                T_initial = z -> 300.0,
                ξ_initial = [z -> 1000.0],
                tspan,
                times...,
            )
        end

        @testset "the face balances hold at both faces during a solve, the heated one converged" begin
            problem = _heated_solid(
                HeatFluxBC(50000.0) + RadiativeBC(ε = 0.9, T_inf = 300.0);
                tspan = (0.0, 1.0),
            )
            solution = Pyrolysis.Integration.solve(problem)

            @test solution.retcode == :Success

            # The boundary state the residual leaves on its kernel records a
            # converged balance at the heated face; the adiabatic face needs
            # no balance, its temperature is the cell's.
            def = solution.definition
            thermal = Pyrolysis.Kernels.find_kernel(
                Pyrolysis.Kernels.ThermalBoundaryKernel,
                def.kernels,
            )
            top = only(Pyrolysis.Meshes.boundary_faces(def.setup.topology, :top))
            @test thermal.states[top].converged
            @test solution.boundary_T[:bottom] == solution.T[1, :]

            # Surface temperature should be different from cell temperature
            T_surface = solution.surface_T[end]
            T_cell_top = solution.T[end, end]
            @test T_surface != T_cell_top

            # With 50 kW/m² heating, surface should be hotter than cell center
            @test T_surface > T_cell_top

            # Energy balance at the faces: the applied flux is the conducted one.
            @test solution.extras.surface_heat_flux_net[end] ≈
                  solution.extras.surface_heat_flux_conduction[end] rtol = 1e-6
            @test solution.extras.back_heat_flux[end] == 0.0
        end

        @testset "Solution output uses boundary temperature" begin
            # Verify that solution.surface_T uses T_surface, not T_cell
            problem = _heated_solid(
                HeatFluxBC(75000.0);
                tspan = (0.0, 5.0),
                output_times = [0.0, 5.0],
            )
            solution = Pyrolysis.Integration.solve(problem)

            @test solution.retcode == :Success

            # The reported surface temperature is the top face's balance, not
            # the top cell's temperature.
            T_surf_output = solution.surface_T[end]
            @test T_surf_output == solution.boundary_T[:top][end]

            # Surface should be noticeably hotter than top cell
            ΔT = T_surf_output - solution.T[end, end]
            @test ΔT > 5.0  # Expect meaningful difference with 75 kW/m²
        end

    end

end  # main testset

# Per-component surface absorptivity through the surface Newton solve. A single
# `RadiativeFluxBC(flux = q)` with `absorptivity = nothing` absorbs a_eff·q and
# emits at ε_eff. When every component's absorptivity lies on one line with its
# emissivity, a_j = a + b·ε_j, the mixture is a_eff = a + b·ε_eff and the single
# BC is interchangeable with the two-piece split `RadiativeFluxBC(b·q)`
# (Kirchhoff) + `RadiativeFluxBC(a·q, absorptivity = 1)`.
@testset "surface Newton: a_eff absorption vs. two-piece split" begin
    a_line, b_line = 0.30, 0.55
    ε_v, ε_c = 0.90, 0.80
    gas = GaseousComponent(:g, M = 0.03, c = 1000.0, k = 0.03, λ = 1e-5)
    virgin = SolidComponent(
        :v,
        ρ = 500.0,
        c = 1500.0,
        k = 0.15,
        ε = ε_v,
        absorptivity = a_line + b_line * ε_v,
    )
    char = SolidComponent(
        :c,
        ρ = 150.0,
        c = 1700.0,
        k = 0.08,
        ε = ε_c,
        absorptivity = a_line + b_line * ε_c,
    )
    mat = Material(name = :line, components = (virgin, char, gas), reactions = ())
    virgin_k = SolidComponent(:v, ρ = 500.0, c = 1500.0, k = 0.15, ε = ε_v)
    char_k = SolidComponent(:c, ρ = 150.0, c = 1700.0, k = 0.08, ε = ε_c)
    mat_k =
        Material(name = :kirchhoff, components = (virgin_k, char_k, gas), reactions = ())

    ξ = (300.0, 40.0, 0.5)
    T_cell = 650.0
    ε_eff = mixture_emissivity(ξ, T_cell, mat)
    a_eff = mixture_absorptivity(ξ, T_cell, mat)
    @test ε_eff ≈ mixture_emissivity(ξ, T_cell, mat_k)
    @test a_eff ≈ a_line + b_line * ε_eff rtol = 1e-14
    @test mixture_absorptivity(ξ, T_cell, mat_k) === mixture_emissivity(ξ, T_cell, mat_k)

    @testset "effective_absorptivity on the line" begin
        effective_absorptivity = Pyrolysis.Properties.effective_absorptivity
        effective_emissivity = Pyrolysis.Properties.effective_emissivity
        @test effective_emissivity(ξ, T_cell, mat) ≈ ε_eff rtol = 1e-14
        @test effective_absorptivity(ξ, T_cell, mat) ≈ a_eff rtol = 1e-14
        @test effective_absorptivity(ξ, T_cell, mat_k) ===
              effective_emissivity(ξ, T_cell, mat_k)
    end

    q_inc = 50e3
    k_eff, Δz_half, t = 0.12, 2.5e-4, 0.0
    T_inf, h = 300.0, 10.0
    J0 = ntuple(_ -> 0.0, Val(3))
    single =
        RadiativeFluxBC(flux = q_inc) +
        RadiativeBC(T_inf = T_inf) +
        ConvectiveBC(h = h, T_inf = T_inf)
    split =
        RadiativeFluxBC(flux = b_line * q_inc) +
        RadiativeFluxBC(flux = a_line * q_inc, absorptivity = 1.0) +
        RadiativeBC(T_inf = T_inf) +
        ConvectiveBC(h = h, T_inf = T_inf)

    s_single =
        _solve_bt(T_cell, k_eff, Δz_half, single, t, ε_eff; material = mat, J = J0, a_eff)
    s_split = _solve_bt(T_cell, k_eff, Δz_half, split, t, ε_eff; material = mat_k, J = J0)
    @test s_single.converged && s_split.converged
    @test s_single.T_surface ≈ s_split.T_surface rtol = 1e-9
    @test s_single.q_applied ≈ s_split.q_applied rtol = 1e-9

    # The converged balance decomposes as absorbed(a_eff) + emitted(ε_eff) + convective.
    T_s = s_single.T_surface
    q_expected =
        a_eff * q_inc + ε_eff * STEFAN_BOLTZMANN * (T_inf^4 - T_s^4) + h * (T_inf - T_s)
    @test s_single.q_applied ≈ q_expected rtol = 1e-12
    @test s_single.q_conducted ≈ k_eff * (T_s - T_cell) / Δz_half rtol = 1e-9

    # Default a_eff = ε_eff: the Kirchhoff material absorbs ε_eff·q instead.
    s_k = _solve_bt(T_cell, k_eff, Δz_half, single, t, ε_eff; material = mat_k, J = J0)
    T_k = s_k.T_surface
    @test s_k.q_applied ≈
          ε_eff * q_inc + ε_eff * STEFAN_BOLTZMANN * (T_inf^4 - T_k^4) + h * (T_inf - T_k) rtol =
        1e-12
    @test T_k != T_s

    # The diagnostic decomposition applies the same a_eff.
    hf_face = compute_boundary_heat_fluxes(
        single,
        T_s,
        t,
        ε_eff,
        SURFACE_ABSORPTION,
        Pyrolysis.BoundaryConditions.BlowingFace((0.0, 0.0, 0.0), T_s, mat);
        a_eff = a_eff,
    )
    @test hf_face.absorbed ≈ a_eff * q_inc
    @test hf_face.net ≈ s_single.q_applied rtol = 1e-12

    @testset "compute_boundary_heat_fluxes(::CombinedThermalBC) without a face" begin
        hf = compute_boundary_heat_fluxes(
            single,
            T_s,
            t,
            ε_eff,
            SURFACE_ABSORPTION;
            a_eff = a_eff,
        )
        @test hf.absorbed ≈ a_eff * q_inc
        @test hf.net ≈ s_single.q_applied rtol = 1e-12
        @test hf_face.net ≈ hf.net
    end
end

# The surface Newton starts from the cell temperature, never from the stored
# BoundaryState, and stores only a finite iterate: a poisoned state cannot leak
# into a later solve.
@testset "surface Newton ignores the stored state" begin
    solve! = Pyrolysis.BoundaryConditions.solve_boundary_temperature!
    bc =
        HeatFluxBC(50e3) +
        RadiativeBC(T_inf = 300.0) +
        ConvectiveBC(h = 10.0, T_inf = 300.0)
    T_cell, k_eff, Δz_half, ε = 650.0, 0.2, 1e-3, 0.9
    fresh = _solve_bt(T_cell, k_eff, Δz_half, bc, 0.0, ε)
    @test fresh.converged

    poisoned = BoundaryState(NaN, NaN, NaN, NaN, NaN, NaN, false, 0)
    result = solve!(
        poisoned,
        T_cell,
        k_eff,
        Δz_half,
        bc,
        _BC_TEST_MAT,
        _BC_ZERO_J,
        0.0,
        ε;
        trial = trial_temperature,
    )
    @test result.T_s == fresh.T_surface
    @test poisoned.T_surface == fresh.T_surface
    @test poisoned.q_conducted == fresh.q_conducted
    @test poisoned.converged

    # A non-finite cell temperature gives a non-finite result, which is returned
    # and not stored.
    bad = @test_logs min_level = Logging.Warn solve!(
        poisoned,
        NaN,
        k_eff,
        Δz_half,
        bc,
        _BC_TEST_MAT,
        _BC_ZERO_J,
        0.0,
        ε;
        trial = trial_temperature,
    )
    @test isnan(bad.T_s)
    @test poisoned.T_surface == fresh.T_surface
    @test poisoned.converged

    @testset "solve_boundary_temperature! without a trial clamp: a NaN cell warns" begin
        bad = @test_logs (:warn, r"did not converge") match_mode = :any solve!(
            poisoned,
            NaN,
            k_eff,
            Δz_half,
            bc,
            _BC_TEST_MAT,
            _BC_ZERO_J,
            0.0,
            ε,
        )
        @test isnan(bad.T_s)
        @test poisoned.T_surface == fresh.T_surface
        @test poisoned.converged
    end
end

# An exhausted iteration budget is recorded on the state; the trial-clamped
# solve reports it there, not in the log.
@testset "surface Newton: an exhausted iteration budget is recorded" begin
    bc = RadiativeBC(T_inf = 300.0) + ConvectiveBC(h = 10.0, T_inf = 300.0)
    # One iteration cannot converge the quartic radiation balance.
    state = @test_logs min_level = Logging.Warn _solve_bt(
        800.0,
        0.15,
        1e-3,
        bc,
        0.0,
        0.9;
        max_iter = 1,
    )
    @test state.converged == false
    @test state.iterations == 1
    @test abs(state.q_applied - state.q_conducted) > 1.0
    # The default budget converges.
    state_ok = _solve_bt(800.0, 0.15, 1e-3, bc, 0.0, 0.9)
    @test state_ok.converged
    @test state_ok.q_applied ≈ state_ok.q_conducted rtol = 1e-8
end

# The reported decomposition's `net` is the flux the surface balance applies,
# for a combined and a standalone condition under every radiation model.
@testset "net is the applied face flux, per BC and per model" begin
    combined =
        RadiativeFluxBC(flux = 50e3) +
        RadiativeBC(T_inf = 300.0) +
        ConvectiveBC(h = 10.0, T_inf = 300.0)
    standalone = RadiativeFluxBC(flux = 50e3)
    ε_eff = 0.9
    T_s = 700.0
    t = 0.0
    face = _bc_face(T_s)

    for bc in (combined, standalone),
        model in (NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT)

        hf = compute_boundary_heat_fluxes(bc, T_s, t, ε_eff, model, face)
        q_applied, _ = evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff, model, face)
        @test hf.net ≈ q_applied atol = 1e-9
    end

    # The two most regime-sensitive cases, pinned explicitly:
    # (a) standalone RadiativeFluxBC under SURFACE_ABSORPTION reports the
    #     ε_eff·q_ext the Newton applies.
    hf = compute_boundary_heat_fluxes(standalone, T_s, t, ε_eff, SURFACE_ABSORPTION, face)
    @test hf.net ≈ ε_eff * 50e3
    @test hf.absorbed ≈ ε_eff * 50e3

    # (b) combined under BEER_LAMBERT: net carries only the rerad +
    #     convection the face actually applies, not α·q_ext, and the face
    #     absorbs nothing (the deposit is volumetric).
    hf_bl = compute_boundary_heat_fluxes(combined, T_s, t, ε_eff, BEER_LAMBERT, face)
    q_rerad = ε_eff * STEFAN_BOLTZMANN * (300.0^4 - T_s^4)
    q_conv = 10.0 * (300.0 - T_s)
    @test hf_bl.absorbed == 0.0
    @test hf_bl.net ≈ q_rerad + q_conv rtol = 1e-12
    @test abs(
        hf_bl.net -
        (hf_bl.absorbed - hf_bl.reradiation - hf_bl.convection + hf_bl.prescribed),
    ) < 1e-9  # components stay consistent with net

    # (c) combined under SURFACE_ABSORPTION includes the absorbed flux.
    hf_sa = compute_boundary_heat_fluxes(combined, T_s, t, ε_eff, SURFACE_ABSORPTION, face)
    @test hf_sa.net ≈ ε_eff * 50e3 + q_rerad + q_conv rtol = 1e-12
end

# The P1 radiation model takes the incident flux in its own solver.
@testset "P1_QUASI_STEADY: RadiativeFluxBC contributes nothing at the face" begin
    T_s, t, ε_eff, I_incident = 500.0, 0.0, 0.9, 50e3
    q, dq = evaluate_bc_flux_and_derivative(
        RadiativeFluxBC(I_incident),
        T_s,
        t,
        ε_eff,
        Pyrolysis.Materials.P1_QUASI_STEADY,
    )
    @test q ≈ 0.0   # No boundary flux contribution
    @test dq ≈ 0.0  # No temperature dependence
    bc_with_α = RadiativeFluxBC(flux = I_incident, absorptivity = 0.7)
    q_p1, _ = evaluate_bc_flux_and_derivative(
        bc_with_α,
        T_s,
        t,
        ε_eff,
        Pyrolysis.Materials.P1_QUASI_STEADY,
    )
    @test q_p1 ≈ 0.0
end

# The face-less forms of the combined condition and of the convective one.
@testset "evaluate_bc_flux_and_derivative without a face: combined and convective" begin
    T_s, t, ε_eff, I_incident = 500.0, 0.0, 0.9, 50e3
    h_conv, T_inf = 10.0, 300.0
    bc_combined = RadiativeFluxBC(I_incident) + ConvectiveBC(h = h_conv, T_inf = T_inf)
    q_rad_expected = I_incident * ε_eff
    q_conv_expected = h_conv * (T_inf - T_s)
    q_no_rad, dq_no_rad =
        evaluate_bc_flux_and_derivative(bc_combined, T_s, t, ε_eff, NO_RADIATION)
    @test q_no_rad ≈ q_rad_expected + q_conv_expected
    @test dq_no_rad ≈ -h_conv
    q_surf, dq_surf =
        evaluate_bc_flux_and_derivative(bc_combined, T_s, t, ε_eff, SURFACE_ABSORPTION)
    @test q_surf ≈ q_rad_expected + q_conv_expected
    @test dq_surf ≈ -h_conv
    q_bl, dq_bl = evaluate_bc_flux_and_derivative(bc_combined, T_s, t, ε_eff, BEER_LAMBERT)
    @test q_bl ≈ q_conv_expected
    @test dq_bl ≈ -h_conv

    # a_eff is forwarded to the members, and the face form with no outflow
    # agrees.
    a_eff = 0.72
    bc = RadiativeFluxBC(I_incident)
    rad = RadiativeBC(T_inf = 300.0)
    q, _ = evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff, NO_RADIATION; a_eff = a_eff)
    q_emit, dq_emit =
        evaluate_bc_flux_and_derivative(rad, T_s, t, ε_eff, NO_RADIATION; a_eff = a_eff)
    q_c, dq_c = evaluate_bc_flux_and_derivative(
        bc + rad,
        T_s,
        t,
        ε_eff,
        NO_RADIATION;
        a_eff = a_eff,
    )
    @test q_c ≈ q + q_emit
    @test dq_c ≈ dq_emit
    mat_face = Material(
        name = :face,
        components = (
            SolidComponent(:s; ρ = 1000.0, c = 1000.0, k = 0.2, ε = 0.9),
            GaseousComponent(:g; M = 0.03, c = 1000.0, k = 0.03, λ = 1e-5),
        ),
        reactions = (),
    )
    q6, _ = evaluate_bc_flux_and_derivative(
        bc + rad,
        T_s,
        t,
        ε_eff,
        NO_RADIATION,
        Pyrolysis.BoundaryConditions.BlowingFace((0.0, 0.0), T_s, mat_face);
        a_eff = a_eff,
    )
    @test q6 ≈ q_c

    # The convective condition's own face-less form.
    q_cv, dq_cv = evaluate_bc_flux_and_derivative(
        ConvectiveBC(h = 10.0, T_inf = 300.0),
        T_s,
        0.0,
        0.9,
    )
    @test q_cv ≈ 10.0 * (300.0 - T_s)
    @test dq_cv ≈ -10.0
end

# The mesh-bound accessors of the stored boundary states.
@testset "Discretization accessors: get_surface_temperature, check_boundary_energy_balance" begin
    get_surface_temperature = Pyrolysis.Discretization.get_surface_temperature
    get_back_face_temperature = Pyrolysis.Discretization.get_back_face_temperature
    check_boundary_energy_balance = Pyrolysis.Discretization.check_boundary_energy_balance


    @testset "get_surface_temperature and get_back_face_temperature" begin
        # Create a mesh with boundary states initialized
        L = 0.01
        n_cells = 10
        mesh = create_uniform_mesh(L, n_cells, Val(2))

        # Manually set boundary states for testing
        T_top = 800.0
        T_bottom = 350.0
        mesh.boundary_states[1] = Pyrolysis.BoundaryConditions.BoundaryState(
            T_top,
            50000.0,
            50000.0,
            750.0,
            0.2,
            50e-6,
            true,
            3,
        )
        mesh.boundary_states[2] = Pyrolysis.BoundaryConditions.BoundaryState(
            T_bottom,
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
        mesh_test.boundary_states[1] = Pyrolysis.BoundaryConditions.BoundaryState(
            600.0,    # T_surface
            10000.0,  # q_applied
            10000.0,  # q_conducted (matching!)
            550.0,    # T_cell
            0.2,      # k_eff
            50e-6,    # Δz_half
            true,     # converged
            3,        # iterations
        )
        mesh_test.boundary_states[2] = Pyrolysis.BoundaryConditions.BoundaryState(
            350.0,
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
        mesh_test.boundary_states[1] = Pyrolysis.BoundaryConditions.BoundaryState(
            600.0,    # T_surface
            10000.0,  # q_applied
            8000.0,   # q_conducted (20% error!)
            550.0,    # T_cell
            0.2,      # k_eff
            50e-6,    # Δz_half
            true,     # converged
            3,        # iterations
        )

        balance_results_bad = check_boundary_energy_balance(mesh_test)
        # 20% imbalance
        @test balance_results_bad[1] ≈ 0.2 atol=0.01
    end

end
