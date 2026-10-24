"""
Material-supplied gas viscosity: the `SutherlandProperty` type and the
`Material.gas_viscosity` field.

- `SutherlandProperty()` carries the canonical air triple and evaluates
  Sutherland's law over the physical range.
- `Material` defaults to air.
- The Darcy mobility reads the curve at the trial temperature and floors it
  at `MU_FLOOR`: `μ` appears only as the denominator of the mobility `κ/μ`,
  so a user curve that goes non-positive outside its fit range must not reach
  the solver.
"""

using Test
using ForwardDiff
using Pyrolysis
using Pyrolysis.Materials:
    Material,
    SolidComponent,
    GaseousComponent,
    Reaction,
    SutherlandProperty,
    PolynomialProperty,
    ConstantProperty
using Pyrolysis.Materials: T_PROPERTY_MIN
using Pyrolysis.Properties: MU_FLOOR
using Pyrolysis.Properties: skeletal_volume_fraction, effective_permeability
using Pyrolysis.Physics: ideal_gas_pressure, background_gas_concentration
using Pyrolysis.Materials: DepletionLimiter
using Pyrolysis.BoundaryConditions:
    BoundaryConditionSet, DirichletThermalBC, DirichletPressureBC, ImpermeableBC

# Minimal Darcy-capable material; `gas_viscosity` overridable per case.
function _visc_material(; kwargs...)
    virgin = SolidComponent(:wood, ρ = 500.0, c = 1500.0, k = 0.15, λ = 1e-5, κ = 1e-14)
    gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5)
    return Material(;
        name = :visc_test,
        components = (virgin, gas),
        reactions = (Reaction(:pyrolysis, 1 => 2, A = 1e10, E = 150e3, h = -300e3),),
        kwargs...,
    )
end

# The derived state of a Darcy column at a uniform temperature `T`.
function _visc_derived(mat, T)
    problem = PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, 4),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (top = DirichletPressureBC(101325.0), bottom = ImpermeableBC()),
        ),
        ξ_initial = [z -> 450.0, z -> 0.5],
        tspan = (0.0, 1.0),
    )
    cp = Pyrolysis.Integration.setup(problem)
    layout = cp.def.layout
    u = copy(cp.u0)
    for i = 1:4
        u[Pyrolysis.States.energy_index(layout, i)] = T
    end
    du = zeros(length(u))
    Pyrolysis.Assembly.residual!(du, u, cp.def, 0.0, cp.ws)
    @assert cp.def.setup.gas isa Pyrolysis.States.DarcyFickTransport
    ds, _ = Pyrolysis.Assembly.buffers(cp.ws, cp.def, Float64)
    return ds
end

@testset "Material gas viscosity" begin

    # ========================================================================
    @testset "SutherlandProperty is Sutherland's law" begin
        p = SutherlandProperty()
        for T in (273.15, 300.0, 500.0, 800.0, 1200.0, 2000.0)
            @test p(T) ≈ 1.716e-5 * (T / 273.15)^1.5 * (273.15 + 110.4) / (T + 110.4) rtol =
                1e-14
        end
        # Documented anchor values.
        @test p(273.15) ≈ 1.716e-5 rtol = 1e-10
        @test p(300.0) ≈ 1.846e-5 rtol = 1e-3
    end

    # ========================================================================
    @testset "property-function contract" begin
        p = SutherlandProperty()
        ξ = (500.0, 1.0)
        @test p(650.0, ξ) == p(650.0)               # 2-arg form
        @test p(ForwardDiff.Dual(650.0, 1.0)) isa ForwardDiff.Dual  # AD-clean
    end

    # ========================================================================
    @testset "Material defaults to air" begin
        mat = _visc_material()
        @test mat.gas_viscosity isa SutherlandProperty
        for T in (300.0, 650.0, 1100.0)
            # The default path is bit-for-bit the air Sutherland curve, which
            # every calibrated case relies on.
            @test mat.gas_viscosity(T) == SutherlandProperty()(T)
        end
    end

    # ========================================================================
    @testset "material curve is actually consulted" begin
        mat = _visc_material(gas_viscosity = ConstantProperty(3.0e-5))
        for T in (300.0, 900.0)
            @test mat.gas_viscosity(T) == 3.0e-5
            @test mat.gas_viscosity(T) != SutherlandProperty()(T)
        end

        # A bare Function and a Real both go through `to_property`.
        @test _visc_material(gas_viscosity = 2.5e-5).gas_viscosity(400.0) == 2.5e-5
        mat_f = _visc_material(gas_viscosity = T -> 1.0e-5 + 1.0e-8 * T)
        @test mat_f.gas_viscosity(400.0) ≈ 1.4e-5 rtol = 1e-12
    end

    # ========================================================================
    @testset "polynomial curve (the surrogate's native form)" begin
        # Ascending coefficients, as `evalpoly` and the Cantera surrogate agree.
        coeffs = (
            -2.590750377731154e-7,
            3.664256350094585e-8,
            -6.54432524120188e-14,
            -6.352113427713809e-15,
            2.3401486488384493e-18,
        )
        # Constructing at all is the no-false-rejection gate: a real fitted
        # correlation must survive `_validate_gas_viscosity`'s probe range.
        mat = _visc_material(gas_viscosity = PolynomialProperty(coeffs))
        # Reproduces the surrogate's reported Cantera blend values. These
        # coefficients are a quartic FIT to those values (max residual ~5e-4
        # over 280-1200 K), so the tolerance is the fit residual, not roundoff.
        @test mat.gas_viscosity(300.0) ≈ 1.0579843085895027e-5 rtol = 1e-3
        @test mat.gas_viscosity(1000.0) ≈ 3.2307444332042786e-5 rtol = 1e-3
        # Below air at every pyrolysis-relevant temperature (the physical point).
        for T in (300.0, 500.0, 700.0, 1000.0)
            @test mat.gas_viscosity(T) < SutherlandProperty()(T)
        end
    end

    # ========================================================================
    @testset "guards: μ only ever divides" begin
        # Newton trial states can present nonphysical T; the curve is read on
        # the property bracket instead of throwing DomainError on T^1.5.
        mat = _visc_material()
        @test mat.gas_viscosity(-50.0) == mat.gas_viscosity(T_PROPERTY_MIN)
        @test isfinite(mat.gas_viscosity(0.0))

        # The Darcy mobility is κ_eff/μ with μ floored at MU_FLOOR. A user
        # curve that goes non-positive must be floored, not passed through to
        # `κ/μ` where it would flip the sign of the Darcy velocity. This one
        # stays positive across the constructor's probe range (≤1800 K) and
        # only turns negative beyond it — the realistic failure, a fit used
        # outside its validity range by a Newton trial state that overshoots.
        neg = _visc_material(gas_viscosity = T -> 1.0e-5 - 4.0e-9 * T)
        for (m, T, μ) in (
            (mat, 650.0, mat.gas_viscosity(650.0)),
            (neg, 1000.0, neg.gas_viscosity(1000.0)),   # still physical
            (neg, 3000.0, MU_FLOOR),                    # floored
        )
            ds = _visc_derived(m, T)
            for i in eachindex(ds.mob)
                ξ = (ds.ξ⁺[1, i], ds.ξ⁺[2, i])
                κ = effective_permeability(ξ, T, m, ds.φ[i])
                @test κ > 0
                @test ds.mob[i] ≈ κ / μ rtol = 1e-12
            end
        end
        @test neg.gas_viscosity(3000.0) < 0.0

        # The floor sits far below any physical gas viscosity, so it never
        # binds on a legitimate state.
        for T in (200.0, 300.0, 650.0, 1500.0)
            @test mat.gas_viscosity(T) > 1e6 * MU_FLOOR
        end
    end

    # ========================================================================
    # End-to-end: the field must reach the solved pressure field, not merely
    # the property function. V-GT-02's exact solution is LINEAR in μ,
    #
    #   P²(z) = P∞² + (μ·R·T·ġ/(κ·M))·(L² − z²),
    #
    # so doubling the material's viscosity must double the steady overpressure².
    # If `gas_viscosity` were ignored anywhere between `Material` and the Darcy
    # mobility, the two solves would coincide and this fails hard.
    @testset "end-to-end: doubling μ doubles the Darcy overpressure²" begin
        κ_perm, L, T_iso, M_gas, P∞, ġ = 1e-12, 0.01, 500.0, 0.028, 101325.0, 5.0
        R = 8.314462618

        function solve_with(μ_curve)
            solid = SolidComponent(
                :solid,
                ρ = 1000.0,
                c = 1000.0,
                k = 0.5,
                κ = κ_perm,
                λ = 0.0,
                ρ_intrinsic = 1e5,
            )
            gas = GaseousComponent(:gas, M = M_gas, c = 1000.0, k = 0.03, λ = 0.0)
            rxn = Reaction(:generate, 1 => 2; A = ġ, E = 0.0, h = 0.0, n = 0.0)
            mat = Material(
                name = :VGT02_visc,
                components = (solid, gas),
                reactions = (rxn,),
                depletion_limiter = DepletionLimiter(enabled = false),
                gas_viscosity = μ_curve,
            )
            ξ_g0 = background_gas_concentration(mat, 2, (1000.0, 0.0), T_iso; P = P∞)
            mesh = create_uniform_mesh(L, 40)
            bc = BoundaryConditionSet(
                thermal = (
                    top = DirichletThermalBC(T_iso),
                    bottom = DirichletThermalBC(T_iso),
                ),
                mass = (top = DirichletPressureBC(P∞), bottom = ImpermeableBC()),
            )
            prob = PyrolysisProblem(
                mesh = mesh,
                material = mat,
                bc_set = bc,
                T_initial = z -> T_iso,
                ξ_initial = [z -> 1000.0, z -> ξ_g0],
                tspan = (0.0, 1.0),
                output_times = [0.0, 1.0],
            )
            sol = Pyrolysis.Integration.solve(prob; abstol = 1e-10, reltol = 1e-8)
            @test sol.retcode == :Success
            z = Pyrolysis.Outputs.get_cell_positions(sol)
            P = map(eachindex(z)) do i
                ξi = (sol.ξ[1, i, end], sol.ξ[2, i, end])
                φi = 1 - skeletal_volume_fraction(ξi, T_iso, mat)
                ideal_gas_pressure(ξi, T_iso, φi, mat)
            end
            return z, P, mat
        end

        air = SutherlandProperty()
        doubled = SutherlandProperty(2 * air.μ_ref, air.T_ref, air.S)   # 2·μ(T) ∀T

        z1, P1, mat1 = solve_with(air)
        z2, P2, mat2 = solve_with(doubled)

        # Each solve matches its OWN exact profile (the case's usual accuracy).
        for (z, P, mat) in ((z1, P1, mat1), (z2, P2, mat2))
            μ = mat.gas_viscosity(T_iso)
            coeff = μ * R * T_iso * ġ / (κ_perm * M_gas)
            err = maximum(eachindex(z)) do i
                abs(P[i]^2 - (P∞^2 + coeff * (L^2 - z[i]^2))) / (coeff * L^2)
            end
            @test err < 5e-3
        end

        # And the overpressure² ratio is 2 — the signature of μ actually
        # reaching the mobility κ/μ.
        over1 = P1 .^ 2 .- P∞^2
        over2 = P2 .^ 2 .- P∞^2
        @test all(over1 .> 0)
        @test maximum(abs, over2 ./ over1 .- 2.0) < 1e-2
    end

    # ========================================================================
    @testset "constructor rejects an unusable curve" begin
        # Caught at construction, not discovered as a sign-flipped Darcy
        # velocity three hours into a solve.
        @test_throws ArgumentError _visc_material(gas_viscosity = ConstantProperty(-1.0e-5))
        @test_throws ArgumentError _visc_material(gas_viscosity = ConstantProperty(0.0))
        @test_throws ArgumentError _visc_material(gas_viscosity = T -> NaN)
    end
end
