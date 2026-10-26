"""
Property laws are evaluated on the bracket [T_PROPERTY_MIN, T_PROPERTY_MAX]
(argument clamped) and validated on it at construction: ρ, c, k positive,
ε within [0, 1]. The error names the first offending temperature and the
opt-out; inside the bracket the clamp is the identity with derivative 1.
"""

using Test
using ForwardDiff
using Pyrolysis
using Pyrolysis.Materials:
    T_PROPERTY_MIN, T_PROPERTY_MAX, property_temperature, identically_zero

@testset "property-law bracket" begin
    @testset "argument clamp is the identity inside the bracket" begin
        law = CallableProperty(T -> 0.2156 * (T / 300.0)^-0.745)   # throws for T < 0 unclamped
        for T in (T_PROPERTY_MIN, 300.0, 714.29, 2999.5, T_PROPERTY_MAX)
            @test law(T) === 0.2156 * (T / 300.0)^-0.745
        end
        @test law(-8.84) == law(T_PROPERTY_MIN)
        @test law(1e5) == law(T_PROPERTY_MAX)
        @test ForwardDiff.derivative(law, 500.0) ≈
              -0.745 * 0.2156 / 300.0 * (500.0 / 300.0)^-1.745
        @test ForwardDiff.derivative(law, -8.84) == 0.0
        @test property_temperature(500.0) === 500.0
        # Every temperature-dependent law type takes the clamp.
        @test LinearProperty(1.0, 1.0)(-100.0) == 1.0 + T_PROPERTY_MIN
        @test PolynomialProperty((1.0, 1.0))(1e6) == 1.0 + T_PROPERTY_MAX
        @test ArrheniusProperty(1.0, 1e4)(-5.0) ==
              ArrheniusProperty(1.0, 1e4)(T_PROPERTY_MIN)
        @test SutherlandProperty()(-5.0) == SutherlandProperty()(T_PROPERTY_MIN)
        @test ConstantProperty(3.0)(-5.0) == 3.0
    end

    @testset "zero-crossing conductivity is rejected with the crossing temperature" begin
        k_linear = T -> T <= 394.0 ? 0.16 : (0.3 - 4.2e-4 * T)   # zero at 714.29 K
        err = try
            SolidComponent(:PMMA, ρ = 1210.0, c = 1500.0, k = k_linear)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("thermal conductivity", err.msg)
        @test occursin("T = 715.0 K", err.msg)        # first 1 K grid point past the crossing
        @test occursin("validate_properties = false", err.msg)
        comp = SolidComponent(
            :PMMA,
            ρ = 1210.0,
            c = 1500.0,
            k = k_linear,
            validate_properties = false,
        )
        @test comp.thermal_conductivity(800.0) < 0
    end

    @testset "other laws and component types" begin
        @test_throws ArgumentError SolidComponent(
            :s,
            ρ = 1000.0,
            c = T -> -2290.0 + 11.2 * T,
            k = 0.2,
        )
        @test_throws ArgumentError SolidComponent(
            :s,
            ρ = T -> 1000.0 - 0.5 * T,
            c = 1000.0,
            k = 0.2,
        )
        @test_throws ArgumentError SolidComponent(
            :s,
            ρ = 1000.0,
            c = 1000.0,
            k = 0.2,
            ε = (1.2, 0.0),
        )
        @test_throws ArgumentError SolidComponent(
            :s,
            ρ = 1000.0,
            c = 1000.0,
            k = 0.2,
            ε = (0.5, -3e-4),
        )
        @test_throws ArgumentError LiquidComponent(
            :w,
            ρ = 1000.0,
            c = 4180.0,
            k = T -> 0.6 - 3e-4 * T,
        )
        @test_throws ArgumentError GaseousComponent(
            :g,
            M = 0.03,
            c = T -> 100.0 - 0.1 * T,
            k = 0.03,
        )
        @test_throws ArgumentError GaseousComponent(
            :g,
            M = 0.03,
            c = 1000.0,
            k = 0.03,
            ε = 1.5,
        )
        # Laws positive on the whole bracket build normally.
        s = SolidComponent(
            :ok,
            ρ = 500.0,
            c = T -> 1100.0 + 0.6 * T,
            k = T -> 0.1 + 6e-4 * T,
            ε = (0.95, -1.8e-4),
        )
        @test s.thermal_conductivity(T_PROPERTY_MAX) > 0
        g = GaseousComponent(:g, M = 0.028, c = T -> 950.0 + 0.3 * T, k = 0.03)
        @test g.heat_capacity(T_PROPERTY_MIN) > 0
        w = LiquidComponent(:w, ρ = 1000.0, c = 4180.0, k = 0.6)
        @test w.density(300.0) == 1000.0
    end

    @testset "state-dependent laws are not sampled" begin
        h = StateDependentProperty((T, ξ) -> ξ[1] > 0 ? 1.0 : -1.0)
        @test SolidComponent(:s, ρ = 1000.0, c = 1000.0, k = 0.2, ε = h) isa SolidComponent
    end
end

@testset "a law vanishes by its structure, never by samples" begin
    @test identically_zero(ConstantProperty(0.0))
    @test !identically_zero(ConstantProperty(1e-5))
    @test identically_zero(LinearProperty(0.0, 0.0))
    @test !identically_zero(LinearProperty(0.0, 1e-8))
    @test identically_zero(PolynomialProperty((0.0, 0.0, 0.0)))
    @test !identically_zero(PolynomialProperty((0.0, 0.0, 1e-12)))
    @test identically_zero(ArrheniusProperty(0.0, 5e4, 300.0))
    @test identically_zero(TableProperty((300.0, 600.0), (0.0, 0.0)))
    # Zero at 200, 300, 1000 and 3000 K and positive between them.
    @test !identically_zero(
        TableProperty((200.0, 300.0, 600.0, 1000.0, 3000.0), (0.0, 0.0, 1e-5, 0.0, 0.0)),
    )
    @test !identically_zero(CallableProperty(T -> 0.0))
end

# Boundary-condition parameters: a number or a function of time, and of
# temperature for the emissivity-like ones; a number takes the argument's type.
@testset "parameter evaluation: numbers, functions and property laws" begin
    evaluate_param = Pyrolysis.Materials.evaluate_param
    evaluate_param_T = Pyrolysis.Materials.evaluate_param_T
    @test evaluate_param(25e3, 0.0) === 25e3
    @test evaluate_param(300.0, 10.0) == 300.0
    @test evaluate_param(t -> 100.0 * t, 5.0) == 500.0
    @test evaluate_param_T(0.9, 500.0) === 0.9
    @test evaluate_param_T(T -> 0.8 + 1e-4 * T, 500.0) ≈ 0.85
    @test evaluate_param_T(LinearProperty(0.8, 1e-4), 500.0) ==
          LinearProperty(0.8, 1e-4)(500.0)
    @test evaluate_param_T(ConstantProperty(0.7), 500.0) == 0.7
    # A Dual argument: the number is a constant of it.
    d = ForwardDiff.Dual(5.0, 1.0)
    @test evaluate_param(25e3, d) isa ForwardDiff.Dual
    @test ForwardDiff.partials(evaluate_param(25e3, d))[1] == 0.0
    @test evaluate_param_T(0.9, d) isa ForwardDiff.Dual
    @test ForwardDiff.partials(evaluate_param_T(0.9, d))[1] == 0.0
    @test ForwardDiff.derivative(t -> evaluate_param(t -> 100.0 * t, t), 5.0) == 100.0
    # Every method infers its result type.
    @test @inferred(evaluate_param(100.0, 1.0)) isa Float64
    @test @inferred(evaluate_param(t -> 100.0 * t, 1.0)) isa Float64
    @test @inferred(evaluate_param_T(0.9, 500.0)) isa Float64
    @test @inferred(evaluate_param_T(T -> 0.8 + 1e-4 * T, 500.0)) isa Float64
    @test @inferred(evaluate_param_T(LinearProperty(0.8, 1e-4), 500.0)) isa Float64
    @test @inferred(evaluate_param_T(ConstantProperty(0.7), 500.0)) isa Float64
end
