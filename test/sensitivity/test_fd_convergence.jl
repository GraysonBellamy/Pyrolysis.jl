using Test

include(joinpath(@__DIR__, "..", "gates", "fd_convergence.jl"))

# The finite-difference convergence rule of `gates/fd_convergence.jl` on the
# synthetic objective J(ε) = exp(a ε), whose derivative at ε = 0 is `a` and
# whose central difference is sinh(a h)/h = a (1 + x²/6 + x⁴/120 + …) with
# x = a h, so that Richardson's extrapolant of h and 2h is a (1 − x⁴/30 + …);
# plus a noise `η(ε)` of known values: the values the gate's runs give at
# ε = 0 and ±h.
const STEPS = (1e-3, 2e-3, 4e-3, 8e-3)
values(a; η = ε -> 0.0) =
    Dict(ε => exp(a * ε) + η(ε) for ε in (0.0, (s * h for h in STEPS for s in (1, -1))...))
# A noise of `δ` at ±h with the sign of ε, which moves D(h) by δ/h.
kick(pairs...) = ε -> sum(sign(ε) * δ for (h, δ) in pairs if abs(ε) == h; init = 0.0)

@testset "the central difference, the spread and the extrapolant" begin
    J = values(3.0)
    @test central(J, 2e-3) ≈ sinh(3 * 2e-3) / 2e-3 rtol = 1e-12
    @test spread([1.0, 1.04, 1.02]) ≈ 0.04
    @test spread([-1.0, -1.04]) ≈ 0.04
    @test spread([1.0, -1.0]) == Inf
    @test spread([0.0, 1.0]) == Inf
    @test extrapolant([1.0, 2.0], [1e-3, 2e-3], 1) ≈ (4 * 1.0 - 2.0) / 3
    @test extrapolant([1.0, 2.0], [1e-3, 3e-3], 1) ≈ (9 * 1.0 - 2.0) / 8
end

@testset "a smooth objective converges and Richardson removes the truncation" begin
    for a in (50.0, -50.0)
        c = convergence(values(a), STEPS)
        # The truncation is 2.7 % at the largest step and 0.04 % at the smallest.
        @test c.D[end] / a - 1 ≈ (0.4^2 / 6) rtol = 0.02
        @test c.window == [1e-3, 2e-3, 4e-3]
        @test c.value ≈ a * (1 - 0.05^4 / 30) rtol = 1e-9
        @test c.error ≈ abs(c.extrapolants[1] - c.extrapolants[2])
        @test c.error <= 1e-5 * abs(a)
        @test usable_steps(c, c.value) == collect(STEPS)
        # The literal reading agrees here: the truncation is within 5 % at the
        # largest step.
        @test c.agreeing == collect(STEPS)
        @test c.literal_value ≈ c.value
    end
    # Truncation of 7 % at the largest step: the literal reading finds no
    # agreement among the three largest, Richardson's converges.
    a = 80.0
    c = convergence(values(a), STEPS)
    @test spread(c.D[2:4]) > AGREE
    @test isempty(c.agreeing) && isnan(c.literal_value)
    @test c.window == [1e-3, 2e-3, 4e-3]
    @test c.value ≈ a * (1 - 0.08^4 / 30) rtol = 1e-8
    @test usable_steps(c, c.value) == [1e-3, 2e-3, 4e-3]
end

@testset "noise: its steps are left out of the window and are unusable" begin
    a = 50.0
    # D(1e-3) moved by 10 % of the derivative.
    c = convergence(values(a; η = kick(1e-3 => 0.1a * 1e-3)), STEPS)
    @test c.D[1] ≈ sinh(a * 1e-3) / 1e-3 + 0.1a rtol = 1e-12
    @test c.window == [2e-3, 4e-3, 8e-3]
    @test c.value ≈ a * (1 - 0.1^4 / 30) rtol = 1e-8
    @test usable_steps(c, c.value) == [2e-3, 4e-3, 8e-3]
    # Noise at the three largest steps: no window agrees.
    erratic = convergence(
        values(a; η = kick(2e-3 => 0.2a * 2e-3, 4e-3 => -0.2a * 4e-3, 8e-3 => 0.2a * 8e-3)),
        STEPS,
    )
    @test isempty(erratic.window) && isnan(erratic.value) && isnan(erratic.error)
    @test isempty(usable_steps(erratic, erratic.value))
end

@testset "the tolerances and the verdict" begin
    a = 50.0
    tight = convergence(values(a), STEPS)
    w = tight.window
    # A looser tolerance 1 % off at every step agrees on the window; 10 % off
    # does not; noise outside the window does not count, inside it does.
    @test tolerances_agree(convergence(values(1.01a), STEPS), tight, w)
    @test !tolerances_agree(convergence(values(1.1a), STEPS), tight, w)
    noisy_large = convergence(values(a; η = kick(8e-3 => 0.3a * 8e-3)), STEPS)
    noisy_small = convergence(values(a; η = kick(1e-3 => 0.3a * 1e-3)), STEPS)
    @test tolerances_agree(noisy_large, tight, w)
    @test !tolerances_agree(noisy_small, tight, w)
    @test !tolerances_agree(tight, tight, Float64[])
    # The verdict, loosest first, and the usable steps per tolerance.
    loose =
        convergence(values(a; η = kick(1e-3 => 0.3a * 1e-3, 2e-3 => 0.3a * 2e-3)), STEPS)
    mid = convergence(values(1.01a), STEPS)
    q = qualification([loose, mid, tight])
    @test q.pass && q.converges && q.across
    @test q.value == tight.value && q.error == tight.error && q.window == w
    @test q.usable == [[4e-3, 8e-3], collect(STEPS), collect(STEPS)]
    @test q.literal.pass && q.literal.value ≈ q.value
    # The two tightest tolerances disagreeing fail the verdict.
    @test !qualification([loose, convergence(values(1.1a), STEPS), tight]).pass
    # Truncation past the literal reading's reach: Richardson's verdict
    # passes, the literal one fails.
    steep = [convergence(values(80.0), STEPS), convergence(values(80.0), STEPS)]
    q = qualification(steep)
    @test q.pass && !q.literal.pass && !q.literal.agrees
    # No window at the tightest tolerance, or one tolerance only, fails it.
    erratic = convergence(
        values(a; η = kick(2e-3 => 0.2a * 2e-3, 4e-3 => -0.2a * 4e-3, 8e-3 => 0.2a * 8e-3)),
        STEPS,
    )
    q = qualification([tight, erratic])
    @test !q.pass && !q.converges && isnan(q.value)
    @test all(isempty, q.usable)
    @test !qualification([tight]).pass
    @test_throws ArgumentError convergence(values(a), (1e-3, 2e-3))
end
