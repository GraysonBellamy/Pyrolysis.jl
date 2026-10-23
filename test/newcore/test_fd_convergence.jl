using Test

include(joinpath(@__DIR__, "gates", "fd_convergence.jl"))

# W10's convergence rule (plan §5.3 item 16) on the synthetic objective
# J(ε) = exp(a ε), whose derivative at ε = 0 is `a` and whose central
# difference is sinh(a h)/h = a (1 + (a h)²/6 + …), plus a noise `η(ε)` of
# known values: the values the gate's runs give at ε = 0 and ±h.
const STEPS = (1e-3, 2e-3, 4e-3, 8e-3)
values(a; η = ε -> 0.0) =
    Dict(ε => exp(a * ε) + η(ε) for ε in (0.0, (s * h for h in STEPS for s in (1, -1))...))
# A noise of `δ` at ±h with the sign of ε, which moves D(h) by δ/h.
kick(pairs...) = ε -> sum(sign(ε) * δ for (h, δ) in pairs if abs(ε) == h; init = 0.0)

@testset "the central difference and the spread" begin
    J = values(3.0)
    @test central(J, 2e-3) ≈ sinh(3 * 2e-3) / 2e-3 rtol = 1e-12
    @test spread([1.0, 1.04, 1.02]) ≈ 0.04
    @test spread([-1.0, -1.04]) ≈ 0.04
    @test spread([1.0, -1.0]) == Inf
    @test spread([0.0, 1.0]) == Inf
end

@testset "a smooth objective converges at every step and Richardson removes the truncation" begin
    for a in (50.0, -50.0)
        c = convergence(values(a), STEPS)
        # The truncation is 2.7 % at the largest step and 0.04 % at the smallest.
        @test c.D[end] / a - 1 ≈ (0.4^2 / 6) rtol = 0.02
        @test c.agreeing == collect(STEPS)
        @test c.value ≈ a rtol = 1e-6
        @test c.error ≈ abs(c.D[1] - c.D[2]) / 3
        @test smallest_usable(c, c.value) == 1e-3
    end
    # Truncation beyond the agreement at the largest step: the three largest
    # disagree, so the rule reports no converged value.
    c = convergence(values(80.0), STEPS)
    @test spread(c.D[2:4]) > AGREE
    @test isempty(c.agreeing) && isnan(c.value)
    @test isnan(smallest_usable(c, c.value))
end

@testset "noise at the smallest step: its difference is unusable, the rest converge" begin
    a = 50.0
    # D(1e-3) moved by 10 % of the derivative.
    c = convergence(values(a; η = kick(1e-3 => 0.1a * 1e-3)), STEPS)
    @test c.D[1] ≈ sinh(a * 1e-3) / 1e-3 + 0.1a rtol = 1e-12
    @test c.agreeing == [2e-3, 4e-3, 8e-3]
    # Richardson from 2e-3 and 4e-3 leaves (a h)⁴/30 at h = 2e-3.
    @test c.value ≈ a * (1 - 0.1^4 / 30) rtol = 1e-8
    @test smallest_usable(c, c.value) == 2e-3
    # Noise that breaks the three largest: no converged value.
    erratic = convergence(
        values(a; η = kick(2e-3 => 0.2a * 2e-3, 4e-3 => -0.2a * 4e-3, 8e-3 => 0.2a * 8e-3)),
        STEPS,
    )
    @test isempty(erratic.agreeing) && isnan(erratic.value)
    # A step that falls between agreeing ones is not taken past a disagreeing
    # smaller one: the set extends down only while it keeps agreeing.
    gap = convergence(values(a; η = kick(1e-3 => -0.2a * 1e-3)), STEPS)
    @test gap.agreeing == [2e-3, 4e-3, 8e-3]
end

@testset "the tolerances and the verdict" begin
    a = 50.0
    tight = convergence(values(a), STEPS)
    # A looser tolerance 1 % off at every step agrees; 10 % off does not; noise
    # at its smallest step alone does not count against it.
    @test tolerances_agree(convergence(values(1.01a), STEPS), tight)
    @test !tolerances_agree(convergence(values(1.1a), STEPS), tight)
    @test tolerances_agree(
        convergence(values(a; η = kick(1e-3 => 0.3a * 1e-3)), STEPS),
        tight,
    )
    @test !tolerances_agree(
        convergence(values(a; η = kick(8e-3 => 0.1a * 8e-3)), STEPS),
        tight,
    )
    # The verdict, loosest first: the loosest tolerance's noise leaves its
    # smallest steps unusable against the tightest's value.
    loose =
        convergence(values(a; η = kick(1e-3 => 0.3a * 1e-3, 2e-3 => 0.3a * 2e-3)), STEPS)
    mid = convergence(values(a; η = kick(1e-3 => 0.1a * 1e-3)), STEPS)
    q = qualification([loose, mid, tight])
    @test q.pass && q.agrees && q.across
    @test q.value == tight.value
    @test q.usable == [4e-3, 2e-3, 1e-3]
    # The two tightest tolerances disagreeing fail the verdict.
    @test !qualification([loose, convergence(values(1.1a), STEPS), tight]).pass
    # No converged value at the tightest tolerance fails it, whatever the rest.
    bad = convergence(values(80.0), STEPS)
    q = qualification([tight, bad])
    @test !q.pass && !q.agrees && isnan(q.value)
    @test all(isnan, q.usable)
    @test_throws ArgumentError convergence(values(a), (1e-3, 2e-3))
end
