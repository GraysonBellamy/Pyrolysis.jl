@testset "extended caloric map (constant-cp instance)" begin
    delta, floor = 1e-8, 1e-12
    cp = [1000.0, 1500.0]
    for m in (-1.0, 0.0, delta / 10, delta, 10delta, 1.0)
        @test positive_mass(m, delta) >= 0
        @test isfinite(ForwardDiff.derivative(x -> positive_mass(x, delta), m))
        if m <= 0
            @test positive_mass(m, delta) == 0
        elseif m <= delta
            @test abs(positive_mass(m, delta) - m) <= delta
        else
            @test abs(positive_mass(m, delta) - m) <= delta^2 / m + eps(m)
        end
    end
    @test ForwardDiff.derivative(m -> positive_mass(m, delta), 0.0) == 0
    # Analytic derivatives are independent of AD through the candidate map.
    m = [1e-4, 0.2]
    H = 10.0
    # Both masses sit above delta, where the positive part is exactly m.
    C = sum(cp[j] * m[j] for j = 1:2)
    expected_derivative = [1 / C; -H .* cp ./ C^2]
    actual_derivative = ForwardDiff.gradient(
        v -> temperature_from_inventory(v[1], v[2:3], cp; delta, floor),
        [H; m],
    )
    @test actual_derivative ≈ expected_derivative rtol = 1e-12
    for m in ([0.0, 0.0], [-1.0, -2.0], [1.0, -1.0], [1.0, 0.0], [delta, delta]),
        H in (-100.0, 0.0, 100.0)

        f(v) = temperature_from_inventory(v[1], v[2:3], cp; delta, floor)
        v = [H; m]
        @test isfinite(f(v))
        @test all(isfinite, ForwardDiff.gradient(f, v))
        @test ForwardDiff.gradient(f, v)[1] > 0
    end
    # An absent species contributes exactly no effective heat capacity.
    @test temperature_from_inventory(1.0, [1.0, 0.0], cp; delta, floor) ==
          temperature_from_inventory(1.0, [1.0], cp[1:1]; delta, floor)
    # Independent physical-energy comparison bounds the stated regularization.
    for m in (delta / 2, delta, 10delta)
        dH = (m - positive_mass(m, delta)) * cp[1] * 100
        @test 0 <= dH <= delta * cp[1] * 100
    end
    # No inversion exists for the rejected unextended zero-mass relation H=0*T.
    @test all(T -> 0.0 * T != 1.0, (-1e6, 0.0, 1e6))

    # Region of exact agreement: every species absent or at least delta. A
    # rational positive part that only approached m failed this by 0.01 K.
    for present in (delta, 100delta, 1.0)
        masses = [present, 0.0]
        target_T = 400.0
        physical_H = dot(masses, cp) * (target_T - 300)
        recovered_T = temperature_from_inventory(physical_H, masses, cp; delta, floor)
        present == 100delta &&
            metric("caloric: T error with one absent species", recovered_T - target_T)
        @test isapprox(recovered_T, target_T; atol = 1e-12, rtol = 0)
    end
end
