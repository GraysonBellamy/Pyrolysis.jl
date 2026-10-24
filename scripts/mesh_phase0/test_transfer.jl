@testset "merge conservation, porosity and pressure" begin
    cp = [1.0, 1.0, 1.0]
    densities = [T -> 300 / T, T -> 1.0]
    @testset "original pore-free counterexample" begin
        a = (; m = [1.0, 0.0, 0.0], T = 300.0, V = 1.0)
        b = (; m = [0.0, 1.0, 0.0], T = 500.0, V = 1.0)
        merged = merge_state(a, b, cp, densities)
        @test merged.T ≈ 400 atol = 1e-12
        @test 1 - merged.Vs / 2 ≈ -1 / 6 atol = 1e-12
        @test merged.V ≈ 7 / 3 atol = 1e-12
        @test merged.phi ≈ 0 atol = 1e-12
        @test merged.m == a.m + b.m
        @test merged.H ≈ 200 atol = 1e-12
    end
    @testset "gas-containing merge" begin
        for phi in (0.01, 0.05, 0.1), gas_mass in (1e-6, 1e-3)
            a = (; m = [1 - phi, 0.0, gas_mass], T = 300.0, V = 1.0)
            b = (; m = [0.0, 1 - phi, gas_mass], T = 500.0, V = 1.0)
            merged = merge_state(a, b, cp, densities)
            @test merged.T ≈ 400 atol = 1e-12
            # Pore volume is the sum of the inputs' pore volumes; created volume
            # is the skeletal thermal expansion at the merged temperature only.
            @test merged.V - merged.Vs ≈ 2phi atol = 1e-12
            @test merged.V ≈ 7 / 3 - phi / 3 atol = 1e-12
            @test merged.phi ≈ 2phi / (7 / 3 - phi / 3) atol = 1e-12
            @test isfinite(merged.pressure) && merged.pressure > 0
            @test merged.m == a.m + b.m
            expected_H = dot(b.m, cp) * 200
            @test dot(merged.m, cp) * (merged.T - 300) ≈ expected_H atol = 1e-12
            expected_pressure = (2gas_mass / 0.028) * 8.314462618 * 400 / (2phi)
            @test merged.pressure ≈ expected_pressure rtol = 1e-12
            @test merged.created_volume ≈ 1 / 3 - phi / 3 atol = 1e-12
        end
    end
    @testset "one pore-free gas-free donor plus one porous gas-containing cell" begin
        # Both inputs are admissible. A rule keyed to the smaller input porosity
        # would reserve no pore volume here; the pore-preserving rule keeps the
        # donor's 0.05 and the retained gas stays at a finite pressure.
        a = (; m = [1.0, 0.0, 0.0], T = 300.0, V = 1.0)
        b = (; m = [0.0, 0.95, 1e-3], T = 500.0, V = 1.0)
        @test 1 - a.m[1] / densities[1](a.T) == 0
        @test 1 - b.m[2] / densities[2](b.T) > 0
        merged = merge_state(a, b, cp, densities)
        @test merged.Vs > a.V + b.V
        @test merged.V - merged.Vs ≈ 0.05 atol = 1e-12
        @test merged.phi > 0
        metric("merge: pressure with a pore-free donor", merged.pressure)
        @test isfinite(merged.pressure) && merged.pressure > 0
        # Unequal masses, so the merged temperature is below 400 K.
        @test merged.T ≈ 300 + 0.951 * 200 / 1.951 atol = 1e-12
        @test merged.pressure ≈ (1e-3 / 0.028) * 8.314462618 * merged.T / 0.05 rtol = 1e-12
        @test merged.created_volume ≈ merged.Vs - (1.0 + 0.95) atol = 1e-12
    end
end
