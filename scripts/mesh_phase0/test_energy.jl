@testset "temperature versus enthalpy: closed stationary cells" begin
    # Nondimensional temperatures: cp1=T, cp2=1, unit masses. No reactions,
    # gas flux or mesh motion. Analytically d(T1^2/2 + T2)/dt = 0.
    function rhs_temperature!(du, u, p, t)
        q = u[2] - u[1]
        du[1] = q / u[1]
        du[2] = -q
    end
    function rhs_inventory!(du, u, p, t)
        q = u[2] - sqrt(2u[1])
        du[1] = q
        du[2] = -q
    end
    errors = Float64[]
    for dt in (0.2, 0.1, 0.05)
        st = solve(
            ODEProblem(rhs_temperature!, [1.0, 3.0], (0.0, 1.0)),
            KenCarp4();
            adaptive = false,
            dt,
            reltol = 1e-12,
            abstol = 1e-14,
        )
        sh = solve(
            ODEProblem(rhs_inventory!, [0.5, 3.0], (0.0, 1.0)),
            KenCarp4();
            adaptive = false,
            dt,
            reltol = 1e-12,
            abstol = 1e-14,
        )
        @test SciMLBase.successful_retcode(st)
        @test SciMLBase.successful_retcode(sh)
        error_T = abs(st.u[end][1]^2 / 2 + st.u[end][2] - 3.5)
        error_H = maximum(abs(sum(u) - 3.5) for u in sh.u)
        push!(errors, error_T)
        # Test the invariant/order, not the historical probe's exact numbers.
        @test error_H <= 100eps(3.5)
        @test error_T > 1e-8
        metric("energy: dt=$dt, T error", error_T)
        metric("energy: dt=$dt, H error", error_H)
    end
    orders = log2.(errors[1:2] ./ errors[2:3])
    @test all(p -> 3.5 < p < 4.5, orders)
    metric("energy: observed orders", orders)
end
