# =============================================================================
# V-RAD-07 reference: in-depth heating with emission by an independent method
# of lines
# =============================================================================
#
# The slab of V-RAD-07 solved without the package: a gray solid (δ = 1 cm,
# κ = 1000 m⁻¹, ρ = 1000 kg/m³, c = 1400 J/kg·K, k = 0.2 W/m·K, ε = 1) at
# T₀ = 300 K, a collimated heater q₀ = 25 kW/m² and surroundings at T₀ through
# its exposed face, convection h = 20 W/m²·K to T₀ there, an adiabatic back
# that reflects all radiation (a mirror).
#
# The discretization is the solver's, written out independently: cell-centred
# finite volumes on a uniform mesh of N cells; conduction k ΔT/Δ across
# interior faces; at the exposed face the surface temperature of the balance
# h(T∞ − T_s) = k(T_s − T₁)/(Δ/2), which is linear and so taken in closed
# form; the radiation rows of exact/two_stream.jl (gray two-stream,
# D = 3/2, exact across each cell for its emissive power σTᵢ⁴, the heater a
# beam e^{−τ}, what reaches the mirror joining q⁺). The rows are eliminated at
# every evaluation by the oracle's direct sweep, so the reference is an
# ordinary differential system in the cell temperatures, where the solver
# carries the fluxes as algebraic unknowns. Integrated by Kvaerno5 (an
# ESDIRK method, not the solver's Rosenbrock default) at reltol 1e-12 with a
# finite-difference Jacobian; the integrator's own error is checked against
# FBDF at the same tolerance and printed.
#
# Writes vrad07_n<N>.csv beside this file for N = 50, 100, 200 (the case's
# ladder, every 5 s to 180 s) and N = 1600 (the continuum proxy of its order
# check, at 30, 60, 90 and 180 s): a names row and a units row, then one row
# per saved time: the time, the surface temperature, and the cell
# temperatures from the back (cell 1, z = 0) to the exposed face (cell N).
#
# Run from the repository root:
#   julia --project=test --startup-file=no test/verification/reference/radiation/vrad07_reference.jl

using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK: Kvaerno5
using Printf

const EXACT = joinpath(@__DIR__, "..", "..", "exact")
include(joinpath(EXACT, "gauss_legendre.jl"))
include(joinpath(EXACT, "two_stream.jl"))

const σ = 5.670374419e-8   # W/m²·K⁴ (CODATA 2018)
const δ, κ, ρ, c, k = 0.01, 1000.0, 1000.0, 1400.0, 0.2
const h, q0, T0 = 20.0, 25.0e3, 300.0
const T_END = 180.0

# The cell temperatures' rates, cells numbered from the exposed face (the
# oracle's order). Returns the surface temperature when asked.
function rates!(dT, T, N)
    Δ = δ / N
    T_s = (h * T0 + (k / (Δ / 2)) * T[1]) / (h + k / (Δ / 2))
    r = two_stream_solve(
        fill(κ * Δ, N),
        σ .* T .^ 4;
        D = 1.5,
        ε = 1.0,
        Q = σ * T0^4,
        ε_b = 0.0,
        beams = [TwoStreamBeam(1.0, q0)],
    )
    for i = 1:N
        top = i == 1 ? h * (T0 - T_s) : k * (T[i-1] - T[i]) / Δ
        bottom = i == N ? 0.0 : k * (T[i] - T[i+1]) / Δ
        dT[i] = (top - bottom + r.S[i]) / (ρ * c * Δ)
    end
    return T_s
end

surface(T, N) = (h * T0 + (2k / (δ / N)) * T[1]) / (h + 2k / (δ / N))

function reference(N, saves, method)
    problem =
        ODEProblem((dT, T, p, t) -> (rates!(dT, T, N); nothing), fill(T0, N), (0.0, T_END))
    sol = solve(problem, method; reltol = 1e-12, abstol = 1e-10, saveat = saves)
    sol.retcode == ReturnCode.Success || error("N = $N: $(sol.retcode)")
    return sol
end

function write_reference(N, saves)
    sol = reference(N, saves, Kvaerno5(autodiff = AutoFiniteDiff()))
    check = reference(N, saves, FBDF(autodiff = AutoFiniteDiff()))
    spread = maximum(maximum(abs, a .- b) for (a, b) in zip(sol.u, check.u))
    path = joinpath(@__DIR__, "vrad07_n$(N).csv")
    open(path, "w") do io
        println(io, join(["Time"; "T_s"; ["T_$i" for i = 1:N]], ","))
        println(io, join(["s"; fill("K", N + 1)], ","))
        for (t, T) in zip(sol.t, sol.u)
            # Written from the back (z = 0) to the exposed face, the solver's order.
            values = [t; surface(T, N); reverse(T)]
            println(io, join([@sprintf("%.15g", v) for v in values], ","))
        end
    end
    i_peak = argmax(sol.u[end])
    @printf(
        "N = %4d: %s, %d saves; Kvaerno5 − FBDF ≤ %.2e K; at %g s T_s = %.6f K, peak %.6f K in cell %d from the exposed face\n",
        N,
        basename(path),
        length(sol.t),
        spread,
        sol.t[end],
        surface(sol.u[end], N),
        sol.u[end][i_peak],
        i_peak,
    )
    return nothing
end

for N in (50, 100, 200)
    write_reference(N, collect(0.0:5.0:T_END))
end
write_reference(1600, [30.0, 60.0, 90.0, T_END])
