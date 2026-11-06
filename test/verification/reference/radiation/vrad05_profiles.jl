# =============================================================================
# V-RAD-05 frames: burning PMMA profiles under radiative transfer
# =============================================================================
#
# Runs the MaCFP-2 PMMA gasification target with the GIDAZE+ property set
# (the deck of cases/V-RAD-05.jl): 6 mm on 120 cells, a collimated heater of
# 10 and 25 kW/m², surroundings at 293 K, no convection, a mirror back, under
# RADIATIVE_TRANSFER on a moving mesh with deactivation, saved every second.
# At the first save at which 50 % and then 10 % of the PMMA mass is left it
# writes the live cells, from the back to the exposed face, to
# vrad05_q<q>_left<percent>.csv beside this file: the lower and upper node
# positions [m], the temperature [K] and the PMMA concentration [kg/m³].
#
# Run from the repository root:
#   julia --project=test --startup-file=no test/verification/reference/radiation/vrad05_profiles.jl

using Printf

const VERIFICATION = joinpath(@__DIR__, "..", "..")
include(joinpath(VERIFICATION, "casedefs.jl"))
include(joinpath(VERIFICATION, "cases", "V-RAD-05.jl"))

const L, N = 6e-3, 120

for (q_kW, t_end) in ((10, 8000.0), (25, 2000.0))
    problem = _vrad_gidaze_problem(
        collect(range(0.0, L; length = N + 1)),
        z -> _VRAD_GIDAZE_T0,
        z -> _VRAD_GIDAZE_ρ;
        q = 1e3 * q_kW,
        tspan = (0.0, t_end),
        output_times = collect(0.0:1.0:t_end),
    )
    sol = solve(
        problem;
        radiation_model = RADIATIVE_TRANSFER,
        use_ale = true,
        deactivation = true,
    )
    mass = [sum(sol.m[1, :, j]) for j in eachindex(sol.t)]
    for left in (50, 10)
        j = findfirst(≤(left / 100 * mass[1]), mass)
        j === nothing && error("$(q_kW) kW/m²: $(left) % of the mass is never reached")
        live = findall(!isnan, sol.T[:, j])
        @assert live == 1:length(live)
        path = joinpath(@__DIR__, "vrad05_q$(q_kW)_left$(left).csv")
        open(path, "w") do io
            println(io, "z_lower,z_upper,T,xi_PMMA")
            println(io, "m,m,K,kg/m3")
            for i in live
                values = (sol.z[i, j], sol.z[i+1, j], sol.T[i, j], sol.ξ[1, i, j])
                println(io, join([@sprintf("%.17g", v) for v in values], ","))
            end
        end
        @printf(
            "%2d kW/m², %2d %% left: t = %g s, mass left %.4f, %d live cells, thickness %.4f mm, T_s %.2f K -> %s\n",
            q_kW,
            left,
            sol.t[j],
            mass[j] / mass[1],
            length(live),
            1e3 * (sol.z[length(live)+1, j] - sol.z[1, j]),
            sol.surface_T[j],
            basename(path),
        )
    end
end
