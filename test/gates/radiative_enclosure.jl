# A gray slab in a radiative enclosure at its own temperature stays there:
# under RADIATIVE_TRANSFER the rows emit and absorb with the same emissivity
# and absorption coefficient, so a slab at T_h that faces surroundings at T_h
# exchanges nothing, at any optical thickness. Three enclosures:
#
#   E1  diffuse σT_h⁴ onto the exposed face, an open back facing T_h;
#   E2  the same exposed face, a mirror back (adiabatic);
#   E3  as E2 behind an opaque 10 µm back layer of another material.
#
# The exposed face's incoming radiation is a diffuse `RadiativeFluxBC(σT_h⁴)`
# beside a `RadiativeBC` facing 1 K, so the flux and the surroundings both
# enter the rows; the 1 K surroundings add σ·1⁴, which settles the slab
# 1/(4T_h³) = 7.3×10⁻¹⁰ K above T_h behind a mirror and half that with E1's
# open back. Each slab of κL ∈ {0.03, …, 10} is integrated to 30 of its
# radiative time constants ρcL/(4εσT_h³) in both energy forms; every cell
# at every saved time must stay within TOLERANCE of T_h.
#
#   julia --project=test test/gates/radiative_enclosure.jl

using Test
using Printf
using Pyrolysis
using Pyrolysis.Materials: RADIATIVE_TRANSFER, STEFAN_BOLTZMANN

const TOLERANCE = 1e-8      # K
const T_H = 700.0           # K
const ρ = 1210.0            # kg/m³
const κ = 2870.0            # 1/m
const ε = 0.96
const BACK_LAYER = 1e-5     # m
const N_CELLS = 40

function enclosure(κL, case)
    L = κL / κ
    slab = SolidComponent(:slab, ρ = ρ, c = 1500.0, k = 0.2, ε = ε, α = κ / ρ)
    components = (slab,)
    nodes = collect(range(0.0, L; length = N_CELLS + 1))
    fields = Any[z->ρ]
    if case === :E3
        back = SolidComponent(:back, ρ = 1000.0, c = 1000.0, k = 1.0, ε = ε, α = 1e6)
        components = (slab, back)
        nodes = [0.0; BACK_LAYER .+ nodes]
        inside(z) = z > BACK_LAYER - 1e-12
        fields = Any[z->inside(z) ? ρ : 0.0, z->inside(z) ? 0.0 : 1000.0]
    end
    material = Material(name = :slab, components = components, reactions = ())
    top = RadiativeFluxBC(STEFAN_BOLTZMANN * T_H^4) + RadiativeBC(T_inf = 1.0)
    bottom = case === :E1 ? RadiativeBC(T_inf = T_H) : AdiabaticBC()
    t_end = 30 * ρ * 1500.0 * max(L, 1e-4) / (4ε * STEFAN_BOLTZMANN * T_H^3)
    return PyrolysisProblem(
        mesh = line_mesh(nodes),
        material = material,
        bc_set = BoundaryConditionSet(thermal = (top = top, bottom = bottom)),
        T_initial = z -> T_H,
        ξ_initial = fields,
        tspan = (0.0, t_end),
        output_times = [0.0; t_end .* exp10.(-6:0)],
    )
end

println("case  κL     form          max|T − T_h| [K]  retcode  accepted")
@testset "radiative enclosure: $(case), κL = $(κL), $(energy) form" for case in
                                                                        (:E1, :E2, :E3),
    κL in (0.03, 0.1, 0.3, 1.0, 3.0, 10.0),
    energy in (:temperature, :enthalpy)

    sol = solve(
        enclosure(κL, case);
        radiation_model = RADIATIVE_TRANSFER,
        energy_state = energy,
    )
    deviation = maximum(abs, sol.T .- T_H)
    @printf(
        "%-4s %5.2f  %-12s  %14.3e    %-8s %d\n",
        case,
        κL,
        energy,
        deviation,
        sol.retcode,
        sol.statistics.accepted
    )
    flush(stdout)
    @test sol.retcode === :Success
    @test deviation <= TOLERANCE
end
