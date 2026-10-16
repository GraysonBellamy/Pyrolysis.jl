"""
Regression tests for boundary heat-flux reporting.

Pre-fix, `surface_heat_flux_net` contradicted the applied flux
regime-dependently: a standalone `RadiativeFluxBC` reported net = 0 under
SURFACE_ABSORPTION although the Newton balance applied α·q_ext, while
`CombinedThermalBC` discarded children's nets and recomputed
net = absorbed − rerad − conv + prescribed — re-including α·q_ext that the
BEER_LAMBERT dispatch routes through the volume (~24–45 kW/m² over-report
in cone configurations).

Post-fix convention: per-BC `net` is the flux the surface energy balance
actually applies at the face (it matches `evaluate_bc_flux_and_derivative`
for every model), `CombinedThermalBC` sums children's nets, and the
Beer–Lambert volumetric deposit is reported separately as
`extras.volumetric_absorption`. The ledger identity `net = q_conduction`
(at the solved T_s) therefore holds in every radiation regime.
"""

using Test
using Pyrolysis
using Pyrolysis.Materials: NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT
using Pyrolysis.Discretization: compute_boundary_heat_fluxes
import Pyrolysis.Discretization: evaluate_bc_flux_and_derivative

@testset "surface gas reporting includes absent-species ambient inflow" begin
    mat = Material(
        name = :gas_reporting,
        components = (
            SolidComponent(
                :solid;
                ρ = 500.0,
                ρ_intrinsic = 1500.0,
                c = 1000.0,
                k = 0.2,
                κ = 1e-12,
            ),
            GaseousComponent(:AIR; M = 0.02897, c = 1006.0, k = 0.026, λ = 1e-5),
            GaseousComponent(:other; M = 0.028, c = 1000.0, k = 0.03, λ = 1e-5),
        ),
        reactions = (),
    )
    mass_bc =
        ConvectiveMassBC(h_m = 0.05, ξ_inf = (0.0, 0.7, 0.0)) +
        DirichletPressureBC(101325.0)
    for ale in (false, true)
        mesh = create_uniform_mesh(
            0.01,
            2,
            Val(3);
            T_initial = 300.0,
            ξ_initial = (500.0, 0.7, 0.0),
        )
        problem = PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = BoundaryConditionSet(
                thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
                mass = (top = mass_bc, bottom = mass_bc),
            ),
            tspan = (0.0, 1.0),
        )
        def, ws, u = initialize_state(problem; use_ale = ale)
        for air in (-1e-12, 0.0, 0.2), other in (0.0, 2.0)
            u[1:2] .= 300.0
            u[3:4] .= 500.0
            u[5:6] .= air
            u[7:8] .= other
            residual!(similar(u), u, def, 0.0, ws)
            ξ = (500.0, air, other)
            for (side, face, sign) in ((:top, 3, 1.0), (:bottom, 1, -1.0))
                reported = Pyrolysis.Solver._surface_gas_flux(
                    mat,
                    ξ,
                    300.0,
                    mass_bc,
                    0.0025,
                    0.0,
                    true,
                    Val(side),
                )
                @test collect(reported) ≈ sign .* ws.rhs_cache.J_gas[:, face] atol=1e-14 rtol=1e-12
                air == 0.0 && @test reported[2] ≈ -0.035
            end
        end
    end
end

function _reporting_material()
    virgin = SolidComponent(
        :virgin,
        ρ = 600.0,
        c = 1500.0,
        k = 0.15,
        ε = 0.95,
        α = 200.0,
        λ = 1e-5,
    )
    gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5)
    rxn = Reaction(:pyro, 1 => 2, A = 1e8, E = 150e3, h = -300e3)
    return Material(name = :wood, components = (virgin, gas), reactions = (rxn,))
end

@testset "net is the applied face flux, per BC and per model" begin
    combined =
        RadiativeFluxBC(flux = 50e3) +
        RadiativeBC(T_inf = 300.0) +
        ConvectiveBC(h = 10.0, T_inf = 300.0)
    standalone = RadiativeFluxBC(flux = 50e3)
    ε_eff = 0.9
    T_s = 700.0
    t = 0.0

    for bc in (combined, standalone),
        model in (NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT)

        hf = compute_boundary_heat_fluxes(bc, T_s, t, ε_eff, model)
        q_applied, _ = evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff, model)
        @test hf.net ≈ q_applied atol = 1e-9
    end

    # Worst pre-fix cases, pinned explicitly:
    # (a) standalone RadiativeFluxBC under SURFACE_ABSORPTION reported 0
    #     although the Newton applied ε_eff·q_ext.
    hf = compute_boundary_heat_fluxes(standalone, T_s, t, ε_eff, SURFACE_ABSORPTION)
    @test hf.net ≈ ε_eff * 50e3
    @test hf.absorbed ≈ ε_eff * 50e3

    # (b) combined under BEER_LAMBERT over-reported by α·q_ext: net must
    #     carry only the rerad + convection the face actually applies, and
    #     the face absorbs nothing (the deposit is volumetric).
    hf_bl = compute_boundary_heat_fluxes(combined, T_s, t, ε_eff, BEER_LAMBERT)
    σ = 5.670374419e-8
    q_rerad = ε_eff * σ * (300.0^4 - T_s^4)
    q_conv = 10.0 * (300.0 - T_s)
    @test hf_bl.absorbed == 0.0
    @test hf_bl.net ≈ q_rerad + q_conv rtol = 1e-12
    @test abs(
        hf_bl.net -
        (hf_bl.absorbed - hf_bl.reradiation - hf_bl.convection + hf_bl.prescribed),
    ) < 1e-9  # components stay consistent with net

    # (c) combined under SURFACE_ABSORPTION includes the absorbed flux.
    hf_sa = compute_boundary_heat_fluxes(combined, T_s, t, ε_eff, SURFACE_ABSORPTION)
    @test hf_sa.net ≈ ε_eff * 50e3 + q_rerad + q_conv rtol = 1e-12
end

@testset "solved extras close the surface ledger ($(label))" for (label, model) in (
    ("BEER_LAMBERT", BEER_LAMBERT),
    ("SURFACE_ABSORPTION", SURFACE_ABSORPTION),
)
    mat = _reporting_material()
    mesh = create_uniform_mesh(0.005, 8, Val(2))
    bc = BoundaryConditionSet(
        thermal = (
            top = RadiativeFluxBC(flux = 50e3) +
                  RadiativeBC(T_inf = 300.0) +
                  ConvectiveBC(h = 10.0, T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
        mass = (top = DiffusiveMassBC(), bottom = ImpermeableBC()),
    )
    problem = PyrolysisProblem(
        mesh = mesh,
        material = mat,
        bc_set = bc,
        tspan = (0.0, 10.0),
        output_times = collect(0.0:2.0:10.0),
    )
    sol = solve(problem; radiation_model = model)
    @test sol.retcode == :Success

    ex = sol.extras
    for i in eachindex(sol.t)
        # Ledger identity: applied net = conduction at the solved T_s.
        # Re-including α·q_ext under BEER_LAMBERT would break this by ~45 kW/m².
        @test ex.surface_heat_flux_net[i] ≈ ex.surface_heat_flux_conduction[i] atol =
            max(1.0, 1e-6 * abs(ex.surface_heat_flux_net[i]))
    end

    if model == BEER_LAMBERT
        # The volumetric deposit is reported separately, bounded by ε·q_ext,
        # and strictly positive for this optically thick sample.
        @test all(ex.volumetric_absorption .> 0.0)
        @test all(ex.volumetric_absorption .<= ex.surface_emissivity .* 50e3 .+ 1e-9)
        # The face never absorbs incident radiation under BEER_LAMBERT.
        @test all(iszero, ex.surface_heat_flux_absorbed)
    else
        # SURFACE_ABSORPTION: absorption happens at the face and is inside
        # net; there is no volumetric deposit.
        @test all(iszero, ex.volumetric_absorption)
        @test all(ex.surface_heat_flux_absorbed .> 0.0)
    end
end
