using Test, Pyrolysis
using Pyrolysis.Materials: R_GAS

function exchange_problem()
    mat=Material(
        name = :exchange_ledger,
        components = (
            SolidComponent(
                :matrix;
                ρ = 500.0,
                ρ_intrinsic = 1000.0,
                c = 1000.0,
                k = 0.2,
                λ = 1e-5,
                κ = 1e-13,
            ),
            GaseousComponent(:product; M = 0.044, c = 1200.0, k = 0.03),
            GaseousComponent(:air; M = 0.029, c = 1000.0, k = 0.03),
        ),
    )
    PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, 4),
        material = mat,
        T_initial = z->500.0,
        ξ_initial = [
            z->500.0,
            z->0.5*0.7*1e5*0.044/(R_GAS*500),
            z->0.5*0.3*1e5*0.029/(R_GAS*500),
        ],
        bc_set = BoundaryConditionSet(
            thermal = (
                top = ConvectiveBC(h = 20.0, T_inf = 300.0, blowing = CouetteFilm()),
                bottom = AdiabaticBC(),
            ),
            mass = (
                top = AmbientGasBC(
                    GasReservoir(P = 1e5, T = 300.0, mole_fractions = (air = 1.0,)),
                ),
                bottom = ImpermeableBC(),
            ),
        ),
        tspan = (0.0, 0.1),
    )
end
@testset "saved pore state and exact species ledger ($energy)" for energy in
                                                                   (:temperature, :enthalpy)
    rec=ConservationRecorder()
    sol=solve(
        exchange_problem();
        energy_state = energy,
        callback = rec,
        saveat = 0.01,
        reltol = 1e-8,
        pressure_tolerance = 0.01,
    )
    @test completed(sol)
    @test sol.gas.species==[:product, :air]
    @test get_porosity(sol)≈fill(0.5, size(sol.T))
    p=get_partial_pressures(sol)
    @test dropdims(sum(p; dims = 1); dims = 1)≈get_pressure(sol)
    @test dropdims(sum(get_gas_composition(sol); dims = 1); dims = 1)≈ones(size(sol.T))
    for g = 1:2
        @test p[g, :, :]≈get_intrinsic_gas_densities(sol)[g, :, :] .* R_GAS .* sol.T ./
                         (g==1 ? 0.044 : 0.029)
    end
    top=get_boundary_exchange(sol, :top)
    @test top.energy≈top.enthalpy .- top.heat atol=1e-10
    @test vec(sum(top.fluxes; dims = 1))≈sol.boundary_mass_flux[:top]
    @test maximum(top.residual)<3e-12
    @test maximum(top.iterations)<60
    audit=gas_species_audit(sol, rec)
    @test audit.method===:recorder
    @test audit.species.product.escape>0
    @test audit.species.air.ingress>0
    for species in values(audit.species)
        @test abs(species.defect)<1e-7*max(species.initial, species.ingress, species.escape)
        @test species.production==0
        @test species.vented==0
    end
    @test abs(conservation_audit(sol, rec).energy.relative)<1e-6
end

function receding_gas_ledger()
    ρ=1000.0
    T0=300.0
    Tp=700.0
    E=1e5
    v=0.002
    A=v^2*E*(Tp-T0)*exp(E/(R_GAS*Tp))/(1e-6*R_GAS*Tp^2)
    mat=Material(
        name = :receding_species,
        components = (
            SolidComponent(:fuel; ρ = ρ, φ = 0.3, c = 1000.0, k = 1.0, λ = 1e-5),
            GaseousComponent(:product; M = 0.05, c = 1000.0, k = 0.02),
            GaseousComponent(:air; M = 0.029, c = 1000.0, k = 0.02),
        ),
        reactions = (Reaction(:ablate, 1=>2; A = A, E = E, h = 1e6),),
    )
    T=z->T0+1e5*z
    p=PyrolysisProblem(
        mesh = create_uniform_mesh(0.003, 6),
        material = mat,
        T_initial = T,
        ξ_initial = [
            z->ρ,
            z->0.0,
            z->background_gas_concentration(mat, :air, (ρ, 0.0, 0.0), T(z)),
        ],
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(2.8e6), bottom = AdiabaticBC()),
            mass = (
                top = AmbientGasBC(
                    GasReservoir(P = 1e5, T = 300.0, mole_fractions = (air = 1.0,)),
                ),
                bottom = ImpermeableBC(),
            ),
        ),
        tspan = (0.0, 0.8),
        output_times = collect(0.0:0.02:0.8),
    )
    return p
end
@testset "species and energy through moving and receding faces ($energy, $moving)" for energy in
                                                                                       (
        :temperature,
        :enthalpy,
    ),
    moving in (false, true)

    rec=ConservationRecorder()
    sol=solve(
        receding_gas_ledger();
        energy_state = energy,
        use_ale = moving,
        deactivation = true,
        reltol = 1e-8,
        abstol = 1e-10,
        pressure_tolerance = 0.01,
        dtmin = 1e-14,
        callback = rec,
    )
    @test completed(sol)
    @test !isempty(sol.events)
    @test sol.statistics.pressure_error_evaluations>0
    audit=gas_species_audit(sol, rec)
    @test audit.species.product.production>0
    @test audit.species.air.production==0
    for (name, species) in pairs(audit.species)
        @test abs(species.defect)<2e-5*max(
            species.initial,
            abs(species.production),
            species.ingress,
            species.escape,
            species.vented,
        )
        j=findfirst(c->c.name===name, sol.problem.material.components)
        @test species.vented≈sum(e.lost_mass[j] for e in sol.events) atol=1e-15
    end
    moving || @test audit.species.air.vented>0
    @test abs(conservation_audit(sol, rec).energy.relative)<2e-5
    @test all(x->!isfinite(x) || x>2e-12, get_porosity(sol))
end
