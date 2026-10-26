#=
Pore gas transport through a charring slab.

The same specimen is solved with superficial lumped transfer and a dusty-gas
friction model. Both use the condensed matrix permeability, skeletal pore
volume, an atmospheric air reservoir, and conservative species/enthalpy fluxes.
The two wall-transfer closures have independent constitutive data.

Run with PYROLYSIS_SMOKE=1 for a short, coarse example.
=#
using Pyrolysis
using Printf

const SMOKE = get(ENV, "PYROLYSIS_SMOKE", "0") == "1"

function transport_example(model)
    material = Material(
        name = :CharringSlab,
        gas_transport = model,
        components = (
            SolidComponent(
                :virgin;
                ρ = 500.0,
                ρ_intrinsic = 1500.0,
                c = 2000.0,
                k = 0.15,
                ε = 0.9,
                λ = 5e-7,
                κ = 1e-15,
            ),
            SolidComponent(
                :char;
                ρ = 150.0,
                ρ_intrinsic = 1950.0,
                c = 1500.0,
                k = 0.10,
                ε = 0.95,
                λ = 1e-5,
                κ = 1e-12,
            ),
            GaseousComponent(:volatiles; M = 0.044, c = 1200.0, k = 0.03),
            GaseousComponent(:air; M = 0.02897, c = 1006.0, k = 0.026),
        ),
        reactions = (
            Reaction(
                :pyrolysis,
                1=>(2, 3);
                yields = (0.2, 0.8),
                A = 1e13,
                E = 150e3,
                h = 300e3,
            ),
        ),
        mixing = MaterialMixing(κ = MixingSpec(WEIGHTED, 0.5)),
    )
    ξ_air = background_gas_concentration(material, :air, (500.0, 0.0, 0.0, 0.0), 300.0)
    reservoir = GasReservoir(P = 101325.0, T = 300.0, mole_fractions = (air = 1.0,))
    bcs = BoundaryConditionSet(
        thermal = (
            top = HeatFluxBC(75e3)+RadiativeBC(ε = 0.9, T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
        mass = (
            top = AmbientGasBC(reservoir; film = MassTransferFilm(0.012)),
            bottom = ImpermeableBC(),
        ),
    )
    final_time = SMOKE ? 10.0 : 120.0
    problem = PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, SMOKE ? 10 : 50),
        material = material,
        bc_set = bcs,
        T_initial = z->300.0,
        ξ_initial = [z->500.0, z->0.0, z->0.0, z->ξ_air],
        tspan = (0.0, final_time),
        output_times = collect(range(0.0, final_time; length = SMOKE ? 11 : 121)),
    )
    recorder = ConservationRecorder()
    solution = solve(
        problem;
        use_ale = true,
        min_thickness = .5e-3,
        pressure_tolerance = 0.1,
        reltol = 1e-7,
        callback = recorder,
    )
    completed(solution) || error("transport example ended with $(solution.outcome)")
    pressure = get_pressure(solution)
    exchange = get_boundary_exchange(solution, :top)
    product = findfirst(==(:volatiles), solution.gas.species)
    ledger = gas_species_audit(solution, recorder)
    @printf(
        "%s: peak pressure %.1f Pa, peak product escape %.4g kg/(m² s)\n",
        nameof(typeof(model)),
        maximum(filter(isfinite, pressure)),
        maximum(exchange.fluxes[product, :])
    )
    println("Species balance defects [kg]: ", map(x->x.defect, ledger.species))
    return solution
end

lumped_solution = transport_example(LumpedTransfer())
dusty_solution = transport_example(
    DustyGas(
        tortuosity = 2.0,
        D_ref = 2e-5,
        T_ref = 300.0,
        P_ref = 101325.0,
        pore_diameter = 1e-6,
    ),
)
