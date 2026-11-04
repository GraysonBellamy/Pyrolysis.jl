# A thin ramp-heated gasifying sample starts at exact vacuum, so its gas begins
# below the pressure budget. Trace production must pass through that range
# without a step-count explosion.
using Test, Pyrolysis

virgin = SolidComponent(:virgin; ρ = 1190.0, φ = 0.3, c = 1500.0, k = 0.20, λ = 1e-6)
gas = GaseousComponent(:MMA; M = 0.100, c = 1500.0, k = 0.02)
decomp = Reaction(:pyrolysis, :virgin => :MMA; A = 8.5e12, E = 188e3, h = 870e3)
pmma = Material(
    gas_transport = LumpedTransfer(
        endpoint = LumpedEndpoint(λ = 1e-6, κ = 0.0, max_skeletal_occupancy = 0.001),
    ),
    name = :PMMA,
    components = (virgin, gas),
    reactions = (decomp,),
)

mesh = uniform_mesh(1e-4, 10)                  # 0.1 mm, 10 cells

heating = t -> 300.0 + t / 6                          # 10 K/min
bc_set = BoundaryConditionSet(
    thermal = (top = DirichletThermalBC(heating), bottom = DirichletThermalBC(heating)),
    mass = (top = AmbientGasBC(GasReservoir(P = 0.0, T = 300.0)), bottom = ImpermeableBC()),
)

problem =
    PyrolysisProblem(mesh = mesh, material = pmma, bc_set = bc_set, tspan = (0.0, 2400.0))



@testset "vacuum trace production stays resolved ($energy)" for energy in
                                                                (:temperature, :enthalpy)
    tight=solve(
        problem;
        energy_state = energy,
        pressure_tolerance = 1e-5,
        reltol = 1e-10,
        abstol = 1e-12,
        maxiters = 20000,
        saveat = 20.0,
    )
    @test completed(tight)
    for ptol in (1.0, 0.1, 0.01)
        sol=solve(
            problem;
            energy_state = energy,
            pressure_tolerance = ptol,
            reltol = 1e-8,
            abstol = 1e-10,
            maxiters = 10000,
            saveat = 20.0,
        )
        @test completed(sol)
        @test sol.statistics.accepted<5000
        completed(sol) && completed(tight) || continue
        @test maximum(abs.(sol.T .- tight.T))<0.001
        @test maximum(abs.(sol.extras.total_mass .- tight.extras.total_mass))<1e-6*0.119
        @test maximum(abs.(sol.gas.pressure .- tight.gas.pressure))<20ptol
    end
end

# FBDF reuses its Jacobian and iterates Newton on each step, so this start, with
# its gas below the pressure budget, is where its stopping test has to hold.
@testset "vacuum trace production under FBDF ($energy)" for energy in
                                                            (:temperature, :enthalpy)
    sol=solve(
        problem;
        energy_state = energy,
        integrator_policy = fbdf_policy(StructuredKLU()),
        pressure_tolerance = 1.0,
        reltol = 1e-8,
        abstol = 1e-10,
        maxiters = 10000,
        saveat = 20.0,
    )
    @test completed(sol)
    @test sol.statistics.accepted<1000
    @test sol.statistics.nonlinear_failures<50
end

# The same start under dusty-gas transport with two product gases. Stages of
# early steps reach pore states far from the accepted ones, where a face's
# transport solve has no qualified root; each such step is rejected and
# retried, so the start completes.
@testset "dusty-gas trace production from vacuum" begin
    products = Reaction(
        :pyrolysis,
        1 => (2, 3);
        A = 8.5e12,
        E = 188e3,
        h = 870e3,
        yields = (0.7, 0.3),
    )
    dusty = Material(
        gas_transport = DustyGas(
            tortuosity = 2.0,
            D_ref = 2e-5,
            T_ref = 300.0,
            P_ref = 1e5,
            pore_diameter = 1e-6,
            endpoint = DustyEndpoint(
                pore_diameter = 1e-6,
                permeability = 0.0,
                max_skeletal_occupancy = 0.001,
            ),
        ),
        name = :PMMA,
        components = (virgin, gas, GaseousComponent(:CO2; M = 0.044, c = 1000.0, k = 0.02)),
        reactions = (products,),
    )
    vacuum = PyrolysisProblem(
        mesh = mesh,
        material = dusty,
        bc_set = bc_set,
        tspan = (0.0, 2400.0),
    )
    solutions = map((1.0, 0.01)) do ptol
        solve(
            vacuum;
            pressure_tolerance = ptol,
            reltol = 1e-8,
            abstol = 1e-10,
            maxiters = 10000,
            saveat = 20.0,
        )
    end
    @test all(completed, solutions)
    if all(completed, solutions)
        loose, tight = solutions
        @test maximum(abs.(loose.T .- tight.T))<0.001
        @test maximum(abs.(loose.gas.pressure .- tight.gas.pressure))<20
    end
end
