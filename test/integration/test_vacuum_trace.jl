# A thin ramp-heated gasifying sample starts at exact vacuum. Trace production
# must cross the smooth positive-part region without a step-count explosion.
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

mesh = create_uniform_mesh(1e-4, 10)                  # 0.1 mm, 10 cells

heating = t -> 300.0 + t / 6                          # 10 K/min
bc_set = BoundaryConditionSet(
    thermal = (top = DirichletThermalBC(heating), bottom = DirichletThermalBC(heating)),
    mass = (top = AmbientGasBC(GasReservoir(P = 0.0, T = 300.0)), bottom = ImpermeableBC()),
)

problem =
    PyrolysisProblem(mesh = mesh, material = pmma, bc_set = bc_set, tspan = (0.0, 2400.0))



@testset "vacuum trace production resolves the gas guard ($energy)" for energy in (
    :temperature,
    :enthalpy,
)
    # At the reference enthalpy tolerance, use a nonlinear implicit method:
    # the caloric extension and an initially flat gas guard make the
    # linearly implicit method take very small startup steps. The ordinary
    # trajectory below still qualifies the default policy in both forms.
    reference_policy=energy===:enthalpy ? kencarp4_policy(StructuredKLU()) : nothing
    tight=solve(
        problem;
        energy_state = energy,
        integrator_policy = reference_policy,
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
        @test maximum(abs.(get_pressure(sol) .- get_pressure(tight)))<20ptol
    end
end
