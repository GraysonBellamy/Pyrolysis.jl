using Test, Pyrolysis, ForwardDiff
using Pyrolysis.Integration:
    PressureErrorNorm,
    pressure_error,
    kencarp4_policy,
    rodas5p_policy,
    fbdf_policy,
    default_linear_solver
using Pyrolysis.States: mass_index, energy_index, position_index, area_index
using Pyrolysis.Materials: R_GAS, DepletionLimiter

function pressure_fixture(; cells = 3)
    solid=SolidComponent(
        :matrix;
        ρ = 500.0,
        ρ_intrinsic = T->1000.0 / (1+1e-4*(T-300)),
        c = 1200.0,
        k = 0.2,
        λ = 0.0,
    )
    gas=GaseousComponent(:air; M = 0.029, c = 1000.0, k = 0.03)
    mat=Material(name = :sealed_heating, components = (solid, gas))
    PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, cells),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(1e4), bottom = AdiabaticBC()),
            mass = (top = ImpermeableBC(), bottom = ImpermeableBC()),
        ),
        T_initial = z->300.0,
        ξ_initial = [z->500.0, z->0.5*1e5*0.029/(R_GAS*300)],
        tspan = (0.0, 0.2),
    )
end
@testset "proposed-state pressure error" begin
    for energy in (:temperature, :enthalpy)
        cp=setup(
            pressure_fixture();
            energy_state = energy,
            pressure_tolerance = 0.01,
            reltol = 1e-10,
        )
        norm=PressureErrorNorm(cp)
        layout=cp.def.layout
        # A gas error with the same old mass scaling must be more restrictive
        # after the pore volume is reduced. The scale is recomputed at u.
        u=copy(cp.u0)
        e=zero(u)
        j=mass_index(layout, 1, 2)
        e[j]=1.0
        base=pressure_error(norm, u, u, e, cp.abstol, cp.reltol)
        u[mass_index(layout, 1, 1)]*=1.5
        if energy===:enthalpy
            # Keep T=300 by rebuilding the exact sensible inventory at T_REF.
            u[energy_index(layout, 1)]=cp.u0[energy_index(layout, 1)]*1.5
        end
        compressed=pressure_error(norm, u, u, e, cp.abstol, cp.reltol)
        @test compressed>1.9base
        # Thermal error alone changes pressure even with exact gas mass.
        e.=0
        e[energy_index(layout, 1)]=1e7
        @test pressure_error(norm, cp.u0, cp.u0, e, cp.abstol, cp.reltol)>0
    end
end
@testset "pressure control in supported steppers" begin
    for energy in (:temperature, :enthalpy),
        policy in (kencarp4_policy, rodas5p_policy, fbdf_policy)

        sol=solve(
            pressure_fixture();
            energy_state = energy,
            pressure_tolerance = 0.1,
            reltol = 1e-7,
            saveat = 0.05,
            integrator_policy = policy(default_linear_solver()),
        )
        @test completed(sol)
        @test sol.T[end, end]>300
        @test maximum(abs, (sol.m[2, :, end] ./ sol.m[2, :, 1] .- 1))*1e5 < 0.01
        @test sol.statistics.pressure_error_evaluations>0
        @test all(isfinite, sol.T)
    end
end

# A closed gas inventory compressed by a mass-conserving condensed reaction.
# m_A=m_A0 exp(-kt), V/V0=e^(-kt)+(ρ_A/ρ_B)(1-e^(-kt));
# the pore volume follows the two intrinsic densities independently of γ.
function compressing_pores()
    T=600.0
    P=1e5
    ρ=500.0
    k=0.5
    mat=Material(
        name = :pore_compression,
        components = (
            SolidComponent(:A; ρ = ρ, φ = 0.3, c = 1000.0, k = 0.2, λ = 0.0),
            SolidComponent(:B; ρ = 4ρ, φ = 0.05, c = 1000.0, k = 0.2, λ = 0.0),
            GaseousComponent(:gas; M = 0.03, c = 1000.0, k = 0.03),
        ),
        reactions = (Reaction(:compact, 1=>2; A = k, E = 0.0, h = 0.0),),
        depletion_limiter = DepletionLimiter(enabled = false),
    )
    problem=PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, 3),
        material = mat,
        bc_set = BoundaryConditionSet(),
        T_initial = z->T,
        ξ_initial = [z->ρ, z->0.0, z->0.3*P*0.03/(R_GAS*T)],
        tspan = (0.0, 8.0),
        output_times = collect(0.0:0.5:8.0),
    )
    exact=t->P*0.3/(0.3*exp(-k*t)+0.05/4*(1-exp(-k*t)))
    return problem, exact
end
@testset "pressure control through shrinking pore volume ($energy)" for energy in (
    :temperature,
    :enthalpy,
)
    problem, exact=compressing_pores()
    for policy in (kencarp4_policy, rodas5p_policy, fbdf_policy)
        # FBDF's first-order startup from an exactly absent product can need
        # a step below 1e-12 at this condensed-mass tolerance. This is its
        # ordinary state-error test, not a pressure-domain rejection.
        common=(
            energy_state = energy,
            use_ale = true,
            dtmin = 1e-16,
            integrator_policy = policy(default_linear_solver()),
        )
        tight=solve(
            problem;
            common...,
            pressure_tolerance = 0.001,
            reltol = 1e-10,
            abstol = 1e-12,
        )
        @test completed(tight)
        @test maximum(abs.(get_pressure(tight)[1, :] ./ exact.(tight.t) .- 1))<1e-7
        for ptol in (0.1, 0.01)
            sol=solve(
                problem;
                common...,
                pressure_tolerance = ptol,
                reltol = 1e-9,
                abstol = 1e-11,
            )
            @test completed(sol)
            @test sol.statistics.pressure_error_evaluations>0
            @test maximum(abs.(sol.T .- 600.0))<1e-6
            @test get_porosity(sol)[1, end]<0.1
            # Local budgets do not bound accumulated global error. This
            # explicit comparison qualifies the controller on this trajectory.
            @test maximum(abs.(get_pressure(sol) .- get_pressure(tight)))<20ptol
        end
    end
end
