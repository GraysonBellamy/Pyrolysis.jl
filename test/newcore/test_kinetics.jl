using Test
using ForwardDiff
using Pyrolysis
using Pyrolysis.Materials: R_GAS, STEFAN_BOLTZMANN, T_PROPERTY_MIN, T_PROPERTY_MAX
using Pyrolysis.Materials: DepletionLimiter, OneSidedLimiter
using Pyrolysis.Physics: reaction_rate, one_sided_factor, evaluate_reactions!
using Pyrolysis.Physics: ideal_gas_pressure
using Pyrolysis.Properties: effective_gas_viscosity, effective_conductivity
using Pyrolysis.States: trial_temperature
using Pyrolysis.Kernels: ReactionKernel, Contributions, contribute!, reset!
using Pyrolysis.Integration
using Pyrolysis.Integration: setup

_rxn(order) = Reaction(:r, 1 => 2; A = 1.0, E = 1000 * R_GAS, h = 0.0, n = order)
_rate(x, order, T; limiter = OneSidedLimiter()) =
    reaction_rate(_rxn(order), (x, zero(x)), trial_temperature(T); limiter = limiter)

# Trial concentrations off the physical domain, including the window
# ξ/base_scale ∈ (−745, −726) where the finite base underflows to zero before
# its slope does.
const _trial_concentrations = (-1e6, -1.0, -7.3e-6, -1e-9, 0.0, 1e-12, 1e-9, 1.0)

@testset "one-sided limiter: finite rates and derivatives at absent reactants" begin
    for order in (0.0, 0.3, 0.5, 1.0, 2.0), x in _trial_concentrations
        f(y) = _rate(y, order, 700.0)
        @test isfinite(f(x))
        @test isfinite(ForwardDiff.derivative(f, x))
        x <= 0 && @test f(x) == 0
        x <= 0 && @test ForwardDiff.derivative(f, x) == 0
        # Disabled: the finite base alone, finite at every real concentration.
        g(y) = _rate(y, order, 700.0; limiter = OneSidedLimiter(enabled = false))
        @test isfinite(g(x)) && g(x) >= 0
        @test isfinite(ForwardDiff.derivative(g, x))
    end
    # Finite at zero and negative trial temperatures through the trial clamp.
    for T in (-1e6, -1.0, 0.0, 1e6), x in (0.5, 10.0)
        @test isfinite(_rate(x, 1.0, T))
        @test isfinite(ForwardDiff.derivative(t -> _rate(x, 1.0, t), T))
    end
    # Negative control: a positive base alone consumes an absent reactant.
    @test _rate(0.0, 0.0, 700.0; limiter = OneSidedLimiter(enabled = false)) > 0
end

@testset "one-sided limiter: C¹ onset joined onto tanh" begin
    for join in (0.5, 1.0, 2.0)
        @test one_sided_factor(0.0, join) == 0
        @test ForwardDiff.derivative(x -> one_sided_factor(x, join), 0.0) == 0
        @test one_sided_factor(join, join) ≈ tanh(join) rtol = 1e-15
        @test ForwardDiff.derivative(x -> one_sided_factor(x, join), join - 1e-9) ≈
              1 / cosh(join)^2 rtol = 1e-6
        for x in range(0.0, join; length = 101)
            @test 0 <= one_sided_factor(x, join) <= 1
            @test ForwardDiff.derivative(y -> one_sided_factor(y, join), x) >= -1e-14
        end
        for x in (0.2join, 0.7join, 1.5join)
            h = 1e-6
            fd = (one_sided_factor(x + h, join) - one_sided_factor(x - h, join)) / 2h
            @test ForwardDiff.derivative(y -> one_sided_factor(y, join), x) ≈ fd rtol = 1e-8
        end
        for x in (join, 1.5join, 10.0, 20.0, 50.0)
            @test one_sided_factor(x, join) == (x < 20 ? tanh(x) : 1.0)
        end
    end
end

@testset "calibration guard: rates above the join equal the tanh limiter" begin
    for order in (0.0, 0.3, 0.5, 1.0, 2.0), threshold in (0.1, 1.0, 5.0)
        tanh_limiter = DepletionLimiter(threshold = threshold)
        one_sided = OneSidedLimiter(threshold = threshold)
        for x in (1.0, 1.5, 2.0, 20.0, 1e3) .* threshold, T in (300.0, 700.0, 1500.0)
            r_tanh = reaction_rate(_rxn(order), (x, 0.0), T; limiter = tanh_limiter)
            r_one = reaction_rate(_rxn(order), (x, 0.0), T; limiter = one_sided)
            @test r_one == r_tanh
        end
    end
    # Disabled limiters agree away from zero; the finite base is the only difference.
    for order in (0.5, 1.0), x in (1e-3, 1.0, 100.0)
        r_tanh = reaction_rate(
            _rxn(order),
            (x, 0.0),
            700.0;
            limiter = DepletionLimiter(enabled = false),
        )
        r_one = reaction_rate(
            _rxn(order),
            (x, 0.0),
            700.0;
            limiter = OneSidedLimiter(enabled = false),
        )
        @test r_one ≈ r_tanh rtol = 1e-12
    end
    @test OneSidedLimiter(DepletionLimiter(threshold = 0.3, enabled = false)) ==
          OneSidedLimiter(threshold = 0.3, enabled = false)
    @test OneSidedLimiter(OneSidedLimiter(join = 2.0)) == OneSidedLimiter(join = 2.0)
    @test_throws ArgumentError OneSidedLimiter(threshold = 0.0)
    @test_throws ArgumentError OneSidedLimiter(join = -1.0)
    @test_throws ArgumentError OneSidedLimiter(base_scale = 0.0)
end

# A decomposing slab whose reactant sits above the join everywhere: the
# kernel's sources equal those of the material's own tanh limiter.
@testset "reaction kernel: sources above the join equal the material's limiter" begin
    virgin = SolidComponent(:virgin, ρ = 1100.0, c = 1500.0, k = 0.2)
    char = SolidComponent(:char, ρ = 300.0, c = 1000.0, k = 0.1)
    gas = GaseousComponent(:gas, M = 0.03, c = 1100.0, k = 0.03, λ = 0.0)
    mat = Material(
        name = :slab,
        components = (virgin, char, gas),
        reactions = (
            Reaction(:pyro, 1 => (2, 3); A = 5e9, E = 1.4e5, h = 8e5, yields = (0.3, 0.7)),
            Reaction(:ox, 2 => 3; A = 1e3, E = 8e4, h = -2e6, n = 0.5),
        ),
        depletion_limiter = DepletionLimiter(threshold = 2.0),
    )
    problem = PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, 6, Val(3)),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(1e4), bottom = AdiabaticBC()),
        ),
        T_initial = z -> 600.0 + 2e4 * z,
        ξ_initial = [z -> 1100.0 - 5e4 * z, z -> 3.0 + 1e3 * z, z -> 0.5],
        tspan = (0.0, 1.0),
    )
    cp = setup(problem)
    kernel = only(k for k in cp.def.kernels if k isa ReactionKernel)
    @test kernel.limiter == OneSidedLimiter(threshold = 2.0)
    ds = cp.ws.derived
    c = Contributions{Float64}(3, 2, 6, 7, 7)
    reset!(c)
    Pyrolysis.States.derive!(ds, cp.u0, cp.def.layout, cp.def.setup)
    contribute!(c, kernel, ds, cp.def.setup.topology, mat, 0.0)
    rates, heats = zeros(2), zeros(2)
    for i = 1:6
        ξ = ntuple(j -> ds.ξ[j, i], 3)
        evaluate_reactions!(rates, heats, ξ, ds.T[i], mat)
        ω = (-rates[1], 0.3rates[1] - rates[2], 0.7rates[1] + rates[2])
        for j = 1:3
            @test c.S_mass[j, i] ≈ ω[j] * ds.V_eval[i] rtol = 1e-14
        end
        @test c.S_energy[i] ≈ -(heats[1] * rates[1] + heats[2] * rates[2]) * ds.V_eval[i] rtol =
            1e-14
    end
end

@testset "reaction kernel: a cell where every reaction lacks a reactant contributes exact zeros" begin
    virgin = SolidComponent(:virgin, ρ = 1100.0, c = 1500.0, k = 0.2)
    char = SolidComponent(:char, ρ = 300.0, c = 1000.0, k = 0.1)
    gas = GaseousComponent(:gas, M = 0.03, c = 1100.0, k = 0.03, λ = 0.0)
    reactions = (
        Reaction(:pyro, 1 => (2, 3); A = 5e9, E = 1.4e5, h = 8e5, yields = (0.3, 0.7)),
        Reaction(:ox, 2 => 3; A = 1e3, E = 8e4, h = -2e6, n = 0.5),
    )
    mat = Material(name = :slab, components = (virgin, char, gas), reactions = reactions)
    quenched = Pyrolysis.Kernels._all_quenched
    @test quenched((0.0, 0.0, 1.0), reactions)
    @test quenched((-1e-9, 0.0, 1.0), reactions)
    @test !quenched((0.0, 1e-12, 1.0), reactions)
    @test !quenched((1.0, 0.0, 1.0), reactions)
    @test quenched((1.0, 0.0), ())
    # Cells 1 and 2 hold neither reactant, cell 3 only char, the rest virgin.
    problem = PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, 6, Val(3)),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(1e4), bottom = AdiabaticBC()),
        ),
        T_initial = z -> 700.0,
        ξ_initial = [z -> z > 0.005 ? 1100.0 : 0.0, z -> z < 0.004 ? 0.0 : 5.0, z -> 0.5],
        tspan = (0.0, 1.0),
    )
    cp = setup(problem)
    kernel = only(k for k in cp.def.kernels if k isa ReactionKernel)
    ds = cp.ws.derived
    Pyrolysis.States.derive!(ds, cp.u0, cp.def.layout, cp.def.setup)
    c = Contributions{Float64}(3, 2, 6, 7, 7)
    reset!(c)
    contribute!(c, kernel, ds, cp.def.setup.topology, mat, 0.0)
    rates, heats = zeros(2), zeros(2)
    for i = 1:6
        ξ = ntuple(j -> ds.ξ[j, i], 3)
        evaluate_reactions!(rates, heats, ξ, ds.T_trial[i], mat, kernel.limiter)
        ω = (-rates[1], 0.3rates[1] - rates[2], 0.7rates[1] + rates[2])
        if quenched(ξ, reactions)
            @test all(iszero, rates)
            @test all(j -> c.S_mass[j, i] === 0.0, 1:3)
            @test c.S_energy[i] === 0.0
        else
            @test any(!iszero, rates)
            for j = 1:3
                @test c.S_mass[j, i] ≈ ω[j] * ds.V_eval[i] rtol = 1e-14
            end
        end
    end
    @test count(i -> quenched(ntuple(j -> ds.ξ[j, i], 3), reactions), 1:6) == 2
end

@testset "trial temperature is safe in representative kernels" begin
    gas = GaseousComponent(:gas, M = 0.03, c = 1100.0, k = 0.03, λ = 1e-5)
    solid = SolidComponent(:solid, ρ = 500.0, c = (900.0, 1.5), k = (0.2, 2e-4))
    mat = Material(name = :trial, components = (solid, gas), reactions = ())
    rxn = Reaction(:r, 1 => 2; A = 1e10, E = 1.5e5, h = 0.0, T_min = 400.0)
    kernels(T) = begin
        t = trial_temperature(T)
        [
            reaction_rate(rxn, (100.0 + zero(t), zero(t)), t),
            ideal_gas_pressure((100.0 + zero(t), 0.2 + zero(t)), t, 0.5, mat),
            effective_gas_viscosity(t, mat),
            effective_conductivity((100.0 + zero(t), 0.2 + zero(t)), t, mat),
            STEFAN_BOLTZMANN * t^4,
        ]
    end
    for T in (-1e6, -1.0, 0.0, 199.0, 200.0, 200.5, 300.0, 2999.5, 3000.0, 1e6)
        @test T_PROPERTY_MIN <= trial_temperature(T) <= T_PROPERTY_MAX
        @test all(isfinite, kernels(T))
        @test all(isfinite, ForwardDiff.jacobian(x -> kernels(x[1]), [T]))
    end
    @test trial_temperature(700.0) == 700.0
    @test ForwardDiff.derivative(trial_temperature, 700.0) == 1
    # Negative controls: the raw temperature breaks the Arrhenius factor.
    @test !isfinite(ForwardDiff.derivative(T -> exp(-1000 / T), 0.0))
    @test !isfinite(exp(-1000 / -1.0))
end

# One cell between Dirichlet faces holding an inert carrier, an absent reactant A
# and a present reactant C, both decomposing to a product: A must stay absent.
function _absent_reactant(order)
    comps = (
        SolidComponent(:carrier, ρ = 1000.0, c = 1000.0, k = 0.5),
        SolidComponent(:A, ρ = 1000.0, c = 1000.0, k = 0.5),
        SolidComponent(:C, ρ = 1000.0, c = 1000.0, k = 0.5),
        SolidComponent(:P, ρ = 1000.0, c = 1000.0, k = 0.5),
    )
    rxns = (
        Reaction(:a, 2 => 4; A = 2.0, E = 0.0, h = 0.0, n = order),
        Reaction(:c, 3 => 4; A = 2.0, E = 0.0, h = 0.0, n = 1.0),
    )
    mat = Material(name = :absent, components = comps, reactions = rxns)
    return PyrolysisProblem(
        mesh = create_uniform_mesh(1e-3, 1, Val(4)),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = DirichletThermalBC(700.0), bottom = DirichletThermalBC(700.0)),
        ),
        T_initial = z -> 700.0,
        ξ_initial = [z -> 500.0, z -> 0.0, z -> 500.0, z -> 0.0],
        tspan = (0.0, 2.0),
        output_times = [0.0, 1.0, 2.0],
    )
end

@testset "absent reactant through a reacting step: order $order, $energy form" for order in
                                                                                   (
        0.0,
        0.5,
        1.0,
    ),
    energy in (:temperature, :enthalpy)

    sol = Integration.solve(_absent_reactant(order); energy_state = energy)
    @test sol.retcode === :Success
    @test all(==(0.0), sol.ξ[2, 1, :])
    @test sol.ξ[3, 1, end] < 0.1 * 500.0
    @test sol.ξ[4, 1, end] > 0.9 * 500.0
    @test sum(sol.ξ[:, 1, end]) ≈ 1000.0 rtol = 1e-12
end

# Decay from a small initial concentration: the one-sided limiter brings the
# reactant to rest at zero instead of driving it negative.
@testset "approach to depletion: order $order" for order in (0.0, 0.5, 1.0)
    comps = (
        SolidComponent(:A, ρ = 1000.0, c = 1000.0, k = 0.5),
        SolidComponent(:B, ρ = 1000.0, c = 1000.0, k = 0.5),
    )
    mat = Material(
        name = :decay,
        components = comps,
        reactions = (Reaction(:d, 1 => 2; A = 10.0, E = 0.0, h = 0.0, n = order),),
    )
    problem = PyrolysisProblem(
        mesh = create_uniform_mesh(1e-3, 1, Val(2)),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = DirichletThermalBC(500.0), bottom = DirichletThermalBC(500.0)),
        ),
        T_initial = z -> 500.0,
        ξ_initial = [z -> 5.0, z -> 995.0],
        tspan = (0.0, 50.0),
        output_times = collect(range(0.0, 50.0; length = 26)),
    )
    sol = Integration.solve(problem; abstol = 1e-10, reltol = 1e-8)
    @test sol.retcode === :Success
    @test minimum(sol.ξ[1, 1, :]) >= -1e-9
    @test sol.ξ[1, 1, end] <= 1.0
    @test issorted(sol.ξ[1, 1, :]; rev = true)
    @test maximum(abs.(sol.ξ[1, 1, :] .+ sol.ξ[2, 1, :] .- 1000.0)) < 1e-9
end
