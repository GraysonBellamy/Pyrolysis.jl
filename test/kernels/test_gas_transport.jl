using Test
using ForwardDiff
using SparseArrays
using Pyrolysis
using Pyrolysis.Materials: DepletionLimiter, R_GAS
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.States: StateLayout, state_length
using Pyrolysis.Kernels: GasTransportKernel
using Pyrolysis.Kernels: ThermalBoundaryKernel, find_kernel
using Pyrolysis.Assembly
using Pyrolysis.Assembly: residual!
using Pyrolysis.ColoredJacobian
using Pyrolysis.Integration: setup

# A decomposing slab with a product gas and a background gas: temperature-
# dependent properties, matrix gas-transfer coefficients and, for
# coupled transport, an impermeable virgin solid and a permeable char, so
# the char-free cells near the back carry zero mobility and the face
# porosities differ across the char front.
function _gas_slab(n; model = LumpedTransfer(), darcy::Bool, thermal, mass)
    virgin = SolidComponent(
        :virgin,
        ρ = 600.0,
        c = (1200.0, 2.0),
        k = (0.15, 1e-4),
        λ = 2e-6,
        ρ_intrinsic = 1400.0,
    )
    char = SolidComponent(
        :char,
        ρ = 150.0,
        c = 900.0,
        k = 0.08,
        λ = 1e-5,
        κ = darcy ? 1e-11 : 0.0,
        ρ_intrinsic = 1800.0,
    )
    gas = GaseousComponent(:gas, M = 0.03, c = (1000.0, 0.4), k = 0.03)
    air = GaseousComponent(:air, M = 0.029, c = 1006.0, k = 0.026)
    mat = Material(
        name = :gas_slab,
        gas_transport = model,
        mixing = MaterialMixing(κ = MixingSpec(PARALLEL)),
        components = (virgin, char, gas, air),
        reactions = (
            Reaction(
                :pyro,
                1 => (2, 3);
                A = 5e8,
                E = 1.3e5,
                h = 5e5,
                yields = (0.25, 0.75),
            ),
        ),
        depletion_limiter = DepletionLimiter(threshold = 1.0),
    )
    L = 0.01
    return PyrolysisProblem(
        mesh = create_stretched_mesh(L, n, 1.15),
        material = mat,
        bc_set = BoundaryConditionSet(thermal = thermal, mass = mass),
        T_initial = z -> 500.0 + 3e4 * z,
        ξ_initial = [
            z -> 600.0 * (1 - 0.6 * max(0.0, z / L - 0.4) / 0.6),
            z -> 150.0 * max(0.0, z / L - 0.4) / 0.6,
            z -> 0.3 + 40.0 * z,
            z -> 0.5 - 20.0 * z,
        ],
        tspan = (0.0, 10.0),
    )
end

# States on and off the physical manifold.
function _states(cp)
    u0 = cp.u0
    layout = cp.def.layout
    hot = copy(u0)
    for i = 1:n_cells(layout)
        hot[energy_index(layout, i)] =
            energy_state(layout) isa TemperatureState ? 700.0 + 20i :
            u0[energy_index(layout, i)] * 3
    end
    perturbed = u0 .* (1 .+ 0.2 .* sin.(1:length(u0)))
    zero_state = zeros(length(u0))
    negative = -abs.(u0)
    return (; u0, hot, perturbed, zero_state, negative)
end


const _models=(
    lumped = LumpedTransfer(),
    dusty = DustyGas(
        tortuosity = 2.0,
        D_ref = 2e-5,
        T_ref = 300.0,
        P_ref = 1e5,
        pore_diameter = 1e-6,
    ),
)
const _exchange=(
    vacuum = AmbientGasBC(GasReservoir(P = 0.0, T = 300.0)),
    atmosphere = AmbientGasBC(
        GasReservoir(P = 1e5, T = 300.0, mole_fractions = (air = 1.0,)),
    ),
    film = AmbientGasBC(
        GasReservoir(P = 1e5, T = 300.0, mole_fractions = (air = 1.0,));
        film = MassTransferFilm(1e-3),
    ),
    hydraulic = ConvectivePressureBC(
        GasReservoir(P = 1e5, T = 300.0, mole_fractions = (air = 1.0,));
        resistance = 1e6,
        flow = :darcy_velocity,
    ),
    gradient = NeumannPressureBC(
        gradient = 0.0,
        T = 300.0,
        mole_fractions = (gas = 0.2, air = 0.8),
    ),
    prescribed = MassFluxBC(
        fluxes = (gas = 1e-4, air = -1e-5),
        inflow_temperature = (air = 300.0,),
    ),
)
@testset "coupled transport Jacobian ($energy, $modelname, $boundary)" for energy in (
        :temperature,
        :enthalpy,
    ),
    modelname in keys(_models),
    boundary in keys(_exchange)

    println("qualifying gas boundary: ", energy, ", ", modelname, ", ", boundary)
    flush(stdout)
    cp=setup(
        _gas_slab(
            6;
            model = _models[modelname],
            darcy = true,
            thermal = (
                top = ConvectiveBC(h = 15.0, T_inf = 300.0, blowing = CouetteFilm()),
                bottom = AdiabaticBC(),
            ),
            mass = (top = _exchange[boundary], bottom = ImpermeableBC()),
        );
        energy_state = energy,
    )
    def, ws, cj=cp.def, cp.ws, cp.jacobian
    du=zero(cp.u0)
    J=jac_prototype(cj)
    for (state, u) in pairs(_states(cp))
        @testset "$state" begin
            residual!(du, u, def, 0.0, ws)
            @test all(isfinite, du)
            Jd=dense_jacobian(def, ws, u, 0.0)
            @test all(isfinite, Jd)
            @test isempty(missing_entries(Jd, cj.pattern))
            jacobian!(J, cj, u, 0.0)
            @test compare_entries(Matrix(J), Jd; rtol = 1e-11).passed
        end
    end
    residual!(du, cp.u0, def, 0.0, ws)
    @test (@allocated residual!(du, cp.u0, def, 0.0, ws))==0
end
# A sealed, isothermal coupled column with a pressure peak in the middle:
# the gas flows out of the peak both ways and stays in the column.
function _sealed_column()
    L = 0.01
    solid = SolidComponent(
        :solid,
        ρ = 1000.0,
        c = 1500.0,
        k = 0.3,
        κ = 1e-11,
        λ = 0.0,
        ρ_intrinsic = 1500.0,
    )
    gas = GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.025)
    return PyrolysisProblem(
        mesh = create_uniform_mesh(L, 20),
        material = Material(
            name = :sealed_column,
            components = (solid, gas),
            reactions = (),
        ),
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (top = ImpermeableBC(), bottom = ImpermeableBC()),
        ),
        T_initial = z -> 400.0,
        ξ_initial = [z -> 800.0, z -> 0.1 + 2.0 * exp(-((z / L - 0.5) / 0.1)^2)],
        tspan = (0.0, 1.0),
        output_times = [0.0, 0.5, 1.0],
    )
end

@testset "gas transport: a sealed coupled column keeps its gas and its solid ($energy form)" for energy in
                                                                                                 (
    :temperature,
    :enthalpy,
)
    problem = _sealed_column()
    cp = setup(problem; energy_state = energy)
    layout, topo = cp.def.layout, cp.def.setup.topology
    @test cp.def.setup.gas isa LumpedTransfer
    du = zeros(length(cp.u0))
    residual!(du, cp.u0, cp.def, 0.0, cp.ws)
    ds, c = buffers(cp.ws, cp.def, Float64)
    # A permeable solid under a steep pressure gradient: the gas flows.
    @test all(>(0), ds.mob)
    @test maximum(ds.P) > 10 * minimum(ds.P)
    @test maximum(abs, c.F_mass) > 0
    # The sealed faces carry nothing, the solid does not move, and what the
    # gas rows lose and gain sums to zero.
    for tag in (:top, :bottom)
        @test all(iszero, c.F_mass[:, only(boundary_faces(topo, tag))])
    end
    @test all(i -> du[mass_index(layout, i, 1)] == 0, 1:n_cells(layout))
    rates = [du[mass_index(layout, i, 2)] for i = 1:n_cells(layout)]
    @test abs(sum(rates)) <= 1e-12 * sum(abs, rates)
    # Over a solve the gas inventory stays where it started.
    sol = Pyrolysis.Integration.solve(
        problem;
        energy_state = energy,
        reltol = 1e-8,
        abstol = 1e-10,
    )
    @test sol.retcode === :Success
    gas = [sum(sol.m[2, :, k]) for k in eachindex(sol.t)]
    @test maximum(abs, gas .- gas[1]) <= 1e-12 * gas[1]
end

# Two sealed isothermal cells at 600 K with unequal gas concentrations, the gas
# in a matrix carrying the transfer law `λ`.
function _sealed_pair(λ)
    solid = SolidComponent(:solid; ρ = 600.0, c = 1000.0, k = 0.1, φ = 0.4, λ = λ)
    gas = GaseousComponent(:gas; M = 0.032, c = 1000.0, k = 0.03)
    return PyrolysisProblem(
        mesh = create_uniform_mesh(1.0, 2),
        material = Material(name = :sealed_pair, components = (solid, gas), reactions = ()),
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (top = ImpermeableBC(), bottom = ImpermeableBC()),
        ),
        T_initial = z -> 600.0,
        ξ_initial = [z -> 600.0, z -> z < 0.5 ? 1.0 : 2.0],
        tspan = (0.0, 1.0),
    )
end

@testset "gas transport: a transfer law zero at sampled temperatures still transports" begin
    rates = map((
        TableProperty((200.0, 300.0, 600.0, 1000.0, 3000.0), (0.0, 0.0, 1e-5, 0.0, 0.0)),
        1e-5,
    )) do λ
        cp = setup(_sealed_pair(λ))
        @test cp.def.setup.gas isa LumpedTransfer
        du = zeros(length(cp.u0))
        residual!(du, cp.u0, cp.def, 0.0, cp.ws)
        [du[mass_index(cp.def.layout, i, 2)] for i = 1:2]
    end
    @test rates[1] ≈ rates[2] rtol = 1e-12
    @test rates[1][1] > 0
    @test setup(_sealed_pair(0.0)).def.setup.gas isa LumpedTransfer
end

# The Carman–Kozeny permeability reads the state mapping's smoothly floored
# porosity, so a trial state whose skeleton overfills its cell keeps the
# positive mobility of the floor, joined C¹ to the physical states.
@testset "Carman–Kozeny permeability at the smooth porosity" begin
    ℓ = 1e-4
    solid =
        SolidComponent(:s, ρ = 500.0, c = 1500.0, k = 0.2, λ = 1e-5, ρ_intrinsic = 1000.0)
    gas = GaseousComponent(:gas, M = 0.03, c = 1100.0, k = 0.05)
    mat = Material(;
        name = :ck,
        components = (solid, gas),
        reactions = (),
        mixing = MaterialMixing(κ = MixingSpec(CARMAN_KOZENY; length_scale = ℓ)),
    )
    T = 400.0
    for ξ in ((500.0, 0.1), (900.0, 0.1))
        φ = 1 - ξ[1] / 1000.0
        @test Pyrolysis.Properties.effective_permeability(ξ, T, mat, φ) ≈
              min(ℓ^2 * φ^3 / (180 * (1 - φ)^2), ℓ^2) rtol = 1e-14
    end
    # Overfilled: the floored porosity keeps a positive permeability.
    κ(s) = Pyrolysis.Properties.effective_permeability(
        (s, oftype(s, 0.1)),
        T,
        mat,
        Pyrolysis.States.smooth_floor(1 - s / 1000.0, 1e-12),
    )
    @test κ(1100.0) > 0
    @test κ(1100.0) == κ(1000.0 * (1 + 1e-9))
    @test ForwardDiff.derivative(κ, 1100.0) == 0
    # Inside the floor's band the permeability still falls with the solid.
    band = 1000.0 * (1 - 1.5e-12)
    @test ForwardDiff.derivative(κ, band) < 0
    @test isfinite(ForwardDiff.derivative(κ, band))
end
