using Test
using ForwardDiff
using SparseArrays
using StaticArrays
using Pyrolysis
using Pyrolysis.Materials: DepletionLimiter, R_GAS
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.States: StateLayout, state_length
using Pyrolysis.Physics: interface_flux
using Pyrolysis.Properties: guarded_harmonic_mean
using Pyrolysis.Kernels: GasTransportKernel, contribute!, face_peclet
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

# A sealed, adiabatic, unreacting column of cells between `nodes`: a permeable
# solid of conductivity `k` holding the gases `gases` at the cell temperatures
# `T` and partial pressures `p` (one tuple per cell, in gas order).
function _gas_problem(
    T,
    p;
    nodes = 0.01 .* (0:length(T)),
    gases = (GaseousComponent(:gas; M = 0.028, c = 1000.0, k = 0.02),),
    model = LumpedTransfer(),
    k = 0.01,
    tspan = (0.0, 0.1),
)
    solid = SolidComponent(
        :solid;
        ρ = 1000.0,
        ρ_intrinsic = 2000.0,
        c = 1000.0,
        k = k,
        λ = 0.0,
        κ = 1e-13,
    )
    material = Material(
        name = :gas_cells,
        gas_transport = model,
        components = (solid, gases...),
        reactions = (),
    )
    cell(z) = clamp(searchsortedlast(nodes, z), 1, length(T))
    ξ0 = (500.0, ntuple(_ -> 0.0, length(gases))...)
    concentration(g) =
        z ->
            background_gas_concentration(material, 1 + g, ξ0, T[cell(z)]; P = p[cell(z)][g])
    return PyrolysisProblem(
        mesh = create_mesh_from_nodes(collect(nodes)),
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (top = ImpermeableBC(), bottom = ImpermeableBC()),
        ),
        T_initial = z -> T[cell(z)],
        ξ_initial = [z -> 500.0, (concentration(g) for g in eachindex(gases))...],
        tspan = tspan,
    )
end
_gas_cells(T, p; energy, kwargs...) =
    setup(_gas_problem(T, p; kwargs...); energy_state = energy)

# The residual at `u` and the gas transport kernel's own face fluxes there.
function _transport(cp, u = cp.u0)
    def = cp.def
    du = zeros(length(u))
    residual!(du, u, def, 0.0, cp.ws)
    ds, c = buffers(cp.ws, def, Float64)
    fill!(c.F_mass, 0)
    fill!(c.F_energy, 0)
    kernel = find_kernel(GasTransportKernel, def.kernels)
    contribute!(c, kernel, ds, def.setup.topology, def.setup.material, 0.0)
    return du, copy(c.F_mass), copy(c.F_energy)
end

# The partial derivatives of cell `i`'s energy state in its masses and its
# temperature: the enthalpies and heat capacity of the temperature row, or
# those of the extended caloric relation the inventory form stores.
function _caloric_partials(cp, i)
    ds, s = cp.ws.derived, cp.def.setup
    NC = n_components(cp.def.layout)
    m = SVector{NC}(ntuple(j -> ds.m[j, i], NC))
    energy_state(cp.def.layout) isa TemperatureState &&
        return SVector{NC}(ntuple(j -> ds.h[j, i], NC)), ds.cp_sum[i]
    E_m = ForwardDiff.gradient(x -> extended_energy(s.caloric, s.laws, ds.T[i], x, i), m)
    E_T = ForwardDiff.derivative(T -> extended_energy(s.caloric, s.laws, T, m, i), ds.T[i])
    return E_m, E_T
end

# Every cell's temperature rate in the residual `du`.
function _temperature_rates(cp, du)
    layout = cp.def.layout
    NC = n_components(layout)
    return map(1:n_cells(layout)) do i
        energy_state(layout) isa TemperatureState && return du[energy_index(layout, i)]
        E_m, E_T = _caloric_partials(cp, i)
        ṁ = SVector{NC}(ntuple(j -> du[mass_index(layout, i, j)], NC))
        (du[energy_index(layout, i)] - E_m' * ṁ) / E_T
    end
end

# Every cell's temperature rate from the gas transport alone, with the scale of
# its terms: the enthalpy each face carries in less the enthalpy its masses
# held in the cell.
function _gas_rates(cp, F_mass, q)
    topo, g = cp.def.setup.topology, cp.ws.derived.geometry
    rates = map(1:n_cells(topo)) do i
        E_m, C = _caloric_partials(cp, i)
        rate, scale = 0.0, 0.0
        for f in topo.cells[i].faces
            o, nb = topo.faces[f].cells
            nb == 0 && continue
            σ = o == i ? -1.0 : 1.0
            held = sum(F_mass[j, f] * E_m[j] for j in eachindex(E_m))
            rate += σ * g.face_area[f] * (q[f] - held)
            scale += g.face_area[f] * (abs(q[f]) + abs(held))
        end
        (rate / C, scale / C)
    end
    return first.(rates), last.(rates)
end

# The distance-weighted face temperature of the interior face `f`.
function _face_temperature(cp, f)
    ds, topo = cp.ws.derived, cp.def.setup.topology
    o, nb = topo.faces[f].cells
    d_o = center_to_face_distance(ds.geometry, topo, f, o)
    d_n = center_to_face_distance(ds.geometry, topo, f, nb)
    return (d_n * ds.T[o] + d_o * ds.T[nb]) / (d_o + d_n), d_o, d_n
end

const _HOT = ([600.0, 600.0, 300.0], [(3e5,), (2e5,), (1e5,)])
const _COLD = ([300.0, 300.0, 600.0], [(3e5,), (2e5,), (1e5,)])
const _MIRRORED = ([300.0, 600.0, 600.0], [(1e5,), (2e5,), (3e5,)])

# The conductivity a face adds for the net heat-capacity flux `C` of constant
# heat capacities: the least that keeps both cells' couplings nonnegative.
_added(C, d_o, d_n, k_f) = max(0.0, C * d_o - k_f, -C * d_n - k_f)

# The initial state with the gas of the first cell scaled by 1 + s, `k` its row.
_scaled(cp, k, s) = (u = copy(cp.u0); u[k] *= 1 + s; u)

# The interior face of a column of two cells.
_interior_face(cp) = only(
    f for f in eachindex(cp.def.setup.topology.faces) if
    cp.def.setup.topology.faces[f].cells[2] != 0
)

@testset "gas enthalpy at a plateau edge: no extremum ($energy form)" for energy in (
    :temperature,
    :enthalpy,
)
    # At a cell Péclet number of 45 the edge face is upwind: the gas reaches
    # the cell downstream at the plateau's enthalpy, and conduction no longer
    # reaches back across the face, so the plateau keeps its temperature and
    # only the cell downstream changes.
    for ((T, p), downstream) in ((_HOT, 3), (_COLD, 3), (_MIRRORED, 1))
        cp = _gas_cells(T, p; energy)
        du, F_mass, q = _transport(cp)
        ds, topo, g = cp.ws.derived, cp.def.setup.topology, cp.ws.derived.geometry
        rates = _temperature_rates(cp, du)
        scale = last(_gas_rates(cp, F_mass, q))
        for i in eachindex(T)
            T[i] == maximum(T) && @test rates[i] <= 1e-13 * scale[i]
            T[i] == minimum(T) && @test rates[i] >= -1e-13 * scale[i]
        end
        @test abs(rates[2]) <= 1e-13 * scale[2]
        f = only(
            topo.cells[downstream].faces[k] for
            k = 1:2 if topo.faces[topo.cells[downstream].faces[k]].cells[2] != 0
        )
        upstream = 2
        C = last(_caloric_partials(cp, downstream))
        inflow = abs(F_mass[2, f]) * (ds.h[2, upstream] - ds.h[2, downstream])
        @test rates[downstream] ≈ g.face_area[f] * inflow / C rtol = 1e-12
        T == _HOT[1] && @test rates[3] ≈ 2.690230604 atol = 1e-9
    end
    # Under a uniform pressure no gas moves and only conduction remains.
    cp = _gas_cells(_HOT[1], fill((2e5,), 3); energy)
    du, F_mass, q = _transport(cp)
    @test maximum(abs, _gas_rates(cp, F_mass, q)[1]) <= 1e-12
end

@testset "gas enthalpy: the kernel's species fluxes are the mass law's ($energy form)" for energy in
                                                                                           (
    :temperature,
    :enthalpy,
)
    T, p = _HOT
    cp = _gas_cells(T, p; energy)
    du, F_mass, q = _transport(cp)
    ds, topo = cp.ws.derived, cp.def.setup.topology
    kernel = find_kernel(GasTransportKernel, cp.def.kernels)
    for f = 1:n_faces(topo)
        o, nb = topo.faces[f].cells
        nb == 0 && continue
        d_o = center_to_face_distance(ds.geometry, topo, f, o)
        d_n = center_to_face_distance(ds.geometry, topo, f, nb)
        side(i) = (
            SVector(ds.p[2, i]),
            (K = SVector(ds.K[2, i]), t = ds.molecular_resistance[i], mob = ds.mob[i]),
        )
        (pl, cl), (pr, cr) = side(o), side(nb)
        n = interface_flux(pl, pr, ds.T_trial[o], ds.T_trial[nb], d_o, d_n, cl, cr).n
        @test F_mass[1, f] == 0
        @test F_mass[2, f] == n[1] * kernel.masses[1]
    end
end

@testset "gas enthalpy at a uniform temperature: no transport heating ($energy form)" for energy in
                                                                                          (
    :temperature,
    :enthalpy,
)
    # Two gases whose partial pressures vary from cell to cell of a stretched
    # column, all of them well above the caloric relation's mass scale, so its
    # partial derivatives are the enthalpies and the heat capacity.
    gases = (
        GaseousComponent(:a; M = 0.028, c = (1000.0, 0.3), k = 0.02),
        GaseousComponent(:b; M = 0.044, c = 1200.0, k = 0.02),
    )
    p = [(2e5, 0.5e5), (1e5, 1e5), (0.4e5, 2.5e5)]
    for model in values(_models)
        cp = _gas_cells(
            fill(600.0, 3),
            p;
            energy,
            gases,
            model,
            nodes = [0.0, 0.005, 0.02, 0.03],
        )
        du, F_mass, q = _transport(cp)
        ds, s = cp.ws.derived, cp.def.setup
        if energy === :enthalpy
            for i = 1:3
                m = ntuple(j -> ds.m[j, i], 3)
                @test all(j -> m[j] >= s.caloric.δ[i], 2:3)
                @test regularization(s.caloric, m, i) == 0
            end
        end
        @test count(>(0), abs.(F_mass)) >= 4
        gas, scale = _gas_rates(cp, F_mass, q)
        @test all(abs.(gas) .<= 1e-13 .* scale)
        @test maximum(abs, _temperature_rates(cp, du)) <= 1e-13 * maximum(scale)
    end
end

@testset "gas enthalpy with counterflowing species ($energy form)" for energy in (
    :temperature,
    :enthalpy,
)
    # Two gases of one molar mass at one total pressure trade places across a
    # temperature jump: equal and opposite mass fluxes, so no net mass moves,
    # but their heat capacities differ, so the net heat-capacity flux
    # C = Σⱼ Fⱼ cⱼ moves heat, from the cold cell to the hot one. In a
    # conducting matrix the face carries both gases at the face temperature
    # and the cold cell heats; in a nearly insulating one the face adds the
    # conductivity that keeps both couplings nonnegative, and the cold cell,
    # which C leaves, is no longer reached across the face.
    gases = (
        GaseousComponent(:a; M = 0.028, c = 1000.0, k = 0.02),
        GaseousComponent(:b; M = 0.028, c = 1500.0, k = 0.02),
    )
    for k in (0.01, 1e-4)
        cp = _gas_cells(
            [600.0, 300.0],
            [(1.8e5, 0.2e5), (0.2e5, 1.8e5)];
            energy,
            gases,
            model = _models.dusty,
            k,
        )
        du, F_mass, q = _transport(cp)
        ds, laws = cp.ws.derived, cp.def.setup.laws
        f = _interior_face(cp)
        o, nb = cp.def.setup.topology.faces[f].cells
        Fa, Fb = F_mass[2, f], F_mass[3, f]
        @test Fa > 0 > Fb
        @test abs(Fa + Fb) <= 1e-10 * abs(Fa)
        T_f, d_o, d_n = _face_temperature(cp, f)
        k_f = guarded_harmonic_mean(ds.k[o], ds.k[nb], d_o, d_n)
        D = _added(1000.0 * Fa + 1500.0 * Fb, d_o, d_n, k_f)
        @test (D == 0) == (k == 0.01)
        h = component_enthalpies(laws, T_f)
        @test q[f] ≈ Fa * h[2] + Fb * h[3] + D * (ds.T[o] - ds.T[nb]) / (d_o + d_n) rtol =
            1e-12
        @test abs(q[f]) > 1e-3 * abs(Fa) * h[2]
        rates = _temperature_rates(cp, du)
        scale = last(_gas_rates(cp, F_mass, q))
        @test rates[o] < 0
        if D > 0
            @test abs(rates[nb]) <= 1e-13 * scale[nb]
        else
            @test rates[nb] > 0
        end
    end
end

@testset "gas enthalpy: an exactly zero flux carries exactly nothing ($energy form)" for energy in
                                                                                         (
    :temperature,
    :enthalpy,
)
    # Empty pores across a temperature jump.
    for (gases, p, model) in (
        (
            (GaseousComponent(:gas; M = 0.028, c = 1000.0, k = 0.02),),
            [(0.0,), (0.0,)],
            LumpedTransfer(),
        ),
        (
            (
                GaseousComponent(:a; M = 0.028, c = 1000.0, k = 0.02),
                GaseousComponent(:b; M = 0.044, c = 1500.0, k = 0.02),
            ),
            [(0.0, 0.0), (0.0, 0.0)],
            _models.dusty,
        ),
    )
        cp = _gas_cells([600.0, 300.0], p; energy, gases, model)
        du, F_mass, q = _transport(cp)
        @test all(iszero, F_mass)
        @test all(iszero, q)
    end
end

@testset "gas enthalpy on unequal cells, resolved and upwind ($energy form)" for energy in (
    :temperature,
    :enthalpy,
)
    # Cells of 5, 15 and 10 mm, with the flow both ways. The face carries the
    # enthalpy at the distance-weighted face temperature plus
    # D (T_o − T_n)/(d_o + d_n), D = max(0, F c d_o − k_f, −F c d_n − k_f):
    # nothing more in a conducting matrix (k = 1), the excess of the upstream
    # half cell's heat-capacity flux over conduction in a poorly conducting one.
    for k in (1.0, 0.01),
        (T, p, nodes) in
        ((_HOT..., [0.0, 0.005, 0.02, 0.03]), (_MIRRORED..., [0.0, 0.01, 0.025, 0.03]))

        cp = _gas_cells(T, p; energy, nodes, k)
        du, F_mass, q = _transport(cp)
        ds, topo, laws = cp.ws.derived, cp.def.setup.topology, cp.def.setup.laws
        for f = 1:n_faces(topo)
            o, nb = topo.faces[f].cells
            nb == 0 && continue
            @test all(isfinite, F_mass[:, f])
            @test sign(F_mass[2, f]) == sign(ds.p[2, o] - ds.p[2, nb])
            T_f, d_o, d_n = _face_temperature(cp, f)
            interpolated = sum(F_mass[:, f] .* component_enthalpies(laws, T_f))
            k_f = guarded_harmonic_mean(ds.k[o], ds.k[nb], d_o, d_n)
            D = _added(1000.0 * F_mass[2, f], d_o, d_n, k_f)
            if k == 1.0
                @test D == 0
                @test q[f] == interpolated
            else
                @test D > 0
                @test q[f] ≈ interpolated + D * (ds.T[o] - ds.T[nb]) / (d_o + d_n) rtol =
                    1e-12
            end
        end
        rates = _temperature_rates(cp, du)
        scale = last(_gas_rates(cp, F_mass, q))
        for i in eachindex(T)
            T[i] == maximum(T) && @test rates[i] <= 1e-13 * scale[i]
            T[i] == minimum(T) && @test rates[i] >= -1e-13 * scale[i]
        end
    end
end

# A heat capacity that peaks a hundredfold over 60 μK.
const _PEAKED = TableProperty((600.00030, 600.00033, 600.00036), (1000.0, 100000.0, 1000.0))

@testset "gas enthalpy across the knots of a tabulated heat capacity ($energy form)" for energy in
                                                                                         (
    :temperature,
    :enthalpy,
)
    # Two cells 0.5 mK apart straddle the peak, and the face is upwind: the
    # added conductivity follows the enthalpy secants across the knots, so
    # the hot cell is not heated.
    gases = (GaseousComponent(:gas; M = 0.028, c = _PEAKED, k = 0.02),)
    cp = _gas_cells([600.0005, 600.0], [(3e5,), (1e5,)]; energy, gases)
    du, F_mass, q = _transport(cp)
    ds, topo, laws = cp.ws.derived, cp.def.setup.topology, cp.def.setup.laws
    f = _interior_face(cp)
    o, nb = topo.faces[f].cells
    T_f, d_o, d_n = _face_temperature(cp, f)
    h(T) = component_enthalpies(laws, T)[2]
    a(T) = (d_o + d_n) * F_mass[2, f] * (h(T_f) - h(T)) / (ds.T[nb] - ds.T[o])
    k_f = guarded_harmonic_mean(ds.k[o], ds.k[nb], d_o, d_n)
    D =
        face_peclet(find_kernel(GasTransportKernel, cp.def.kernels), ds, topo, f).diffusion *
        k_f
    @test D ≈ max(0.0, a(ds.T[o]) - k_f, a(ds.T[nb]) - k_f) rtol = 1e-9
    @test D > 0
    rates = _temperature_rates(cp, du)
    scale = last(_gas_rates(cp, F_mass, q))
    @test rates[o] <= 1e-13 * scale[o]
    @test rates[nb] >= -1e-13 * scale[nb]
end

@testset "gas enthalpy: the added conductivity is continuous across the secant limit" begin
    # Temperature differences just under and just over the resolution below
    # which the secants give way to their limit, with the face upwind: for a
    # heat capacity linear in T, and for the peaked table with a knot inside
    # the owner's half interval.
    τ = Pyrolysis.Kernels.SECANT_RESOLUTION
    for (c, T_o) in (((1000.0, 0.4), 600.0), (_PEAKED, 600.00033 + 1.5e-9))
        gases = (GaseousComponent(:gas; M = 0.028, c, k = 0.02),)
        added = map((1 - 1e-3, 1 + 1e-3)) do r
            cp = _gas_cells(
                [T_o, T_o - r * τ * T_o],
                [(3e5,), (1e5,)];
                energy = :temperature,
                gases,
            )
            _transport(cp)
            kernel = find_kernel(GasTransportKernel, cp.def.kernels)
            face_peclet(kernel, cp.ws.derived, cp.def.setup.topology, _interior_face(cp)).diffusion
        end
        @test all(>(0), added)
        @test added[1] ≈ added[2] rtol = 1e-3
    end
end

@testset "gas enthalpy conserves energy across each face ($energy form)" for energy in (
    :temperature,
    :enthalpy,
)
    # Two cells: what one gains the other loses. Three cells: a sealed,
    # adiabatic column holds its total enthalpy.
    for (T, p) in (([600.0, 300.0], [(3e5,), (1e5,)]), _HOT, _COLD)
        cp = _gas_cells(T, p; energy)
        du, F_mass, q = _transport(cp)
        ds, layout = cp.ws.derived, cp.def.layout
        rates = _temperature_rates(cp, du)
        H = map(eachindex(T)) do i
            E_m, C = _caloric_partials(cp, i)
            C * rates[i] + sum(E_m[j] * du[mass_index(layout, i, j)] for j = 1:2)
        end
        scale = sum(abs, cp.ws.derived.geometry.face_area) * maximum(abs, q)
        @test abs(sum(H)) <= 1e-13 * scale
        length(T) == 2 && @test H[1] ≈ -H[2] rtol = 1e-12
    end
end

@testset "gas enthalpy derivatives across a flow reversal and the upwind switch ($energy form)" for energy in
                                                                                                    (
    :temperature,
    :enthalpy,
)
    # The gas of the hot cell, scaled by 1 + s, sets the flow. Dense and
    # colored Jacobians agree with central differences on either side of each
    # switch, and the face's added conductivity is zero on one side of the
    # upwind switch and positive on the other.
    function check(cp, k, s)
        def, ws = cp.def, cp.ws
        state = _scaled(cp, k, s)
        residual(x) = (du = zeros(length(x)); residual!(du, x, def, 0.0, ws); du)
        J = jac_prototype(cp.jacobian)
        Jd = dense_jacobian(def, ws, state, 0.0)
        jacobian!(J, cp.jacobian, state, 0.0)
        @test compare_entries(Matrix(J), Jd; rtol = 1e-11).passed
        h = 1e-7
        up, down = copy(state), copy(state)
        up[k] += h * cp.u0[k]
        down[k] -= h * cp.u0[k]
        fd = (residual(up) - residual(down)) / (2h)
        ad = Jd[:, k] * cp.u0[k]
        @test maximum(abs, ad - fd) <= 1e-6 * maximum(abs, ad)
        return ad
    end
    # A resolved reversal: the carried enthalpy is h(T_f) on both sides, so
    # only the mass law's own upwind switch acts, and each energy row's
    # derivative jumps by its mass rate's jump times h(T_f), less the cell's
    # own enthalpy in the temperature form.
    cp = _gas_cells([600.0, 300.0], [(2e5,), (2e5,)]; energy)
    layout = cp.def.layout
    k = mass_index(layout, 1, 2)
    slopes = map((-1e-3, 1e-3)) do s
        @test sign(_transport(cp, _scaled(cp, k, s))[2][2, _interior_face(cp)]) == sign(s)
        check(cp, k, s)
    end
    _transport(cp)
    ds = cp.ws.derived
    h_f = component_enthalpies(
        cp.def.setup.laws,
        first(_face_temperature(cp, _interior_face(cp))),
    )[2]
    for i = 1:2
        below, above = (slope[mass_index(layout, i, 2)] for slope in slopes)
        jump = (above - below) * h_f
        energy === :temperature &&
            (jump = (above - below) * (h_f - ds.h[2, i]) / ds.cp_sum[i])
        e = energy_index(layout, i)
        @test slopes[2][e] - slopes[1][e] ≈ jump rtol = 0.02
    end
    # The upwind switch: from a pressure ratio of three, where the face is
    # upwind, down to nearly equal pressures, where it is resolved.
    cp = _gas_cells([600.0, 300.0], [(3e5,), (1e5,)]; energy)
    k = mass_index(cp.def.layout, 1, 2)
    kernel = find_kernel(GasTransportKernel, cp.def.kernels)
    f = _interior_face(cp)
    function added(s)
        _transport(cp, _scaled(cp, k, s))
        return face_peclet(kernel, cp.ws.derived, cp.def.setup.topology, f).diffusion
    end
    lo, hi = -0.66, 0.0
    @test added(lo) == 0 && added(hi) > 0
    for _ = 1:50
        mid = (lo + hi) / 2
        added(mid) > 0 ? (hi = mid) : (lo = mid)
    end
    for s in (lo - 2e-3, hi + 2e-3, 0.0)
        @test (added(s) > 0) == (s > lo)
        check(cp, k, s)
    end
end

@testset "gas enthalpy: face Péclet number and added conductivity" begin
    # On the hot plateau the hot/cold face has a cell Péclet number of 44.9,
    # and the face adds max(0, F c d_up/k_f − 1) face conductivities: Pe/2 − 1
    # on equal cells. In a conducting matrix nothing is added.
    for k in (0.01, 1.0), nodes in (0.01 .* (0:3), [0.0, 0.005, 0.02, 0.03])
        cp = _gas_cells(_HOT...; energy = :temperature, nodes, k)
        du, F_mass, q = _transport(cp)
        ds, topo = cp.ws.derived, cp.def.setup.topology
        kernel = find_kernel(GasTransportKernel, cp.def.kernels)
        for f = 1:n_faces(topo)
            o, nb = topo.faces[f].cells
            r = face_peclet(kernel, ds, topo, f)
            if nb == 0
                @test isnan(r.peclet) && isnan(r.diffusion)
                continue
            end
            d_o = center_to_face_distance(ds.geometry, topo, f, o)
            d_n = center_to_face_distance(ds.geometry, topo, f, nb)
            k_f = guarded_harmonic_mean(ds.k[o], ds.k[nb], d_o, d_n)
            C = 1000.0 * F_mass[2, f]
            @test r.peclet ≈ abs(C) * (d_o + d_n) / k_f rtol = 1e-14
            @test r.diffusion ≈ _added(C, d_o, d_n, k_f) / k_f rtol = 1e-12
            k == 1.0 && @test r.diffusion == 0
            if k == 0.01 && nodes isa AbstractRange && Set((o, nb)) == Set((2, 3))
                @test r.peclet ≈ 44.9127 atol = 1e-4
                @test r.diffusion ≈ r.peclet / 2 - 1 rtol = 1e-12
            end
        end
    end
end

@testset "gas enthalpy over a solve: the plateau stays within its extremes ($energy form)" for energy in
                                                                                               (
    :temperature,
    :enthalpy,
)
    # The sealed hot plateau over the second in which its pressures equalize:
    # gas pours from the plateau into the cold cell, and no temperature leaves
    # the initial range by more than the integration tolerance allows.
    for (T, p) in (_HOT, _COLD)
        problem = _gas_problem(T, p; tspan = (0.0, 1.0))
        sol = Pyrolysis.Integration.solve(
            problem;
            energy_state = energy,
            reltol = 1e-8,
            abstol = 1e-10,
            saveat = 0.01,
        )
        @test sol.retcode === :Success
        pressure = get_pressure(sol)
        @test maximum(pressure[:, end]) - minimum(pressure[:, end]) <
              0.1 * (maximum(pressure[:, 1]) - minimum(pressure[:, 1]))
        @test maximum(sol.T) <= 600 + 1e-4
        @test minimum(sol.T) >= 300 - 1e-4
    end
end

@testset "gas enthalpy: warmed residual allocates nothing ($energy form)" for energy in (
    :temperature,
    :enthalpy,
)
    cp = _gas_cells(_HOT...; energy)
    du = zeros(length(cp.u0))
    residual!(du, cp.u0, cp.def, 0.0, cp.ws)
    @test (@allocated residual!(du, cp.u0, cp.def, 0.0, cp.ws)) == 0
end
