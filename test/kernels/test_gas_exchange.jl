using Test, Pyrolysis, StaticArrays, ForwardDiff
using Pyrolysis.Materials: R_GAS
using Pyrolysis.Assembly: residual!, buffers
using Pyrolysis.States: mass_index, energy_index
using Pyrolysis.Meshes: n_faces

function gas_problem(;
    model = LumpedTransfer(),
    mass = ImpermeableBC(),
    thermal = AdiabaticBC(),
    T = 500.0,
    P = 1e5,
    x = (0.7, 0.3),
    cells = 2,
    κ = 0.0,
)
    solid=SolidComponent(
        :matrix;
        ρ = 500.0,
        ρ_intrinsic = 1000.0,
        c = 1000.0,
        k = 0.2,
        λ = 1e-5,
        κ = κ,
    )
    a=GaseousComponent(:a; M = 0.028, c = 1000.0, k = 0.03)
    b=GaseousComponent(:b; M = 0.044, c = 1200.0, k = 0.03)
    mat=Material(name = :gas_exchange, components = (solid, a, b), gas_transport = model)
    ξ=(500.0, 0.5*P*x[1]*0.028/(R_GAS*T), 0.5*P*x[2]*0.044/(R_GAS*T))
    PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, cells),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = thermal, bottom = AdiabaticBC()),
            mass = (top = mass, bottom = ImpermeableBC()),
        ),
        T_initial = z->T,
        ξ_initial = [z->ξ[j] for j = 1:3],
        tspan = (0.0, 0.01),
    )
end
function exchange_snapshot(problem; energy = :temperature)
    cp=setup(problem; energy_state = energy)
    du=zero(cp.u0)
    residual!(du, cp.u0, cp.def, 0.0, cp.ws)
    ds, c=buffers(cp.ws, cp.def, Float64)
    return cp, ds, c, du
end

@testset "a trace reservoir retains its specified pressure" begin
    P=1e-5
    T=500.0
    x=(0.7, 0.3)
    reservoir=GasReservoir(P = P, T = T, mole_fractions = (a = x[1], b = x[2]))
    for energy in (:temperature, :enthalpy)
        cp, ds, c, du=exchange_snapshot(
            gas_problem(P = 0.0, T = T, mass = AmbientGasBC(reservoir));
            energy,
        )
        f=n_faces(cp.def.setup.topology)
        @test c.p_surface[2:3, f]≈P .* collect(x) rtol=2e-14
        # Wall-only half-cell transport to an initially empty pore. The
        # known reservoir pressure must not be smoothed with a cell tolerance.
        expected=-1e-5 .* P .* collect(x) .* [0.028, 0.044] ./ (R_GAS*T*(0.01/4))
        @test c.F_mass[2:3, f]≈expected rtol=2e-14
    end
end

@testset "trace pressure bias tracks the current pore volume and temperature" begin
    cp=setup(gas_problem(); pressure_tolerance = 1.0)
    layout=cp.def.layout
    u=copy(cp.u0)
    du=zero(u)
    tstar=(1+sqrt(33.0))/16
    for T in (250.0, 500.0, 1500.0),
        φ in (0.05, 0.5, 0.9),
        fraction in (0.0, tstar, 1.0, 2.0)

        V=cp.def.setup.reference.volume[1]
        for i = 1:2
            u[energy_index(layout, i)]=T
            u[mass_index(layout, i, 1)]=(1-φ)*1000V
            for (j, M) in ((2, 0.028), (3, 0.044))
                u[mass_index(layout, i, j)]=fraction*M*φ*V/(2R_GAS*T)
            end
        end
        residual!(du, u, cp.def, 0.0, cp.ws)
        ds, _=buffers(cp.ws, cp.def, Float64)
        @test all(p->-32eps()<=fraction-p<0.26, ds.P)
        fraction==0 && @test all(iszero, ds.P)
        fraction>=1 && @test ds.P≈fill(fraction, 2) rtol=3e-14
        @test all(isfinite, du)
    end
end

@testset "assembled reservoir equilibrium and incoming energy" begin
    for energy in (:temperature, :enthalpy)
        res=GasReservoir(P = 1e5, T = 300.0, mole_fractions = (a = 0.7, b = 0.3))
        cp, ds, c, du=exchange_snapshot(gas_problem(mass = AmbientGasBC(res)); energy)
        @test maximum(abs, c.F_mass)<1e-15
        @test ds.P≈fill(1e5, 2)
        @test maximum(abs, du)<1e-8
        cp, ds, c, du=exchange_snapshot(
            gas_problem(P = 5e4, mass = AmbientGasBC(res));
            energy,
        )
        f=n_faces(cp.def.setup.topology)
        @test all(c.F_mass[2:3, f] .< 0)
        h=Pyrolysis.States.component_enthalpies(cp.def.setup.laws, 300.0)
        @test c.F_energy[f]≈sum(c.F_mass[:, f] .* h) rtol=1e-12
        @test c.boundary_heat[f]==0
        @test c.surface_residual[f]<3e-12
    end
end

@testset "film and hydraulic surface balances" begin
    for model in (
        LumpedTransfer(),
        DustyGas(
            tortuosity = 2.0,
            D_ref = 2e-5,
            T_ref = 300.0,
            P_ref = 1e5,
            pore_diameter = 1e-6,
        ),
    )
        res=GasReservoir(P = 1e5, T = 300.0, mole_fractions = (a = 0.4, b = 0.6))
        for P in (8e4, 1e5, 1.2e5), g in (0.0, 1e-3, Inf)
            bc=AmbientGasBC(res; film = MassTransferFilm(g))
            cp, ds, c, du=exchange_snapshot(
                gas_problem(;
                    model,
                    P,
                    mass = bc,
                    thermal = ConvectiveBC(
                        h = 20.0,
                        T_inf = 300.0,
                        blowing = CouetteFilm(),
                    ),
                ),
            )
            f=n_faces(cp.def.setup.topology)
            @test c.surface_residual[f]<3e-12
            N=c.F_mass[2:3, f]
            ps=c.p_surface[2:3, f]
            Ys=ps .* [0.028, 0.044]
            Ys./=sum(Ys)
            Y∞=[0.4*0.028, 0.6*0.044]
            Y∞./=sum(Y∞)
            m=sum(N)
            if isinf(g)
                @test Ys≈Y∞ rtol=1e-10
            else
                expected=g==0 ? (m>=0 ? Ys : Y∞)*m :
                         Ys*m+(abs(m/g)<1e-8 ? g : m/expm1(m/g))*(Ys-Y∞)
                @test N≈expected rtol=1e-8 atol=1e-13
            end
            hs=Pyrolysis.States.component_enthalpies(cp.def.setup.laws, c.T_face[f])
            Hs=sum(c.F_mass[:, f] .* hs)
            qsurf=ds.k[end]*(c.T_face[f]-ds.T[end])/(0.01/4)
            @test c.F_energy[f]≈Hs-qsurf rtol=1e-10 atol=1e-9
        end
        for flow in (:darcy_velocity, :molar, :mass), P in (8e4, 1.2e5)
            R=flow===:darcy_velocity ? 1e6 : 1e4
            bc=ConvectivePressureBC(res; resistance = R, flow = flow)
            cp, ds, c, du=exchange_snapshot(gas_problem(; model, P, mass = bc, κ = 1e-12))
            f=n_faces(cp.def.setup.topology)
            rate=flow===:darcy_velocity ? c.darcy_velocity[f] :
                 flow===:mass ? sum(c.F_mass[:, f]) :
                 sum(c.F_mass[2:3, f] ./ [0.028, 0.044])
            @test sum(c.p_surface[:, f])-1e5≈R*rate rtol=1e-8 atol=1e-6
        end
    end
end

@testset "assembled wall-only / Darcy-only layers transmit total flux" begin
    a=SolidComponent(:wall; ρ = 500.0, φ = 0.5, c = 1000.0, k = 0.2, λ = 1e-5, κ = 0.0)
    b=SolidComponent(
        :hydraulic;
        ρ = 500.0,
        φ = 0.5,
        c = 1000.0,
        k = 0.2,
        λ = 0.0,
        κ = 1e-12,
    )
    gas=GaseousComponent(:gas; M = 0.03, c = 1000.0, k = 0.03)
    mat=Material(name = :switch_layers, components = (a, b, gas), gas_viscosity = 1e-5)
    L=0.01
    T=500.0
    PL=2e5
    PR=1e5
    d=L/4
    p=PyrolysisProblem(
        mesh = create_uniform_mesh(L, 2),
        material = mat,
        T_initial = z->T,
        ξ_initial = [
            z->z<L/2 ? 500.0 : 0.0,
            z->z<L/2 ? 0.0 : 500.0,
            z->0.5*0.03*(z<L/2 ? PL : PR)/(R_GAS*T),
        ],
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
        ),
        tspan = (0.0, 0.01),
    )
    cp, ds, c, du=exchange_snapshot(p)
    # Quadratic from equality of the two half-path fluxes, solved entirely
    # in high precision from the original physical coefficients.
    reference=setprecision(256) do
        a=BigFloat(1e-12)/(BigFloat(1e-5)*BigFloat(d))
        b=BigFloat(1e-5)/BigFloat(d)-a*BigFloat(PR)
        cc=-BigFloat(1e-5)*BigFloat(PL)/BigFloat(d)
        ps=(-b+sqrt(b*b-4a*cc))/(2a)
        Float64(
            BigFloat(0.03)*BigFloat(1e-5)*(BigFloat(PL)-ps)/(
                BigFloat(R_GAS)*BigFloat(T)*BigFloat(d)
            ),
        )
    end
    @test reference>0
    @test c.F_mass[3, 2]≈reference rtol=1e-13
    @test du[mass_index(cp.def.layout, 1, 3)]≈-du[mass_index(cp.def.layout, 2, 3)] rtol=1e-14
end
