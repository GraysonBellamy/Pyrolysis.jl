using Test, Pyrolysis, StaticArrays, ForwardDiff
using Pyrolysis.Materials: R_GAS
using Pyrolysis.Assembly: residual!, buffers
using Pyrolysis.States: mass_index, energy_index
using Pyrolysis.Meshes: n_faces
using Pyrolysis.Kernels: find_kernel
using Pyrolysis.BoundaryConditions: molar_film_flux
using Pyrolysis.ColoredJacobian:
    dense_jacobian, missing_entries, compare_entries, jac_prototype, jacobian!

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
        mesh = uniform_mesh(0.01, cells),
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
        for flow in (:darcy_velocity, :mass), P in (8e4, 1.2e5)
            R=flow===:darcy_velocity ? 1e6 : 1e4
            bc=ConvectivePressureBC(res; resistance = R, flow = flow)
            cp, ds, c, du=exchange_snapshot(gas_problem(; model, P, mass = bc, κ = 1e-12))
            f=n_faces(cp.def.setup.topology)
            rate=flow===:darcy_velocity ? c.darcy_velocity[f] : sum(c.F_mass[:, f])
            @test sum(c.p_surface[:, f])-1e5≈R*rate rtol=1e-8 atol=1e-6
        end
        # The resistance may be a law of the surface temperature, like a conductance.
        for P in (8e4, 1.2e5)
            Rlaw=T->1e4*(T/300.0)^1.7
            bc=ConvectivePressureBC(res; resistance = Rlaw, flow = :mass)
            cp, ds, c, du=exchange_snapshot(
                gas_problem(;
                    model,
                    P,
                    T = 600.0,
                    mass = bc,
                    κ = 1e-12,
                    thermal = ConvectiveBC(h = 20.0, T_inf = 300.0),
                ),
            )
            f=n_faces(cp.def.setup.topology)
            @test c.T_face[f]!=300.0
            @test sum(c.p_surface[:, f])-1e5≈Rlaw(c.T_face[f])*sum(c.F_mass[:, f]) rtol=1e-8 atol=1e-6
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
        mesh = uniform_mesh(L, 2),
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

@testset "the analogy film is the twin of the convective film" begin
    res=GasReservoir(P = 1e5, T = 300.0, mole_fractions = (a = 0.4, b = 0.6))
    Y∞=[0.4*0.028, 0.6*0.044]
    Y∞./=sum(Y∞)
    h=20.0
    # Constant heat capacities, so c_p,∞ does not depend on the film temperature.
    cp∞=Y∞[1]*1000.0+Y∞[2]*1200.0
    for P in (8e4, 1.2e5), Le in (1.0, 2.5), combined in (false, true)
        conv=ConvectiveBC(h = h, T_inf = 300.0, blowing = CouetteFilm())
        thermal=combined ? HeatFluxBC(1e4)+conv : conv
        g=h/(cp∞*Le^(2/3))
        cpA, _, cA, _=exchange_snapshot(
            gas_problem(;
                P,
                mass = AmbientGasBC(res; film = AnalogyFilm(Le = Le)),
                thermal,
            ),
        )
        cpE, _, cE, _=exchange_snapshot(
            gas_problem(; P, mass = AmbientGasBC(res; film = MassTransferFilm(g)), thermal),
        )
        f=n_faces(cpA.def.setup.topology)
        @test cA.F_mass[:, f]≈cE.F_mass[:, f] rtol=1e-12
        @test cA.p_surface[:, f]≈cE.p_surface[:, f] rtol=1e-12
        @test cA.F_energy[f]≈cE.F_energy[f] rtol=1e-12
        @test cA.T_face[f]≈cE.T_face[f] rtol=1e-12
        @test cA.surface_residual[f]<3e-12
    end
    # The default of AmbientGasBC is the analogy film.
    @test AmbientGasBC(res).film isa AnalogyFilm
    cpD, _, cD, _=exchange_snapshot(
        gas_problem(;
            P = 8e4,
            mass = AmbientGasBC(res),
            thermal = ConvectiveBC(h = h, T_inf = 300.0),
        ),
    )
    cpA, _, cA, _=exchange_snapshot(
        gas_problem(;
            P = 8e4,
            mass = AmbientGasBC(res; film = AnalogyFilm()),
            thermal = ConvectiveBC(h = h, T_inf = 300.0),
        ),
    )
    f=n_faces(cpD.def.setup.topology)
    @test cD.F_mass[:, f]==cA.F_mass[:, f]
    # Without a convective condition, or against a vacuum, the twin film is
    # absent and the reservoir is in direct contact, exactly as `film = nothing`.
    for (P, Pr, thermal) in (
        (8e4, 1e5, AdiabaticBC()),
        (1.2e5, 1e5, HeatFluxBC(1e4)),
        (8e4, 0.0, ConvectiveBC(h = h, T_inf = 300.0)),
    )
        r=GasReservoir(P = Pr, T = 300.0, mole_fractions = (a = 0.4, b = 0.6))
        cpA, _, cA, _=exchange_snapshot(
            gas_problem(; P, mass = AmbientGasBC(r; film = AnalogyFilm()), thermal),
        )
        cpN, _, cN, _=exchange_snapshot(
            gas_problem(; P, mass = AmbientGasBC(r; film = nothing), thermal),
        )
        f=n_faces(cpA.def.setup.topology)
        @test cA.F_mass[:, f]==cN.F_mass[:, f]
        @test cA.p_surface[:, f]==cN.p_surface[:, f]
        @test cA.F_energy[f]==cN.F_energy[f]
        @test cA.T_face[f]==cN.T_face[f]
    end
end

# A purged face carrying the pooling model: matrix, one volatile and nitrogen,
# the volatile at a pore-transport molar mass `Mpore` the model need not share.
function pool_problem(;
    volatiles = GasSpecies(:pyrolysate; M = 0.042, c = 2000.0, k = 0.05, V_D = 65.4),
    Mpore = 0.070,
    P = 1.002e5,
    xv = 0.9,
    T = 870.0,
    κ = 2e-12,
    film = AnalogyFilm(),
    blowing = PoolingFilm(NITROGEN; volatiles),
    thermal = HeatFluxBC(5e3)+ConvectiveBC(h = 11.2, T_inf = 360.0, blowing = blowing),
    reservoir = (N2 = 1.0,),
    cells = 2,
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
    vol=GaseousComponent(:volatiles; M = Mpore, c = 2000.0, k = 0.05)
    n2=GaseousComponent(:N2; M = 0.028014, c = 1100.0, k = 0.05)
    mat=Material(name = :pooled_face, components = (solid, vol, n2))
    ξ=(500.0, 0.5*P*xv*Mpore/(R_GAS*T), 0.5*P*(1-xv)*0.028014/(R_GAS*T))
    PyrolysisProblem(
        mesh = uniform_mesh(0.01, cells),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = thermal, bottom = AdiabaticBC()),
            mass = (
                top = AmbientGasBC(
                    GasReservoir(P = 1e5, T = 360.0, mole_fractions = reservoir);
                    film = film,
                ),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z->T,
        ξ_initial = [z->ξ[j] for j = 1:3],
        tspan = (0.0, 0.01),
    )
end

@testset "a pooling face exchanges gas through the model's own molar film" begin
    BC=Pyrolysis.BoundaryConditions
    surrogate=GasSpecies(:pyrolysate; M = 0.042, c = 2000.0, k = 0.05, V_D = 65.4)
    mma=GasSpecies(:MMA; M = 0.100, c = 2000.0, k = 0.03, V_D = 110.2)
    ψs=Float64[]
    ingress=false
    for (volatiles, Mpore) in ((surrogate, 0.070), (mma, 0.100)),
        P in (1.0005e5, 1.002e5, 1.005e5, 1.01e5)

        model=PoolingFilm(NITROGEN; volatiles)
        cp, ds, c, du=exchange_snapshot(
            pool_problem(; volatiles, Mpore, P, blowing = model),
        )
        f=n_faces(cp.def.setup.topology)
        mat=cp.def.setup.material
        @test c.surface_residual[f]<3e-12
        @test all(isfinite, du)
        N=c.F_mass[2:3, f]
        @test N[1]>0
        ingress|=N[2]<0                 # nitrogen diffusing in against the outflow
        # The solved surface state, in the material's basis.
        w=c.p_surface[2:3, f] .* [Mpore, 0.028014]
        Ys=w/sum(w)
        Y=(0.0, Ys[1], Ys[2])
        Ts=c.T_face[f]
        # The film the model states, with moles counted at its molar masses,
        # carries exactly the pore fluxes: mass and moles balance.
        Mf=collect(BC.film_molar_masses(model, mat)[2:3])
        @test Mf==[volatiles.molar_mass, NITROGEN.molar_mass]
        n=Ys ./ Mf
        g=BC.pooling_film_conductance(model, Y, mat, Ts, 360.0)
        F=molar_film_flux(n/sum(n), [0.0, 1.0], sum(N ./ Mf), g/NITROGEN.molar_mass, Mf)
        @test F≈N rtol=1e-9 atol=1e-15
        @test sum(F)≈sum(N) rtol=1e-9
        @test sum(F ./ Mf)≈sum(N ./ Mf) rtol=1e-9
        # Heat and gas exchange share one ψ.
        p=BC.pooling_factor(model, Y, mat, Ts, 360.0)
        st=find_kernel(Pyrolysis.Kernels.ThermalBoundaryKernel, cp.def.kernels).states[f]
        @test st.pooling_factor≈p.ψ rtol=1e-12
        @test st.face_fraction≈p.Y rtol=1e-12
        volatiles===mma && push!(ψs, p.ψ)
    end
    # The heavy gas crosses the switch over the sweep and the Newton holds,
    # with nitrogen entering against the outflow at the low end.
    @test minimum(ψs)<0.8<0.95<maximum(ψs)
    @test ingress

    # Nitrogen alone at the reservoir state: nothing crosses the face, and ψ
    # is the switch's value at Λ = 1.
    cp, ds, c, du=exchange_snapshot(pool_problem(; P = 1e5, xv = 0.0, T = 360.0))
    f=n_faces(cp.def.setup.topology)
    @test maximum(abs, c.F_mass[:, f])<1e-14
    st=find_kernel(Pyrolysis.Kernels.ThermalBoundaryKernel, cp.def.kernels).states[f]
    @test c.T_face[f]>360.0
    @test st.pooling_factor≈0.999003997 atol=1e-9

    # The face's own film is the default; everything that would leave its gas
    # exchange without the model's ψ, or ambiguous, is refused.
    conv=ConvectiveBC(
        h = 11.2,
        T_inf = 360.0,
        blowing = PoolingFilm(NITROGEN; volatiles = surrogate),
    )
    @test_throws ArgumentError exchange_snapshot(
        pool_problem(; film = MassTransferFilm(9.41e-3)),
    )
    @test_throws ArgumentError exchange_snapshot(
        pool_problem(; film = AnalogyFilm(Le = 1.3)),
    )
    @test_throws ArgumentError exchange_snapshot(
        pool_problem(; thermal = conv+ConvectiveBC(h = 2.0, T_inf = 360.0)),
    )
    @test_throws ArgumentError exchange_snapshot(
        pool_problem(; reservoir = (volatiles = 0.1, N2 = 0.9)),
    )

    # The resolved film is inferred through combined conditions.
    K=Pyrolysis.Kernels
    th=HeatFluxBC(1e4)+conv+RadiativeBC(T_inf = 360.0)
    @test @inferred(K._resolve_film(AnalogyFilm(), th)) isa K.PoolingGasFilm
    @test @inferred(
        K._resolve_film(
            AnalogyFilm(),
            HeatFluxBC(1e4)+ConvectiveBC(h = 11.2, T_inf = 360.0),
        )
    ) isa AnalogyFilm
    @test @inferred(K._resolve_film(AnalogyFilm(), HeatFluxBC(1e4))) === nothing
end

@testset "the pooling film's Jacobian" begin
    for P in (1.002e5, 1.01e5), energy in (:temperature, :enthalpy)
        mma=GasSpecies(:MMA; M = 0.100, c = 2000.0, k = 0.03, V_D = 110.2)
        cp=setup(
            pool_problem(; volatiles = mma, Mpore = 0.100, P, cells = 3);
            energy_state = energy,
        )
        def, ws, cj=cp.def, cp.ws, cp.jacobian
        u=copy(cp.u0)
        Jd=dense_jacobian(def, ws, u, 0.0)
        @test all(isfinite, Jd)
        @test isempty(missing_entries(Jd, cj.pattern))
        J=jac_prototype(cj)
        jacobian!(J, cj, u, 0.0)
        @test compare_entries(Matrix(J), Jd; rtol = 1e-12).passed
        # Automatic derivatives through the implicit surface root against
        # central differences of the residual.
        du⁺=zero(u)
        du⁻=zero(u)
        Jf=similar(Jd)
        for j in eachindex(u)
            h=1e-6*max(abs(u[j]), 1e-6)
            v=copy(u)
            v[j]+=h
            residual!(du⁺, v, def, 0.0, ws)
            v[j]-=2h
            residual!(du⁻, v, def, 0.0, ws)
            Jf[:, j]=(du⁺-du⁻)/(2h)
        end
        @test compare_entries(Jf, Jd; rtol = 1e-5).passed
    end
end
