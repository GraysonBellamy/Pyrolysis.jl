using Test, Pyrolysis
using Pyrolysis.Materials: R_GAS
using Pyrolysis.Properties: transport_matrix_state, transport_coefficients

function endpoint_problem(n, width; dusty = false)
    ep=dusty ?
       DustyEndpoint(
        pore_diameter = 2e-6,
        permeability = 0.0,
        tortuosity = 2.0,
        max_skeletal_occupancy = 0.25,
    ) : LumpedEndpoint(λ = 2e-5, κ = 0.0, max_skeletal_occupancy = 0.25)
    model=dusty ?
          DustyGas(
        tortuosity = 2.0,
        D_ref = 1e-5,
        T_ref = 300.0,
        P_ref = 1e5,
        pore_diameter = 1e-6,
        endpoint = ep,
        blend = (width/2, width),
    ) : LumpedTransfer(endpoint = ep, blend = (width/2, width))
    mat=Material(
        name = :OccupiedEndpoint,
        gas_transport = model,
        depletion_limiter = DepletionLimiter(enabled = false),
        components = (
            SolidComponent(
                :matrix;
                ρ = 500.0,
                φ = 0.5,
                c = 1000.0,
                k = 0.2,
                γ = 0.5,
                λ = 1e-5,
            ),
            SolidComponent(:residue; ρ = 500.0, φ = 0.5, c = 1000.0, k = 0.2, γ = 0.0),
            GaseousComponent(:gas; M = 0.03, c = 1000.0, k = 0.03),
        ),
        reactions = (
            Reaction(:decay, 1=>(2, 3); A = 2.0, E = 0.0, h = 0.0, yields = (0.2, 0.8)),
        ),
    )
    PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, n),
        material = mat,
        T_initial = z->500.0,
        ξ_initial = [z->500.0, z->0.0, z->0.0],
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (
                top = AmbientGasBC(GasReservoir(P = 0.0, T = 500.0)),
                bottom = ImpermeableBC(),
            ),
        ),
        tspan = (0.0, 6.0),
        output_times = collect(0.0:0.25:6.0),
    )
end
@testset "occupied endpoint on a reacting domain ($moving, $energy, $dusty)" for moving in
                                                                                 (
        false,
        true,
    ),
    energy in (:temperature, :enthalpy),
    dusty in (false, true)

    rec=ConservationRecorder()
    p=endpoint_problem(6, 0.02; dusty)
    sol=solve(
        p;
        use_ale = moving,
        energy_state = energy,
        callback = rec,
        reltol = 1e-8,
        abstol = 1e-10,
        pressure_tolerance = 0.001,
        dtmin = 1e-15,
    )
    @test completed(sol)
    @test minimum(get_porosity(sol))≥0.5-1e-10
    @test maximum(abs.(sol.T .- 500.0))<1e-5
    for i in axes(sol.T, 1)
        ξ=Tuple(sol.ξ[:, i, end])
        q, _=transport_matrix_state(ξ, sol.T[i, end], p.material)
        @test q<0.01
        @test ξ[2] / 1000.0 > 0.09  # surviving γ=0 skeleton occupies pore space
        @test transport_coefficients(
            ξ,
            sol.T[i, end],
            p.material,
            get_porosity(sol)[i, end],
        ).κ==0
    end
    audit=gas_species_audit(sol, rec).species.gas
    @test abs(audit.relative)<2e-5
    @test audit.production≈4*(1-exp(-12)) rtol=1e-6
end

@testset "endpoint blend and mesh sensitivity ($moving)" for moving in (false, true)
    escape=zeros(3, 3)
    for (i, n) in enumerate((6, 12, 24)), (j, width) in enumerate((0.02, 0.01, 0.005))
        rec=ConservationRecorder()
        sol=solve(
            endpoint_problem(n, width);
            use_ale = moving,
            callback = rec,
            reltol = 1e-9,
            abstol = 1e-11,
            pressure_tolerance = 0.0001,
            dtmin = 1e-15,
        )
        @test completed(sol)
        a=gas_species_audit(sol, rec).species.gas
        @test abs(a.relative)<1e-5
        escape[i, j]=a.escape
    end
    println((
        endpoint_meshes = (6, 12, 24),
        blend_upper = (0.02, 0.01, 0.005),
        moving,
        escape,
    ))
    flush(stdout)
    # Spatial refinement contracts the discretization error. The blend window
    # is a constitutive input: these three windows need not approach a common
    # model at a geometric rate. Qualify its sensitivity over the declared
    # range against a 2% escape budget, independently of spatial convergence.
    @test maximum(abs.(escape[3, :]-escape[2, :]))<maximum(abs.(escape[2, :]-escape[1, :]))/2
    @test maximum(escape[3, :])-minimum(escape[3, :])<0.02*minimum(escape[3, :])
end
