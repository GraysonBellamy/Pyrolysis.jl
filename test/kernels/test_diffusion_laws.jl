# Partial-pressure potential checks with independently differentiated physical
# profiles. Concentrations include pore volume; the superficial coefficient
# multiplies the intrinsic partial-pressure potential.

using Test
using Pyrolysis
using Pyrolysis.Materials: R_GAS
using Pyrolysis.Meshes: create_uniform_mesh
using Pyrolysis.States: mass_index
using Pyrolysis.Assembly: buffers, residual!
using Pyrolysis.Integration: setup

const L = 0.01
const D = 2e-5
const P0 = 101325.0
const M0 = 0.028

# Smooth imposed profiles; s is distance / L. The two gases have equal molar
# masses except in :molar_mass, so composition alone need not change ρ_g T.
function fields(case, z)
    s = z / L
    φ = case == :porosity ? 0.4 + 0.4s : 0.6
    P = case == :pressure ? P0 + 1e4s : P0
    T = case == :temperature ? 400.0 + 200s : 500.0
    Y = case in (:composition, :molar_mass) ? 0.2 + 0.4s : 0.3
    masses = case == :molar_mass ? (0.018, 0.044) : (M0, M0)
    M = inv(Y / masses[1] + (1 - Y) / masses[2])
    ρ = P * M / (R_GAS * T)
    return (; φ, P, T, Y, M, ρ, masses)
end

function snapshot(case, n; darcy, coefficient = D, energy = :temperature)
    κ = darcy ? 1e-12 : 0.0
    masses = fields(case, 0.0).masses
    mat = Material(
        name = :diffusion_laws,
        components = (
            SolidComponent(
                :dense,
                ρ = 1200.0,
                ρ_intrinsic = 2000.0,
                c = 1000.0,
                k = 0.2,
                λ = coefficient,
                κ = κ,
            ),
            SolidComponent(
                :porous,
                ρ = 400.0,
                ρ_intrinsic = 2000.0,
                c = 1000.0,
                k = 0.2,
                λ = coefficient,
                κ = κ,
            ),
            GaseousComponent(:a, M = masses[1], c = 1000.0, k = 0.03),
            GaseousComponent(:b, M = masses[2], c = 1000.0, k = 0.03),
        ),
        reactions = (),
        gas_viscosity = 1.8e-5,
    )
    # The solids fill unit bulk volume while their skeletal volumes give φ.
    concentrations(z) =
        let f = fields(case, z), θ = (f.φ - 0.4) / 0.4
            (1200.0 * (1 - θ), 400.0 * θ, f.φ * f.ρ * f.Y, f.φ * f.ρ * (1 - f.Y))
        end
    problem = PyrolysisProblem(
        mesh = create_uniform_mesh(L, n),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (top = ImpermeableBC(), bottom = ImpermeableBC()),
        ),
        T_initial = z -> fields(case, z).T,
        ξ_initial = [z -> concentrations(z)[j] for j = 1:4],
        tspan = (0.0, 1.0),
    )
    cp = setup(problem; energy_state = energy)
    du = zeros(length(cp.u0))
    residual!(du, cp.u0, cp.def, 0.0, cp.ws)
    ds, c = buffers(cp.ws, cp.def, Float64)
    faces = findall(f -> f.cells[2] != 0, cp.def.setup.topology.faces)
    return (; cp, du, ds, c, faces)
end

# Continuum references from x_j = Y_j M/M_j and n_j = -Λ ∂p_j/(RT).
# The separate continuum Darcy term uses the intrinsic donor density in the
# vanishing-mesh limit. A finite composite face is not a sum of independently
# averaged wall and hydraulic channels.
function reference_fluxes(case, z; κ = 0.0)
    f=fields(case, z)
    dY=case in (:composition, :molar_mass) ? 0.4/L : 0.0
    dP=case===:pressure ? 1e4/L : 0.0
    dM=-f.M^2*dY*(inv(f.masses[1])-inv(f.masses[2]))
    x=f.Y*f.M/f.masses[1]
    dx=(f.M*dY+f.Y*dM)/f.masses[1]
    dp=(f.P*dx+x*dP, -f.P*dx+(1-x)*dP)
    v=-κ/1.8e-5*dP
    return ntuple(j->f.masses[j]/(R_GAS*f.T)*(-D*dp[j]+v*f.P*(j==1 ? x : 1-x)), 2)
end

@testset "partial-pressure potentials through assembled kernel ($energy)" for energy in (
    :temperature,
    :enthalpy,
)
    for case in (:composition, :porosity, :pressure, :temperature)
        a=snapshot(case, 8; darcy = false, energy)
        b=snapshot(case, 8; darcy = true, energy)
        expected=hcat(
            [
                collect(reference_fluxes(case, a.ds.geometry.face_center[f][1])) for
                f in a.faces
            ]...,
        )
        @test a.c.F_mass[3:4, a.faces]≈expected atol=1e-12 rtol=1e-9
        for s in (a, b), i = 1:8
            f=fields(case, s.ds.geometry.cell_center[i][1])
            @test s.ds.P[i]≈f.P rtol=1e-10
            @test s.ds.φ[i]≈f.φ rtol=1e-10
        end
        for s in (a, b), j in (3, 4)
            rows=[s.du[mass_index(s.cp.def.layout, i, j)] for i = 1:8]
            @test abs(sum(rows))<=1e-12*max(sum(abs, rows), 1.0)
        end
        if case in (:porosity, :temperature)
            @test maximum(abs, a.c.F_mass[3:4, a.faces])<1e-12
            @test maximum(abs, b.c.F_mass[3:4, b.faces])<1e-12
        elseif case===:composition
            @test b.c.F_mass[3:4, b.faces]≈expected atol=1e-12 rtol=1e-9
        end
    end
end

@testset "unequal masses and total coupled flow converge to continuum" begin
    for (case, darcy, order) in ((:molar_mass, false, 1.8), (:pressure, true, 0.8))
        errors=Float64[]
        for n in (8, 16, 32)
            s=snapshot(case, n; darcy)
            expected=hcat(
                [
                    collect(
                        reference_fluxes(
                            case,
                            s.ds.geometry.face_center[f][1];
                            κ = darcy ? 1e-12 : 0.0,
                        ),
                    ) for f in s.faces
                ]...,
            )
            numerical=s.c.F_mass[3:4, s.faces]
            push!(errors, maximum(abs, numerical-expected))
            if case===:molar_mass
                # Equal wall mobility yields equimolar counterdiffusion, not
                # zero total mass flux when molar masses differ.
                @test maximum(abs, vec(sum(numerical ./ [0.018, 0.044]; dims = 1)))<1e-11
                @test minimum(vec(sum(numerical; dims = 1)))>1e-5
            end
        end
        @test all(>(order), log2.(errors[1:2] ./ errors[2:3]))
    end
end

# A long production transient at a discontinuous pore state. The two regions
# have different intrinsic density, wall coefficient and Darcy mobility, but
# share each species' partial pressure. No equilibrium diffusion is permitted.
@testset "sealed heterogeneous equilibrium ($energy, $dusty)" for energy in
                                                                  (:temperature, :enthalpy),
    dusty in (false, true)

    mat=Material(
        name = :sealed_jump,
        gas_transport = dusty ?
                        DustyGas(
            tortuosity = 2.0,
            D_ref = 2e-5,
            T_ref = 300.0,
            P_ref = 1e5,
        ) : LumpedTransfer(),
        components = (
            SolidComponent(
                :dense;
                ρ = 1200.0,
                ρ_intrinsic = 2000.0,
                c = 1000.0,
                k = 0.2,
                λ = 1e-5,
                κ = 1e-14,
            ),
            SolidComponent(
                :porous;
                ρ = 400.0,
                ρ_intrinsic = 2000.0,
                c = 1000.0,
                k = 0.2,
                λ = 7e-5,
                κ = 4e-13,
            ),
            GaseousComponent(:light; M = 0.004, c = 1000.0, k = 0.03),
            GaseousComponent(:heavy; M = 0.044, c = 1000.0, k = 0.03),
        ),
    )
    φ=z->z<0.005 ? 0.4 : 0.8
    p=PyrolysisProblem(
        mesh = create_stretched_mesh(0.01, 6, 1.2),
        material = mat,
        T_initial = z->500.0,
        ξ_initial = [
            z->z<0.005 ? 1200.0 : 0.0,
            z->z<0.005 ? 0.0 : 400.0,
            z->φ(z)*7e4*0.004/(R_GAS*500),
            z->φ(z)*3e4*0.044/(R_GAS*500),
        ],
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (top = ImpermeableBC(), bottom = ImpermeableBC()),
        ),
        tspan = (0.0, 1e4),
        output_times = [0.0, 1.0, 100.0, 1e4],
    )
    sol=solve(
        p;
        energy_state = energy,
        reltol = 1e-9,
        abstol = 1e-11,
        pressure_tolerance = 0.0001,
    )
    @test completed(sol)
    @test maximum(abs.(get_partial_pressures(sol)[1, :, :] .- 7e4))<1e-7
    @test maximum(abs.(get_partial_pressures(sol)[2, :, :] .- 3e4))<1e-7
    @test maximum(abs.(sol.T .- 500.0))<1e-9
    @test sol.m[:, :, end]≈sol.m[:, :, 1] rtol=1e-13
end
