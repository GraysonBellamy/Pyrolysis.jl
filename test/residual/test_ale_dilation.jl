"""
ALE mesh-motion bookkeeping of the species equations.

Every component — condensed and gaseous — rides the mesh: the mesh velocity
is the skeleton velocity and the Darcy–Fick gas flux is defined relative to
the skeleton. The residual therefore carries no mesh-motion advection for
any species; its only mesh-motion term is the dilation `−ξⱼθᵢ` with the
cell's own strain rate, so that `d(ξⱼV)/dt = ω̇ⱼV − [JⱼA]_faces` for every
component. Pinned here on a sealed, uniformly collapsing column where the
flux term vanishes and the balance is exact cell by cell.
"""

using Test
using Pyrolysis
using Pyrolysis.Residual: residual!
using Pyrolysis.Problem: z_offset

@testset "ALE dilation applies to every component" begin
    ρ_A = 1000.0
    T_iso = 400.0
    A_r = 2.0
    L = 2e-3
    n_cells = 8

    # A → 0.3 C + 0.7 G at a uniform, isothermal rate (E = 0, h = 0), sealed
    # column: no flux crosses any face, so for each cell and component
    #   dξⱼ/dt = ω̇ⱼ − ξⱼθ,   θ = Σₖ γₖ ω̇ₖ/ρₖ = (−r + 0.3 r)/ρ = −0.7 r/ρ.
    solidA = SolidComponent(:A, ρ = ρ_A, c = 1000.0, k = 0.5)
    char = SolidComponent(:C, ρ = ρ_A, c = 1000.0, k = 0.5)
    gas = GaseousComponent(:G, M = 0.05, c = 1000.0, k = 0.02, λ = 1e-8)
    rxn = Reaction(:r, 1 => (2, 3); A = A_r, E = 0.0, h = 0.0, yields = (0.3, 0.7))
    mat = Material(
        name = :dilation_test,
        components = (solidA, char, gas),
        reactions = (rxn,),
        depletion_limiter = DepletionLimiter(enabled = false),
    )
    ξ0 = (0.6 * ρ_A, 0.1 * ρ_A, 3.0)   # some char and held gas already present
    mesh = create_uniform_mesh(L, n_cells, Val(3); T_initial = T_iso, ξ_initial = ξ0)
    bc_set = BoundaryConditionSet(
        thermal = (top = DirichletThermalBC(T_iso), bottom = DirichletThermalBC(T_iso)),
    )
    problem = PyrolysisProblem(
        mesh = mesh,
        material = mat,
        bc_set = bc_set,
        T_initial = z -> T_iso,
        ξ_initial = [z -> ξ0[1], z -> ξ0[2], z -> ξ0[3]],
        tspan = (0.0, 1.0),
    )
    def, ws, u = initialize_state(problem; use_ale = true)
    du = similar(u)
    residual!(du, u, def, 0.0, ws)

    n = ws.layout.n_cells
    r = A_r * ξ0[1]
    θ = -0.7 * r / ρ_A
    ω̇ = (-r, 0.3 * r, 0.7 * r)
    @test θ < 0

    for i = 1:n, j = 1:3
        expected = ω̇[j] - ξ0[j] * θ
        @test du[n+(j-1)*n+i] ≈ expected rtol = 1e-12
    end

    # The held gas is compressed uniformly: no cell is singled out and no
    # gas crosses a face. Equivalent statement: d(ξⱼV)/dt = ω̇ⱼV per cell.
    zoff = z_offset(ws.layout)
    for i = 1:n
        V = ws.state_cache.volume[i]
        dV = (du[zoff+i+1] - du[zoff+i]) * ws.state_cache.A_section[]
        @test dV ≈ θ * V rtol = 1e-12
        for j = 1:3
            @test du[n+(j-1)*n+i] * V + ξ0[j] * dV ≈ ω̇[j] * V rtol = 1e-12
        end
    end
end
