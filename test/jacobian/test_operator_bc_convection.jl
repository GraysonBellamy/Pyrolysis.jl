"""
`BCConvectionOperator` correctness tests (deferred audit Minor 5).

The operator is the single owner of the boundary cell's `S_conv` energy
row — the product `−Σ_g cp_g·J̄_g·(∂T/∂x)_b` whose blended gradient runs
through the surface-Newton `T_s` and whose face-average `J̄` includes both
the BC mass flux and the mode-appropriate interior-face flux.

Three gates:

1. **`apply!` equivalence (Fickian + Darcy).** The operator's boundary-cell
   `dT` contribution matches the residual's full boundary-cell `S_conv/ρcp`
   slice, reconstructed by running the residual's own kernels in
   `compute_rhs!` order: gas fluxes → `apply_mass_bc_flux!` →
   `apply_thermal_bc_flux!` (fills the boundary `q`) →
   `compute_gas_convective_source!` (which blends the boundary-face slope
   `−q/k` at the solved `T_s` into the boundary cells' gradient).

2. **`contribute_to_jacobian!` equivalence.** Per-boundary AD through the
   surface Newton + both face-flux kernels agrees with `ForwardDiff.jacobian`
   over the same reference closure on the declared support (T-row[bdy] ×
   `(T, ξ)` of the boundary cell and its neighbor); off-support reference
   entries are zero.

3. **Load-bearing chain.** The `T_s` chain is non-trivial in the fixture
   (nonlinear surface balance: convective + reradiative on top of the
   prescribed flux), so agreement is not vacuous.
"""

using Test
using SparseArrays
using ForwardDiff
using Pyrolysis
using Pyrolysis: to_problem_def, build_workspace
using Pyrolysis.Materials: n_components, Material
using Pyrolysis.Properties: PropertyCache, update_properties!
using Pyrolysis.Discretization:
    apply_thermal_bc_flux!,
    apply_mass_bc_flux!,
    compute_conductive_fluxes!,
    update_darcy_properties!
using Pyrolysis.Physics:
    compute_gas_fluxes!, compute_gas_fluxes_darcy_fick!, compute_gas_convective_source!
using Pyrolysis.Problem: populate!, build_ad_workspace
using Pyrolysis.Jacobian:
    BCConvectionOperator,
    build_bc_convection_workspace,
    StructuredJacobian,
    build_residual_context

include(joinpath(@__DIR__, "..", "perf", "perf_helpers.jl"))

function _build_bcc_problem(; darcy::Bool)
    virgin = if darcy
        SolidComponent(:wood, ρ = 500.0, c = 1500.0, k = 0.15, ε = 0.9, κ = 1e-14, λ = 1e-5)
    else
        SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86, λ = 1e-5)
    end
    gas = GaseousComponent(:gas, M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)
    rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
    mat = Material(name = :bcc_fixture, components = (virgin, gas), reactions = (rxn,))
    mesh = create_uniform_mesh(0.01, 8, Val(2))
    # Nonlinear surface balance (convective + reradiative on top of the
    # prescribed flux) so the T_s Newton chain is non-trivial; permeable top
    # so the boundary-face mass flux is in play.
    bc = BoundaryConditionSet(
        thermal = (
            top = HeatFluxBC(50e3) +
                  ConvectiveBC(h = 25.0, T_inf = 300.0) +
                  RadiativeBC(T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
        mass = (top = DiffusiveMassBC(), bottom = ImpermeableBC()),
    )
    problem =
        PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
    Pyrolysis.Problem.apply_initial_conditions!(problem.mesh, problem)
    def = to_problem_def(problem; use_ale = false)
    ws = build_workspace(problem, def)
    u = create_state_vector(ws.layout, ws.mesh)
    n = ws.layout.n_cells
    # Steep ramp + gas everywhere: the regime where the boundary S_conv row
    # is load-bearing.
    for i = 1:n
        u[i] = 600.0 + 400.0 * (i - 1) / (n - 1)
        u[n+i] = darcy ? 400.0 : 1100.0
        u[2n+i] = 8.0
    end
    return def, ws, u
end

# Reference: the residual's boundary-cell S_conv/ρcp slice, by running the
# residual's own kernels in compute_rhs! order. Generic over Dual eltypes.
function _reference_bcc_contribution(
    state,
    mesh,
    material::Material{NC},
    prop_cache,
    J_gas_buf,
    q_buf,
    S_conv_buf,
    bc_set,
    top_face_ids,
    bottom_face_ids,
    t;
    darcy::Bool,
) where {NC}

    update_properties!(prop_cache, state, mesh, material)
    if darcy
        update_darcy_properties!(
            prop_cache,
            state,
            mesh,
            material,
            bc_set,
            top_face_ids,
            bottom_face_ids,
        )
        compute_gas_fluxes_darcy_fick!(
            J_gas_buf,
            state,
            mesh,
            prop_cache,
            material,
            prop_cache.darcy.v_gas,
            prop_cache.darcy.φ,
        )
    else
        compute_gas_fluxes!(J_gas_buf, state, mesh, prop_cache, material)
    end
    apply_mass_bc_flux!(
        J_gas_buf,
        state,
        mesh,
        prop_cache,
        bc_set,
        material,
        t,
        top_face_ids,
        bottom_face_ids,
    )
    compute_conductive_fluxes!(q_buf, state, mesh, prop_cache)
    apply_thermal_bc_flux!(
        q_buf,
        state,
        mesh,
        bc_set,
        prop_cache,
        material,
        J_gas_buf,
        t,
        top_face_ids,
        bottom_face_ids,
        Pyrolysis.Materials.NO_RADIATION;
        use_transpiration_bc = false,
    )
    compute_gas_convective_source!(
        S_conv_buf,
        J_gas_buf,
        q_buf,
        prop_cache.k,
        state,
        mesh,
        material,
    )

    n_cells = length(state.T)
    Tu = eltype(state.T)
    dT = zeros(Tu, n_cells)
    dT[1] = S_conv_buf[1] / prop_cache.ρcp[1]
    dT[n_cells] = S_conv_buf[n_cells] / prop_cache.ρcp[n_cells]
    return dT
end

function _reference_bcc_jacobian(u, ws, def, ::Material{NC}, bc_set; darcy::Bool) where {NC}
    layout = ws.layout
    mesh = ws.mesh
    material = def.material
    n_cells = layout.n_cells
    n_total = length(u)
    x0 = copy(u)
    f_buf = zeros(n_cells)
    cfg = ForwardDiff.JacobianConfig(
        nothing,
        f_buf,
        x0,
        ForwardDiff.Chunk{min(12, n_total)}(),
    )
    Tu = eltype(cfg.duals[1])
    ws_ad = build_ad_workspace(def, ws, Tu)

    J_gas_buf = zeros(Tu, NC, length(mesh.faces))
    q_buf = zeros(Tu, length(mesh.faces))
    S_conv_buf = zeros(Tu, n_cells)

    f! =
        let ws_ad = ws_ad,
            layout = layout,
            mesh = mesh,
            material = material,
            def = def,
            bc_set = bc_set,
            top_face_ids = ws.top_face_ids,
            bottom_face_ids = ws.bottom_face_ids,
            darcy = darcy

            (out, x) -> begin
                populate!(ws_ad.state_cache, x, layout)
                Pyrolysis.Problem.update_state_geometry!(
                    ws_ad.state_cache,
                    mesh,
                    Pyrolysis.Problem.mesh_mode(def.mode),
                )
                dT = _reference_bcc_contribution(
                    ws_ad.state_cache,
                    mesh,
                    material,
                    ws_ad.prop_cache,
                    J_gas_buf,
                    q_buf,
                    S_conv_buf,
                    bc_set,
                    top_face_ids,
                    bottom_face_ids,
                    0.0;
                    darcy = darcy,
                )
                @inbounds for i = 1:n_cells
                    out[i] = dT[i]
                end
                return nothing
            end
        end
    return ForwardDiff.jacobian(f!, f_buf, x0, cfg)
end

@testset "BCConvectionOperator" begin

    @testset "interface declarations" begin
        op = BCConvectionOperator()
        @test Pyrolysis.Jacobian.reads(op) == (:T_s,)
        @test Pyrolysis.Jacobian.writes(op) == ()
        @test Pyrolysis.Jacobian.stencil(op) == (cells = -1:1,)
    end

    @testset "apply! matches residual boundary S_conv slice ($(mode))" for (mode, darcy) in
                                                                           (
        ("Fickian", false),
        ("Darcy", true),
    )
        def, ws, u = _build_bcc_problem(; darcy)
        n = ws.layout.n_cells
        NC = n_components(ws.layout)
        material = def.material
        bc_set = def.bcs

        populate!(ws.state_cache, u, ws.layout)
        Pyrolysis.Problem.update_state_geometry!(
            ws.state_cache,
            ws.mesh,
            Pyrolysis.Problem.mesh_mode(def.mode),
        )

        prop_cache_ref = PropertyCache(n, length(ws.mesh.faces))
        J_gas_ref = zeros(NC, length(ws.mesh.faces))
        q_ref = zeros(length(ws.mesh.faces))
        S_conv_ref = zeros(n)
        dT_ref = _reference_bcc_contribution(
            ws.state_cache,
            ws.mesh,
            material,
            prop_cache_ref,
            J_gas_ref,
            q_ref,
            S_conv_ref,
            bc_set,
            ws.top_face_ids,
            ws.bottom_face_ids,
            0.0;
            darcy = darcy,
        )

        dT_op = zeros(n)
        Pyrolysis.Jacobian.apply!(
            dT_op,
            BCConvectionOperator(),
            ws.state_cache,
            ws.mesh,
            material,
            bc_set,
            0.0,
            darcy,
            Pyrolysis.Materials.NO_RADIATION,
        )

        # The row is load-bearing in this fixture (escaping flux against a
        # steep boundary-layer gradient).
        @test abs(dT_ref[n]) > 1e-3
        @test isapprox(dT_op, dT_ref; rtol = 1e-8, atol = 1e-10)
    end

    @testset "contribute_to_jacobian! matches DenseAD on declared support" begin
        def, ws, u = _build_bcc_problem(; darcy = false)
        n = ws.layout.n_cells
        NC = n_components(ws.layout)
        material = def.material
        bc_set = def.bcs

        populate!(ws.state_cache, u, ws.layout)
        Pyrolysis.Problem.update_state_geometry!(
            ws.state_cache,
            ws.mesh,
            Pyrolysis.Problem.mesh_mode(def.mode),
        )

        # Declared support: T-row[bdy] × (T, ξ) of the boundary cell and its
        # interior neighbor.
        n_total = length(u)
        I = Int[]
        J_idx = Int[]
        V = Float64[]
        for (cell_id, nb_id) in ((1, 2), (n, n - 1))
            for col_cell in (cell_id, nb_id)
                push!(I, cell_id)
                push!(J_idx, col_cell)
                push!(V, 1.0)
                for j = 1:NC
                    ξ_off = n + (j - 1) * n
                    push!(I, cell_id)
                    push!(J_idx, ξ_off + col_cell)
                    push!(V, 1.0)
                end
            end
        end
        A = sparse(I, J_idx, V, n_total, n_total)
        fill!(A.nzval, 0.0)
        jac = StructuredJacobian(:eulerian, A)
        ctx = build_residual_context([BCConvectionOperator()])
        ws_op = build_bc_convection_workspace(material)

        Pyrolysis.Jacobian.contribute_to_jacobian!(
            jac,
            BCConvectionOperator(),
            u,
            ctx,
            0.0,
            ws_op,
            ws.state_cache,
            ws.mesh,
            material,
            ws.layout,
            bc_set,
            false,
            Pyrolysis.Materials.NO_RADIATION,
        )

        J_ref = _reference_bcc_jacobian(u, ws, def, material, bc_set; darcy = false)

        J_op = Matrix(jac.A)
        for row in (1, n)
            row_scale = max(maximum(abs, J_ref[row, :]), 1.0)
            for (r, c) in zip(I, J_idx)
                r == row || continue
                @test abs(J_op[r, c] - J_ref[r, c]) <= 1e-6 * row_scale
            end
        end

        # Sanity: the T_s chain is non-trivial — the top T-row's diagonal and
        # neighbor entries are live.
        @test abs(J_ref[n, n]) > 1e-3
        @test abs(J_ref[n, n-1]) > 1e-3

        # Support-completeness: the reference boundary rows reach nothing
        # outside the declared (T, ξ) columns of the boundary cell + neighbor.
        n_phys = n * (1 + NC)
        Pmask = falses(2, n_phys)
        for (r, c) in zip(I, J_idx)
            ri = r == 1 ? 1 : 2
            c <= n_phys && (Pmask[ri, c] = true)
        end
        @test maximum(abs.(J_ref[1, 1:n_phys][.!Pmask[1, :]])) < 1e-9
        @test maximum(abs.(J_ref[n, 1:n_phys][.!Pmask[2, :]])) < 1e-9
    end
end
