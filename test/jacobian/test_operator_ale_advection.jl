"""
`ALEAdvectionOperator` correctness tests.

The operator is **gas-species-only** and mirrors the conservative flux-form
kernel `compute_ale_species_advection!`: production has no energy
mesh-motion advection (solids are mesh-bound; T is advected only by the
mesh−material relative velocity, which is identically zero — the
`w_override` verification hook is the sole exception and is out of
structured-Jacobian scope). Every gate below therefore expects zero T-row
contributions.

Three gates:

1. **`apply!` equivalence.** Calling
   `apply!(dT, dξ, ALEAdvectionOperator(), state, mesh, w_node, material)`
   reproduces the residual's gas flux-form transport bit-for-bit and leaves
   `dT` untouched. The operator wraps the shared
   [`compute_ale_species_advection!`](@ref) function — there is no math
   duplication.

2. **`contribute_to_jacobian!` direct chain rule.** The flux-form term is
   linear in ξ at fixed `w_node`, so the operator writes analytic entries.
   We verify against ForwardDiff over a closure that holds `w_node` fixed,
   so the AD reference and the operator's writes both ignore the indirect
   (cumulative-`w`) chain rule.

3. **`contribute_to_jacobian!` indirect chain rule.** Once we densify the
   pushed [`EmbeddedScaledPrefixSumCoupling`](@ref) entries onto the
   sparse A's pattern (the [`DirectSolve`](@ref) strategy's combined view),
   the full Jacobian matches ForwardDiff over a closure that recomputes
   `w_node` from state inside the AD pass — i.e., includes both direct
   and indirect chain rules.

The operator is wired into `Structured(...)`; these tests exercise it
standalone so local and indirect chain-rule regressions are easier to
attribute.
"""

using Test
using SparseArrays
using ForwardDiff
using Pyrolysis
using Pyrolysis.Materials: n_components, Material
using Pyrolysis.Discretization: compute_ale_species_advection!
using Pyrolysis.Problem: populate!, build_ad_workspace
using Pyrolysis.Residual: _compute_ale_mesh_velocity!
using Pyrolysis.Jacobian:
    ALEAdvectionOperator,
    ALEMeshVelocityOperator,
    build_ale_advection_workspace,
    build_ale_mesh_velocity_workspace,
    StructuredJacobian,
    build_residual_context,
    EmbeddedScaledPrefixSumCoupling,
    densify!

include(joinpath(@__DIR__, "..", "perf", "perf_helpers.jl"))

# Reach into the strain-rate AD helper (not in public exports).
const _strain_rate_derivatives! =
    getfield(Pyrolysis.Jacobian, :compute_strain_rate_derivatives!)

# Boundary flags used consistently for operator + references throughout this
# file (a venting top, impermeable bottom — the common production shape).
const _TOP_OPEN = true
const _BOTTOM_OPEN = false

# ---------------------------------------------------------------------------
# State setters
# ---------------------------------------------------------------------------

function _set_ramp_state!(u, layout, T_lo, T_hi, ξ_vals)
    n = layout.n_cells
    @inbounds for i = 1:n
        α = (i - 1) / max(1, n - 1)
        u[i] = T_lo * (1 - α) + T_hi * α
    end
    @inbounds for (j, ξj) in enumerate(ξ_vals)
        off = n + (j - 1) * n
        for i = 1:n
            u[off+i] = ξj
        end
    end
    return u
end

# ---------------------------------------------------------------------------
# Helper: refresh the state cache + properties + cross-section + geometry.
# Mirrors the prep done before any residual / Jacobian call.
# ---------------------------------------------------------------------------
function _refresh_state!(ws, u, def, material)
    populate!(ws.state_cache, u, ws.layout)
    Pyrolysis.update_cross_section_area!(
        ws.state_cache,
        material,
        ws.mesh,
        u,
        ws.layout,
        Pyrolysis.Problem.NoChi(),
    )
    Pyrolysis.update_state_geometry!(
        ws.state_cache,
        ws.mesh,
        Pyrolysis.Problem.mesh_mode(def.mode),
    )
    Pyrolysis.Properties.update_properties!(
        ws.prop_cache,
        ws.state_cache,
        ws.mesh,
        material,
    )
    return nothing
end

# ---------------------------------------------------------------------------
# Reference: full gas-transport contribution as a function of (u),
# *recomputing* w_node inside the closure so AD captures the indirect chain
# rule through the cumulative mesh velocity. Output is `[zeros(n_cells);
# advection_ξ_flat]` — the T rows are identically zero (no production energy
# mesh-motion advection); `advection_ξ_flat` is the column-major flatten of
# the gas-only species rows.
# ---------------------------------------------------------------------------
function _reference_full_advection_jacobian(u, ws, def, material::Material{NC}) where {NC}

    layout = ws.layout
    n_cells = layout.n_cells
    n_nodes = layout.n_nodes
    n_total = length(u)

    # AD config / workspace.
    out_len = n_cells + NC * n_cells   # T row + all ξ rows (zeros for non-gas)
    x0 = copy(u)
    f_buf = zeros(out_len)
    cfg =
        ForwardDiff.JacobianConfig(nothing, f_buf, x0, ForwardDiff.Chunk{min(8, n_total)}())
    Tu = eltype(cfg.duals[1])
    ws_ad = build_ad_workspace(def, ws, Tu)

    # Tu-typed scratch for w_node and advection.
    w_node_buf = zeros(Tu, n_nodes)
    θ_buf = zeros(Tu, n_cells)
    adv_ξ_buf = zeros(Tu, NC, n_cells)

    f! =
        let ws_ad = ws_ad,
            layout = layout,
            mesh = ws.mesh,
            material = material,
            def = def,
            w_node_buf = w_node_buf,
            θ_buf = θ_buf,
            adv_ξ_buf = adv_ξ_buf,
            out_len = out_len,
            n_cells = n_cells,
            NC = NC

            (out, x) -> begin
                populate!(ws_ad.state_cache, x, layout)
                Pyrolysis.update_cross_section_area!(
                    ws_ad.state_cache,
                    material,
                    mesh,
                    x,
                    layout,
                    Pyrolysis.Problem.NoChi(),
                )
                Pyrolysis.update_state_geometry!(
                    ws_ad.state_cache,
                    mesh,
                    Pyrolysis.Problem.mesh_mode(def.mode),
                )
                Pyrolysis.Properties.update_properties!(
                    ws_ad.prop_cache,
                    ws_ad.state_cache,
                    mesh,
                    material,
                )

                # Recompute w_node from state. Tu-typed so derivatives propagate.
                fill!(w_node_buf, zero(Tu))
                fill!(θ_buf, zero(Tu))
                _compute_ale_mesh_velocity!(
                    w_node_buf,
                    ws_ad.state_cache,
                    mesh,
                    ws_ad.prop_cache,
                    material,
                    0.0;
                    θ_out = θ_buf,
                    θ_A = zero(Tu),
                )

                # Flux-form gas transport (T rows are identically zero).
                fill!(adv_ξ_buf, zero(Tu))
                compute_ale_species_advection!(
                    adv_ξ_buf,
                    ws_ad.state_cache,
                    mesh,
                    material,
                    w_node_buf;
                    θ_A = zero(Tu),
                    top_gas_open = _TOP_OPEN,
                    bottom_gas_open = _BOTTOM_OPEN,
                )

                @inbounds for i = 1:n_cells
                    out[i] = zero(Tu)
                end
                @inbounds for j = 1:NC, i = 1:n_cells
                    out[n_cells+(j-1)*n_cells+i] = adv_ξ_buf[j, i]
                end
                return nothing
            end
        end

    J = ForwardDiff.jacobian(f!, f_buf, x0, cfg)
    return J
end

# ---------------------------------------------------------------------------
# Reference: direct-only chain rule. Holds `w_node` fixed (treated as a
# Float64 input, not Dual) so AD only captures the ξ stencil derivatives.
# Mirrors the operator's `contribute_to_jacobian!` direct entries in `jac.A`.
# ---------------------------------------------------------------------------
function _reference_direct_advection_jacobian(
    u,
    ws,
    def,
    material::Material{NC},
    w_node_fixed::AbstractVector,
) where {NC}

    layout = ws.layout
    n_cells = layout.n_cells
    n_total = length(u)

    out_len = n_cells + NC * n_cells
    x0 = copy(u)
    f_buf = zeros(out_len)
    cfg =
        ForwardDiff.JacobianConfig(nothing, f_buf, x0, ForwardDiff.Chunk{min(8, n_total)}())
    Tu = eltype(cfg.duals[1])
    ws_ad = build_ad_workspace(def, ws, Tu)

    # Promote the fixed w_node to Tu (acts as constants under AD, no Dual
    # perturbation).
    w_node_Tu = [Tu(w_node_fixed[k]) for k in eachindex(w_node_fixed)]

    adv_ξ_buf = zeros(Tu, NC, n_cells)

    f! =
        let ws_ad = ws_ad,
            layout = layout,
            mesh = ws.mesh,
            material = material,
            def = def,
            w_node_Tu = w_node_Tu,
            adv_ξ_buf = adv_ξ_buf,
            out_len = out_len,
            n_cells = n_cells,
            NC = NC

            (out, x) -> begin
                populate!(ws_ad.state_cache, x, layout)
                Pyrolysis.update_cross_section_area!(
                    ws_ad.state_cache,
                    material,
                    ws_ad.mesh,
                    x,
                    layout,
                    Pyrolysis.Problem.NoChi(),
                )
                Pyrolysis.update_state_geometry!(
                    ws_ad.state_cache,
                    ws_ad.mesh,
                    Pyrolysis.Problem.mesh_mode(def.mode),
                )

                fill!(adv_ξ_buf, zero(Tu))
                compute_ale_species_advection!(
                    adv_ξ_buf,
                    ws_ad.state_cache,
                    ws_ad.mesh,
                    material,
                    w_node_Tu;
                    θ_A = zero(Tu),
                    top_gas_open = _TOP_OPEN,
                    bottom_gas_open = _BOTTOM_OPEN,
                )

                @inbounds for i = 1:n_cells
                    out[i] = zero(Tu)
                end
                @inbounds for j = 1:NC, i = 1:n_cells
                    out[n_cells+(j-1)*n_cells+i] = adv_ξ_buf[j, i]
                end
                return nothing
            end
        end

    J = ForwardDiff.jacobian(f!, f_buf, x0, cfg)
    return J
end

# Permissive dense sparse pattern over the physics block for T + all ξ rows.
function _permissive_pattern(n_total, n_cells, NC)
    n_phys = n_cells * (1 + NC)
    I = Int[]
    J = Int[]
    V = Float64[]
    for i = 1:n_cells
        for c = 1:n_phys
            push!(I, i)
            push!(J, c)
            push!(V, 1.0)
        end
        for j = 1:NC
            ξ_row = n_cells + (j - 1) * n_cells + i
            for c = 1:n_phys
                push!(I, ξ_row)
                push!(J, c)
                push!(V, 1.0)
            end
        end
    end
    A = sparse(I, J, V, n_total, n_total)
    fill!(A.nzval, 0.0)
    return A
end

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@testset "ALEAdvectionOperator" begin

    @testset "interface declarations" begin
        op = ALEAdvectionOperator()
        @test Pyrolysis.Jacobian.reads(op) == (:w_cell, :theta)
        @test Pyrolysis.Jacobian.writes(op) == ()
        @test Pyrolysis.Jacobian.stencil(op) == (cells = -1:1,)
    end

    # Fixture: ALE NoChi PMMA, hot+gas state with a temperature ramp so θ
    # (and hence w) is non-trivial.
    def, ws, u, _ = build_ale_no_radiation()
    material = def.material
    layout = ws.layout
    mesh = ws.mesh
    NC = n_components(layout)
    n_cells = layout.n_cells
    n_nodes = layout.n_nodes

    # State: linear T ramp 600 → 1100 K with non-zero virgin + gas fractions
    # so reactions are active and θ is non-zero.
    _set_ramp_state!(u, layout, 600.0, 1100.0, (700.0, 1.5))
    _refresh_state!(ws, u, def, material)

    # Build the upstream ALEMeshVelocityOperator workspace and run apply! +
    # strain-rate derivatives. The advection operator needs `w_node` directly;
    # the indirect-chain-rule couplings flow through
    # `ctx.derived_fields[:w_cell]`, populated by the provider's
    # `contribute_to_structured!`.
    op_mv_ws = build_ale_mesh_velocity_workspace(material, n_cells, n_nodes)
    Pyrolysis.Jacobian.apply!(
        op_mv_ws,
        ALEMeshVelocityOperator(),
        ws.state_cache,
        mesh,
        material,
        0.0,
    )
    _strain_rate_derivatives!(op_mv_ws, ws.state_cache, material)
    ctx = build_residual_context([ALEMeshVelocityOperator(), ALEAdvectionOperator()])

    @testset "apply! reproduces residual gas transport; dT untouched" begin
        # Residual reference (gas species only — production has no energy
        # mesh-motion advection).
        adv_ξ_ref = zeros(NC, n_cells)
        compute_ale_species_advection!(
            adv_ξ_ref,
            ws.state_cache,
            mesh,
            material,
            op_mv_ws.w_node;
            θ_A = 0.0,
            top_gas_open = _TOP_OPEN,
            bottom_gas_open = _BOTTOM_OPEN,
        )

        # Operator output.
        dT_op = zeros(n_cells)
        dξ_op = zeros(NC, n_cells)
        Pyrolysis.Jacobian.apply!(
            dT_op,
            dξ_op,
            ALEAdvectionOperator(),
            ws.state_cache,
            mesh,
            op_mv_ws.w_node,
            material;
            θ_A = 0.0,
            top_gas_open = _TOP_OPEN,
            bottom_gas_open = _BOTTOM_OPEN,
        )

        @test all(iszero, dT_op)
        @test isapprox(dξ_op, adv_ξ_ref; rtol = 1e-12, atol = 1e-12)
    end

    @testset "gas conservation: fully closed column telescopes to zero" begin
        # With both boundaries impermeable, Σᵢ (advection + ξ·V̇/V)·Vᵢ = 0
        # exactly: the flux-form contribution must change total gas mass only
        # through boundary faces, of which there are none open.
        adv = zeros(NC, n_cells)
        compute_ale_species_advection!(
            adv,
            ws.state_cache,
            mesh,
            material,
            op_mv_ws.w_node;
            θ_A = 0.0,
            top_gas_open = false,
            bottom_gas_open = false,
        )
        A_section = ws.state_cache.A_section[]
        for j = 1:NC
            Pyrolysis.Materials.is_gas(material.components[j]) || continue
            # d(total gas mass)/dt from the mesh-motion term:
            # Σᵢ [advectionᵢ·Vᵢ + ξᵢ·V̇ᵢ] with V̇ᵢ = (w_R − w_L)·A.
            total = 0.0
            for i = 1:n_cells
                cell = mesh.cells[i]
                n_l, n_r = cell.node_ids
                V_i = ws.state_cache.volume[i]
                V̇_i = (op_mv_ws.w_node[n_r] - op_mv_ws.w_node[n_l]) * A_section
                total += adv[j, i] * V_i + ws.state_cache.ξ[j, i] * V̇_i
            end
            scale = maximum(abs, view(adv, j, :)) * maximum(ws.state_cache.volume)
            @test abs(total) ≤ 1e-12 * max(1.0, scale)
        end
    end

    @testset "contribute_to_jacobian! direct chain rule matches AD-w-fixed" begin
        n_total = length(u)
        n_phys = n_cells * (1 + NC)
        A = _permissive_pattern(n_total, n_cells, NC)
        jac = StructuredJacobian(:ale, A)

        op = ALEAdvectionOperator()
        ws_ad = build_ale_advection_workspace(material, n_cells, n_nodes)

        # Provider populates ctx.derived_fields[:w_cell] before consumer reads it.
        Pyrolysis.Jacobian.contribute_to_structured!(
            jac,
            ALEMeshVelocityOperator(),
            ws.state_cache,
            ctx,
            0.0,
            op_mv_ws,
            layout,
            Pyrolysis.Problem.NoChi(),
            n_total,
        )
        Pyrolysis.Jacobian.contribute_to_jacobian!(
            jac,
            op,
            ws.state_cache,
            mesh,
            material,
            layout,
            op_mv_ws.w_node,
            ctx,
            ws_ad;
            θ_A = 0.0,
            top_gas_open = _TOP_OPEN,
            bottom_gas_open = _BOTTOM_OPEN,
        )

        # DenseAD reference with w_node held fixed (= direct chain rule only).
        J_ref = _reference_direct_advection_jacobian(u, ws, def, material, op_mv_ws.w_node)

        # Compare jac.A on the (T, ξ) physics block.
        J_op = Matrix(jac.A)
        J_op_phys = J_op[1:n_phys, 1:n_phys]
        J_ref_phys = J_ref[1:n_phys, 1:n_phys]

        max_dev = maximum(abs, J_op_phys .- J_ref_phys)
        max_ref = maximum(abs, J_ref_phys)
        @test max_dev ≤ 1e-10 * max(1.0, max_ref)
    end

    @testset "contribute_to_jacobian! + densified embedded coupling = full AD" begin
        n_total = length(u)
        n_phys = n_cells * (1 + NC)
        A_sp = _permissive_pattern(n_total, n_cells, NC)
        jac = StructuredJacobian(:ale, A_sp)

        op = ALEAdvectionOperator()
        ws_ad = build_ale_advection_workspace(material, n_cells, n_nodes)

        # Provider populates ctx.derived_fields[:w_cell] before consumer reads it.
        Pyrolysis.Jacobian.contribute_to_structured!(
            jac,
            ALEMeshVelocityOperator(),
            ws.state_cache,
            ctx,
            0.0,
            op_mv_ws,
            layout,
            Pyrolysis.Problem.NoChi(),
            n_total,
        )
        Pyrolysis.Jacobian.contribute_to_jacobian!(
            jac,
            op,
            ws.state_cache,
            mesh,
            material,
            layout,
            op_mv_ws.w_node,
            ctx,
            ws_ad;
            θ_A = 0.0,
            top_gas_open = _TOP_OPEN,
            bottom_gas_open = _BOTTOM_OPEN,
        )

        # Combine A + Σⱼ structuredⱼ via densify into a single dense matrix.
        # This is exactly what `DirectSolve.ldiv!` does internally before
        # factoring.
        M = Matrix(jac.A)
        for c in jac.structured
            densify!(M, c)
        end

        # DenseAD reference with w_node recomputed inside the closure (=
        # direct + indirect chain rules together).
        J_ref = _reference_full_advection_jacobian(u, ws, def, material)

        # Compare on the (T, ξ) physics block.
        M_phys = M[1:n_phys, 1:n_phys]
        J_ref_phys = J_ref[1:n_phys, 1:n_phys]

        max_dev = maximum(abs, M_phys .- J_ref_phys)
        max_ref = maximum(abs, J_ref_phys)
        # AD-vs-analytic on the same formula: tight bound.
        @test max_dev ≤ 1e-8 * max(1.0, max_ref)
    end

    @testset "EmbeddedScaledPrefixSumCoupling densify! arithmetic" begin
        # Standalone unit test for the embedded coupling's densify, on a
        # small synthetic 4-cell setup.
        N = 4
        v = [1.0, 2.0, 3.0, 4.0]              # bare prefix-sum scalings
        n_left = [1, 2, 3, 4]                 # cell i's left node
        n_right = [2, 3, 4, 5]                # cell i's right node
        row_scale = [10.0, 20.0, 30.0, 40.0]
        # Expected per (i, j):
        #   factor = 0.5 * ((j ≤ n_left[i] - 1) + (j ≤ n_right[i] - 1))
        #   M[i, j] = row_scale[i] * factor * v[j]
        row_off = 0
        col_off = 0
        n_global = N

        c = EmbeddedScaledPrefixSumCoupling(
            v,
            n_left,
            n_right,
            row_scale,
            row_off,
            col_off,
            N,
            n_global,
        )

        M = zeros(N, N)
        densify!(M, c)

        # Hand-computed expectation:
        # n_left[1]=1 → (j ≤ 0) = 0;  n_right[1]=2 → (j ≤ 1) = 1 only at j=1
        #   factor row 1: [0.5, 0, 0, 0]
        # n_left[2]=2 → (j ≤ 1)=1@j=1; n_right[2]=3 → (j ≤ 2): j=1,2
        #   factor row 2: [(0.5+0.5), 0.5, 0, 0] = [1.0, 0.5, 0, 0]
        # n_left[3]=3 → (j ≤ 2): j=1,2; n_right[3]=4 → (j ≤ 3): j=1,2,3
        #   factor row 3: [1.0, 1.0, 0.5, 0]
        # n_left[4]=4 → (j ≤ 3): j=1,2,3; n_right[4]=5 → (j ≤ 4): j=1,2,3,4
        #   factor row 4: [1.0, 1.0, 1.0, 0.5]
        factor = [
            0.5 0.0 0.0 0.0
            1.0 0.5 0.0 0.0
            1.0 1.0 0.5 0.0
            1.0 1.0 1.0 0.5
        ]
        expected = [row_scale[i] * factor[i, j] * v[j] for i = 1:N, j = 1:N]

        @test isapprox(M, expected; rtol = 1e-14)

        # Same-node embedding (the flux-form consumer's per-node trick):
        # passing n_left twice must give the single-node indicator 𝟙[j ≤ n−1].
        c_same = EmbeddedScaledPrefixSumCoupling(
            v,
            n_left,
            n_left,
            row_scale,
            row_off,
            col_off,
            N,
            n_global,
        )
        M_same = zeros(N, N)
        densify!(M_same, c_same)
        factor_same = [(j <= n_left[i] - 1) ? 1.0 : 0.0 for i = 1:N, j = 1:N]
        expected_same = [row_scale[i] * factor_same[i, j] * v[j] for i = 1:N, j = 1:N]
        @test isapprox(M_same, expected_same; rtol = 1e-14)

        # Offset embedding: place the block at (row_off=2, col_off=3)
        # inside a 6×8 global matrix. Other entries must remain zero.
        c_emb = EmbeddedScaledPrefixSumCoupling(v, n_left, n_right, row_scale, 2, 3, N, 12)
        M_g = zeros(6, 8)
        densify!(M_g, c_emb)
        @test isapprox(M_g[3:6, 4:7], expected; rtol = 1e-14)
        # Outside the block: zero.
        M_g_block = copy(M_g)
        M_g_block[3:6, 4:7] .= 0.0
        @test maximum(abs, M_g_block) == 0.0
    end

    @testset "EmbeddedScaledPrefixSumCoupling apply! round-trip" begin
        # `apply!(y, c, x)` must agree with `densify!(M, c); y .+= M*x` —
        # this is the StructuredSolve-side gate that pins the
        # inline prefix-sum implementation.
        N = 6
        v = [0.7, -1.2, 0.3, 2.5, -0.4, 1.1]
        n_left = collect(1:N)
        n_right = collect(2:(N+1))
        row_scale = [1.5, -0.8, 2.2, 0.4, -1.7, 0.9]

        # Square global matrix at origin.
        c = EmbeddedScaledPrefixSumCoupling(v, n_left, n_right, row_scale, 0, 0, N, N)
        M = zeros(N, N)
        densify!(M, c)

        # Random x; apply! must accumulate (M*x) into y.
        x = randn(N)
        y_apply = zeros(N)
        Pyrolysis.Jacobian.apply!(y_apply, c, x)
        @test isapprox(y_apply, M * x; rtol = 1e-12, atol = 1e-12)

        # Additivity: applying twice into the same buffer doubles the result.
        Pyrolysis.Jacobian.apply!(y_apply, c, x)
        @test isapprox(y_apply, 2 .* (M * x); rtol = 1e-12, atol = 1e-12)

        # Embedded at non-trivial offsets inside a larger global vector.
        row_off = 5
        col_off = 9
        n_global = 20
        c_emb = EmbeddedScaledPrefixSumCoupling(
            v,
            n_left,
            n_right,
            row_scale,
            row_off,
            col_off,
            N,
            n_global,
        )
        M_g = zeros(n_global, n_global)
        densify!(M_g, c_emb)

        x_g = randn(n_global)
        y_g = zeros(n_global)
        Pyrolysis.Jacobian.apply!(y_g, c_emb, x_g)
        @test isapprox(y_g, M_g * x_g; rtol = 1e-12, atol = 1e-12)

        # Bounds errors: y too small at row_off / x too small at col_off.
        @test_throws DimensionMismatch Pyrolysis.Jacobian.apply!(
            zeros(row_off + N - 1),
            c_emb,
            x_g,
        )
        @test_throws DimensionMismatch Pyrolysis.Jacobian.apply!(
            y_g,
            c_emb,
            zeros(col_off + N - 1),
        )
    end
end
