"""
`ConductionOperator` correctness tests.

Two gates (mirroring the ReactionSource tests):

1. **`apply!` equivalence.** Calling
   `Pyrolysis.Jacobian.apply!(dT, ConductionOperator(), state, mesh, material)`
   reproduces `-divergence(q_cond) / ρcp_eff` for every cell, computed by
   directly running the existing `compute_conductive_fluxes!` and
   `divergence_1d!` in isolation (no BC application). Catches drift between
   the operator and the live conduction kernel.

2. **`contribute_to_jacobian!` equivalence.** The per-cell tridiagonal-AD
   path produces the same row-i Jacobian entries as `ForwardDiff.jacobian`
   over the conduction-only residual function — restricted to the operator's
   declared support pattern.

The operator is **not yet wired into `Structured(...)`** — it is exercised
standalone (the swap is deferred until all 10 operators land).
"""

using Test
using SparseArrays
using ForwardDiff
using Pyrolysis
using Pyrolysis: to_problem_def, build_workspace
using Pyrolysis.Materials: n_components, Material
using Pyrolysis.Properties: PropertyCache, update_properties!, effective_heat_capacity
using Pyrolysis.Physics: compute_conductive_fluxes!
using Pyrolysis.Discretization: divergence_1d!
using Pyrolysis.Problem: populate!, build_ad_workspace
using Pyrolysis.Jacobian:
    ConductionOperator,
    build_conduction_workspace,
    StructuredJacobian,
    build_residual_context

include(joinpath(@__DIR__, "..", "perf", "perf_helpers.jl"))

function _set_state!(u, n_cells, T_val, ξ_vals)
    @inbounds for i = 1:n_cells
        u[i] = T_val
    end
    @inbounds for (j, ξj) in enumerate(ξ_vals)
        off = n_cells + (j - 1) * n_cells
        for i = 1:n_cells
            u[off+i] = ξj
        end
    end
    return u
end

# Set a *spatially varying* T profile so conduction is non-zero (the uniform-T
# state setter above makes conduction identically zero).
function _set_state_ramp!(u, n_cells, T_lo, T_hi, ξ_vals)
    for i = 1:n_cells
        α = (i - 1) / max(1, n_cells - 1)
        u[i] = T_lo * (1 - α) + T_hi * α
    end
    for (j, ξj) in enumerate(ξ_vals)
        off = n_cells + (j - 1) * n_cells
        for i = 1:n_cells
            u[off+i] = ξj
        end
    end
    return u
end

# Reference: -div(q_cond) / ρcp_eff per cell, mirroring the conduction-only
# slice of compute_rhs! (BCs *not* applied). Generic on the state's scalar
# eltype — works for both Float64 (residual-slice equivalence) and Dual
# (AD reference for the Jacobian).
function _reference_conduction_dT(
    state,
    mesh,
    material::Material{NC},
    prop_cache,
    q_cond,
    div_q,
) where {NC}
    update_properties!(prop_cache, state, mesh, material)
    fill!(q_cond, zero(eltype(q_cond)))
    compute_conductive_fluxes!(q_cond, state, mesh, prop_cache)
    fill!(div_q, zero(eltype(div_q)))
    divergence_1d!(div_q, q_cond, state, mesh)

    n_cells = length(state.T)
    dT_ref = Vector{eltype(state.T)}(undef, n_cells)
    for i = 1:n_cells
        dT_ref[i] = -div_q[i] / prop_cache.ρcp[i]
    end
    return dT_ref
end

# Reference Jacobian: ∂dT_cond / ∂u via ForwardDiff over a function that
# computes the conduction-only contribution for the whole vector. The function
# accepts the *full* state vector u (length n_total) so the AD output's
# columns map directly to global state indices — no per-block remapping.
function _reference_conduction_jacobian(u, ws, def, ::Material{NC}) where {NC}
    layout = ws.layout
    mesh = ws.mesh
    material = def.material
    n_cells = layout.n_cells
    n_total = length(u)

    # Build a Tu workspace for the AD pass. We use ForwardDiff.JacobianConfig
    # explicitly so build_ad_workspace can be called once with the right Tu.
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

    # Per-call scratch (Tu-typed). Reused via closure capture.
    q_cond_buf = zeros(Tu, length(mesh.faces))
    div_q_buf = zeros(Tu, n_cells)

    f! =
        let ws_ad = ws_ad,
            layout = layout,
            mesh = mesh,
            material = material,
            def = def,
            q_cond_buf = q_cond_buf,
            div_q_buf = div_q_buf

            (out, x) -> begin
                populate!(ws_ad.state_cache, x, layout)
                Pyrolysis.Problem.update_state_geometry!(
                    ws_ad.state_cache,
                    mesh,
                    Pyrolysis.Problem.mesh_mode(def.mode),
                )
                dT_ref = _reference_conduction_dT(
                    ws_ad.state_cache,
                    mesh,
                    material,
                    ws_ad.prop_cache,
                    q_cond_buf,
                    div_q_buf,
                )
                @inbounds for i = 1:n_cells
                    out[i] = dT_ref[i]
                end
                return nothing
            end
        end

    J = ForwardDiff.jacobian(f!, f_buf, x0, cfg)
    return J
end

@testset "ConductionOperator" begin

    @testset "interface declarations" begin
        op = ConductionOperator()
        @test Pyrolysis.Jacobian.reads(op) == ()
        @test Pyrolysis.Jacobian.writes(op) == ()
        @test Pyrolysis.Jacobian.stencil(op) == (cells = -1:1,)
    end

    # ----------------------------------------------------------
    # Eulerian PMMA fixture
    # ----------------------------------------------------------
    def, ws, u, _ = build_baseline_eulerian()
    n = ws.layout.n_cells
    NC = n_components(ws.layout)
    material = def.material

    # Buffers for the residual-slice reference.
    prop_cache_ref = PropertyCache(n, length(ws.mesh.faces))
    q_cond_ref = zeros(length(ws.mesh.faces))
    div_q_ref = zeros(n)

    @testset "apply! matches -div(q_cond)/ρcp on ramp state" begin
        _set_state_ramp!(u, n, 400.0, 1200.0, (700.0, 5.0))
        populate!(ws.state_cache, u, ws.layout)
        Pyrolysis.Problem.update_state_geometry!(
            ws.state_cache,
            ws.mesh,
            Pyrolysis.Problem.mesh_mode(def.mode),
        )

        dT_ref = _reference_conduction_dT(
            ws.state_cache,
            ws.mesh,
            material,
            prop_cache_ref,
            q_cond_ref,
            div_q_ref,
        )

        op = ConductionOperator()
        dT_op = zeros(n)
        Pyrolysis.Jacobian.apply!(dT_op, op, ws.state_cache, ws.mesh, material)

        @test isapprox(dT_op, dT_ref; rtol = 1e-12, atol = 1e-12)
    end

    @testset "apply! is identically zero on uniform state" begin
        # Uniform T → no thermal gradient → conduction zero everywhere.
        _set_state!(u, n, 1000.0, (700.0, 5.0))
        populate!(ws.state_cache, u, ws.layout)

        op = ConductionOperator()
        dT_op = zeros(n)
        Pyrolysis.Jacobian.apply!(dT_op, op, ws.state_cache, ws.mesh, material)

        @test maximum(abs, dT_op) < 1e-12
    end

    @testset "contribute_to_jacobian! matches DenseAD on T-row × stencil" begin
        _set_state_ramp!(u, n, 400.0, 1200.0, (700.0, 5.0))

        # Build a sparse pattern covering the operator's declared support:
        # for each cell i, columns are T at {iL, i, iR} and ξ at {iL, i, iR}
        # (skipping missing-neighbor slots at boundaries). Only T rows.
        n_total = length(u)
        I = Int[]
        J = Int[]
        V = Float64[]
        for i = 1:n
            iL = i == 1 ? 0 : i - 1
            iR = i == n ? 0 : i + 1
            # T columns
            if iL != 0

                push!(I, i)
                push!(J, iL)
                push!(V, 1.0)
            end
            push!(I, i)
            push!(J, i)
            push!(V, 1.0)
            if iR != 0

                push!(I, i)
                push!(J, iR)
                push!(V, 1.0)
            end
            # ξ columns
            for j = 1:NC
                ξ_off = n + (j - 1) * n
                if iL != 0

                    push!(I, i)
                    push!(J, ξ_off + iL)
                    push!(V, 1.0)
                end
                push!(I, i)
                push!(J, ξ_off + i)
                push!(V, 1.0)
                if iR != 0

                    push!(I, i)
                    push!(J, ξ_off + iR)
                    push!(V, 1.0)
                end
            end
        end
        A = sparse(I, J, V, n_total, n_total)
        fill!(A.nzval, 0.0)
        jac = StructuredJacobian(:eulerian, A)

        op = ConductionOperator()
        ws_ad = build_conduction_workspace(material)
        ctx = build_residual_context([op])

        # Populate state cache for the operator.
        populate!(ws.state_cache, u, ws.layout)
        Pyrolysis.Problem.update_state_geometry!(
            ws.state_cache,
            ws.mesh,
            Pyrolysis.Problem.mesh_mode(def.mode),
        )

        Pyrolysis.Jacobian.contribute_to_jacobian!(
            jac,
            op,
            u,
            ctx,
            0.0,
            ws_ad,
            ws.state_cache,
            ws.mesh,
            material,
            ws.layout,
        )

        # DenseAD reference: J_ref is (n_cells × n_total) — rows are T-row
        # Jacobian entries; columns are global state indices.
        J_ref_full = _reference_conduction_jacobian(u, ws, def, material)

        J_op = Matrix(jac.A)
        max_dev = 0.0
        for (r, c) in zip(I, J)
            diff = abs(J_op[r, c] - J_ref_full[r, c])
            max_dev = max(max_dev, diff)
        end
        @test max_dev < 1e-9

        # Support-completeness: any non-zero (T, ξ) entry of J_ref must be
        # inside the operator's declared support pattern. Compare only the
        # (T, ξ) columns (the first n_cells * (1 + NC) columns).
        n_phys = n * (1 + NC)
        Pmask = falses(n, n_phys)
        for (r, c) in zip(I, J)
            if r <= n && c <= n_phys
                Pmask[r, c] = true
            end
        end
        @test maximum(abs.(J_ref_full[1:n, 1:n_phys][.!Pmask])) < 1e-9
    end
end
