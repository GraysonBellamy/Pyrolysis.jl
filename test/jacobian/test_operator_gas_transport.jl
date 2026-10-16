"""
`GasTransportOperator` (Fickian) correctness tests.

Same shape as the Reaction / Conduction tests:

1. **`apply!` equivalence.** The operator's per-cell `(S_conv / ρcp_eff,
   -div_J_gas[gas])` matches the gas-transport slice of `compute_rhs!`,
   reconstructed by directly running `compute_gas_fluxes!`,
   `compute_gas_convective_source!`, and `divergence_species_1d!` (no
   BC application — boundary face fluxes are zero, and the boundary
   cells' `S_conv` is zeroed: that whole energy row is owned by
   `BCConvectionOperator`, matching the operator's contract).

2. **`contribute_to_jacobian!` equivalence.** Per-cell tridiagonal AD
   matches `ForwardDiff.jacobian` over the gas-transport-only residual
   on the operator's declared support; off-pattern entries of the
   reference are zero (support-completeness).
"""

using Test
using SparseArrays
using ForwardDiff
using Pyrolysis
using Pyrolysis: to_problem_def, build_workspace
using Pyrolysis.Materials: n_components, Material, is_gas
using Pyrolysis.Properties: PropertyCache, update_properties!
using Pyrolysis.Physics: compute_gas_fluxes!, compute_gas_convective_source!
using Pyrolysis.Discretization: divergence_species_1d!
using Pyrolysis.Problem: populate!, build_ad_workspace
using Pyrolysis.Jacobian:
    GasTransportOperator, build_gas_transport_workspace, StructuredJacobian

include(joinpath(@__DIR__, "..", "perf", "perf_helpers.jl"))

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

# Reference gas-transport contribution: matches the *interior-flux-only*
# slice of compute_rhs!. Generic on state's scalar eltype.
function _reference_gas_transport(
    state,
    mesh,
    material::Material{NC},
    prop_cache,
    J_gas_buf,
    S_conv_buf,
    div_J_buf,
) where {NC}

    update_properties!(prop_cache, state, mesh, material)
    # Interior-only J_gas (boundary faces zero by compute_gas_fluxes!).
    compute_gas_fluxes!(J_gas_buf, state, mesh, prop_cache, material)
    # Cell-centered convective source from interior J_gas. The q_faces
    # argument feeds only the boundary cells' blended gradient, which this
    # operator no longer owns (BCConvectionOperator does) — pass zeros and
    # zero the boundary cells below to match the operator's contract.
    q_zero = zeros(eltype(S_conv_buf), length(state.A_face))
    compute_gas_convective_source!(
        S_conv_buf,
        J_gas_buf,
        q_zero,
        prop_cache.k,
        state,
        mesh,
        material,
    )
    S_conv_buf[1] = zero(eltype(S_conv_buf))
    S_conv_buf[end] = zero(eltype(S_conv_buf))
    # Species divergence on the same J_gas (no BC overwrite).
    divergence_species_1d!(div_J_buf, J_gas_buf, state, mesh)

    n_cells = length(state.T)
    Tu = eltype(state.T)
    dT = zeros(Tu, n_cells)
    dξ = zeros(Tu, NC, n_cells)
    for i = 1:n_cells
        dT[i] = S_conv_buf[i] / prop_cache.ρcp[i]
        for j = 1:NC
            if is_gas(material.components[j])
                dξ[j, i] = -div_J_buf[j, i]
            end
        end
    end
    return dT, dξ
end

function _reference_gas_transport_jacobian(u, ws, def, ::Material{NC}) where {NC}
    layout = ws.layout
    mesh = ws.mesh
    material = def.material
    n_cells = layout.n_cells
    n_total = length(u)
    n_phys = n_cells * (1 + NC)   # T + ξ block

    x0 = copy(u)
    # Output: stack (dT, dξ_gas_for_each_j_irrespective_of_phase): n_cells * (1 + NC)
    out_dim = n_cells * (1 + NC)
    f_buf = zeros(out_dim)
    cfg = ForwardDiff.JacobianConfig(
        nothing,
        f_buf,
        x0,
        ForwardDiff.Chunk{min(12, n_total)}(),
    )
    Tu = eltype(cfg.duals[1])
    ws_ad = build_ad_workspace(def, ws, Tu)

    J_gas_buf = zeros(Tu, NC, length(mesh.faces))
    S_conv_buf = zeros(Tu, n_cells)
    div_J_buf = zeros(Tu, NC, n_cells)

    f! =
        let ws_ad = ws_ad,
            layout = layout,
            mesh = mesh,
            material = material,
            def = def,
            J_gas_buf = J_gas_buf,
            S_conv_buf = S_conv_buf,
            div_J_buf = div_J_buf

            (out, x) -> begin
                populate!(ws_ad.state_cache, x, layout)
                Pyrolysis.Problem.update_state_geometry!(
                    ws_ad.state_cache,
                    mesh,
                    Pyrolysis.Problem.mesh_mode(def.mode),
                )
                dT, dξ = _reference_gas_transport(
                    ws_ad.state_cache,
                    mesh,
                    material,
                    ws_ad.prop_cache,
                    J_gas_buf,
                    S_conv_buf,
                    div_J_buf,
                )
                @inbounds for i = 1:n_cells
                    out[i] = dT[i]
                end
                @inbounds for j = 1:NC
                    off = n_cells + (j - 1) * n_cells
                    for i = 1:n_cells
                        out[off+i] = dξ[j, i]
                    end
                end
                return nothing
            end
        end

    J = ForwardDiff.jacobian(f!, f_buf, x0, cfg)
    return J
end

@testset "GasTransportOperator (Fickian)" begin

    def, ws, u, _ = build_baseline_eulerian()
    n = ws.layout.n_cells
    NC = n_components(ws.layout)
    material = def.material

    # Reference scratch
    prop_cache_ref = PropertyCache(n, length(ws.mesh.faces))
    J_gas_ref = zeros(NC, length(ws.mesh.faces))
    S_conv_ref = zeros(n)
    div_J_ref = zeros(NC, n)

    @testset "apply! matches gas-transport residual slice — ramp+gas state" begin
        _set_state_ramp!(u, n, 400.0, 1200.0, (700.0, 5.0))
        populate!(ws.state_cache, u, ws.layout)
        Pyrolysis.Problem.update_state_geometry!(
            ws.state_cache,
            ws.mesh,
            Pyrolysis.Problem.mesh_mode(def.mode),
        )

        dT_ref, dξ_ref = _reference_gas_transport(
            ws.state_cache,
            ws.mesh,
            material,
            prop_cache_ref,
            J_gas_ref,
            S_conv_ref,
            div_J_ref,
        )

        op = GasTransportOperator()
        dT_op = zeros(n)
        dξ_op = zeros(NC, n)
        Pyrolysis.Jacobian.apply!(dT_op, dξ_op, op, ws.state_cache, ws.mesh, material)

        @test isapprox(dT_op, dT_ref; rtol = 1e-12, atol = 1e-12)
        @test isapprox(dξ_op, dξ_ref; rtol = 1e-12, atol = 1e-12)
    end

    @testset "contribute_to_jacobian! matches DenseAD on declared support" begin
        _set_state_ramp!(u, n, 400.0, 1200.0, (700.0, 5.0))

        # Build operator's declared support pattern. T row + every ξ row
        # (including solids — solid rows scatter zeros, but the pattern
        # slots exist).
        n_total = length(u)
        I = Int[]
        J_idx = Int[]
        V = Float64[]
        for i = 1:n
            iL = i == 1 ? 0 : i - 1
            iR = i == n ? 0 : i + 1
            # Block rows: T (1) + ξ_j (1..NC). Each row writes 3 T cols + 3 ξ cols × NC.
            rows = (i, ntuple(j -> n + (j - 1) * n + i, NC)...)
            for r in rows
                if iL != 0

                    push!(I, r)
                    push!(J_idx, iL)
                    push!(V, 1.0)
                end
                push!(I, r)
                push!(J_idx, i)
                push!(V, 1.0)
                if iR != 0

                    push!(I, r)
                    push!(J_idx, iR)
                    push!(V, 1.0)
                end
                for k = 1:NC
                    ξ_off = n + (k - 1) * n
                    if iL != 0

                        push!(I, r)
                        push!(J_idx, ξ_off + iL)
                        push!(V, 1.0)
                    end
                    push!(I, r)
                    push!(J_idx, ξ_off + i)
                    push!(V, 1.0)
                    if iR != 0

                        push!(I, r)
                        push!(J_idx, ξ_off + iR)
                        push!(V, 1.0)
                    end
                end
            end
        end
        A = sparse(I, J_idx, V, n_total, n_total)
        fill!(A.nzval, 0.0)
        jac = StructuredJacobian(A)

        op = GasTransportOperator()
        ws_ad = build_gas_transport_workspace(material)

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
            0.0,
            ws_ad,
            ws.state_cache,
            ws.mesh,
            material,
            ws.layout,
        )

        # Reference Jacobian (full physics block).
        J_ref = _reference_gas_transport_jacobian(u, ws, def, material)

        J_op = Matrix(jac.A)
        max_dev = 0.0
        for (r, c) in zip(I, J_idx)
            diff = abs(J_op[r, c] - J_ref[r, c])
            max_dev = max(max_dev, diff)
        end
        @test max_dev < 1e-9

        # Support-completeness on the (T, ξ) physics block.
        n_phys = n * (1 + NC)
        Pmask = falses(n_phys, n_phys)
        for (r, c) in zip(I, J_idx)
            if r <= n_phys && c <= n_phys
                Pmask[r, c] = true
            end
        end
        @test maximum(abs.(J_ref[1:n_phys, 1:n_phys][.!Pmask])) < 1e-9
    end
end
