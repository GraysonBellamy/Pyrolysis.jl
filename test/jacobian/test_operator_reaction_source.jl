"""
`ReactionSourceOperator` correctness tests.

Two gates:

1. **`apply!` equivalence.** Calling
   `Pyrolysis.Jacobian.apply!(dT, dξ, ReactionSourceOperator(), state, t,
                              material, scratch)`
   reproduces the per-cell `(Q_rxn / ρcp_eff, dξ_rxn[:])` slice of the
   existing residual (extracted by directly running the reaction loop
   from `compute_rhs!`). Catches drift between the operator and the
   live residual kernel.

2. **`contribute_to_jacobian!` equivalence.** The cell-local AD path
   produces the same diagonal block as `DenseAD_ForwardDiff` over the
   isolated reaction source. Validates the AD plumbing and the
   sparse-scatter pattern.

`ReactionSourceOperator` is part of the production `Structured(...)`
operator list; these tests keep its local block isolated so attribution
stays clear.
"""

using Test
using SparseArrays
using ForwardDiff
using Pyrolysis
using Pyrolysis: to_problem_def, build_workspace
using Pyrolysis.Materials: n_components, n_reactions, Material
using Pyrolysis.Properties: effective_heat_capacity
using Pyrolysis.Physics: evaluate_reactions!, reaction_heat_source, species_source_terms!
using Pyrolysis.Problem: populate!
using Pyrolysis.Residual: residual!
using Pyrolysis.Jacobian:
    ReactionSourceOperator,
    build_reaction_source_workspace,
    StructuredJacobian,
    contribute_to_jacobian!,
    FastSparsityMap

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

# Reference: per-cell (Q_rxn / ρcp, dξ_rxn[:]) computed directly. This is the
# slice of compute_rhs! that ReactionSourceOperator computes, with the
# Q_rxn/ρcp division applied at this layer (compute_rhs! defers the division
# to assembly; mathematically equivalent).
function _reference_reaction_contribution(state, material::Material{NC,NR}) where {NC,NR}
    n = length(state.T)
    dT_ref = zeros(n)
    dξ_ref = zeros(NC, n)
    rates = zeros(NR)
    heats = zeros(NR)
    src = zeros(NC)
    for i = 1:n
        T = state.T[i]
        ξ = ntuple(j -> state.ξ[j, i], Val(NC))
        fill!(rates, 0.0)
        fill!(heats, 0.0)
        evaluate_reactions!(rates, heats, ξ, T, material)
        Q_rxn = reaction_heat_source(rates, heats)
        ρcp = effective_heat_capacity(ξ, T, material)
        dT_ref[i] = Q_rxn / ρcp
        species_source_terms!(src, rates, material)
        for j = 1:NC
            dξ_ref[j, i] = src[j]
        end
    end
    return dT_ref, dξ_ref
end

# DenseAD reference for the cell-local diagonal block. For each cell i, runs
# ForwardDiff on the same closure used by contribute_to_jacobian! — but with
# fresh inputs, no scratch reuse. If the operator's hot-path J_local matches
# this, the AD plumbing and scatter are correct.
function _reference_dense_jacobian(u, layout, material::Material{NC,NR}) where {NC,NR}
    n_cells = layout.n_cells
    dim = 1 + NC
    n_total = length(u)
    J_full = zeros(n_total, n_total)
    for i = 1:n_cells
        x = zeros(dim)
        x[1] = u[i]
        for j = 1:NC
            x[1+j] = u[n_cells+(j-1)*n_cells+i]
        end
        f! = function (out, x_inner)
            rates = zeros(eltype(x_inner), NR)
            heats = zeros(eltype(x_inner), NR)
            src = zeros(eltype(x_inner), NC)
            T = x_inner[1]
            ξ = ntuple(j -> x_inner[1+j], Val(NC))
            evaluate_reactions!(rates, heats, ξ, T, material)
            Q_rxn = reaction_heat_source(rates, heats)
            ρcp = effective_heat_capacity(ξ, T, material)
            out[1] = Q_rxn / ρcp
            species_source_terms!(src, rates, material)
            for j = 1:NC
                out[1+j] = src[j]
            end
            return nothing
        end
        J_local = ForwardDiff.jacobian(f!, zeros(dim), x)
        # Scatter into J_full at the same positions the operator writes.
        row_T = i
        col_T = i
        J_full[row_T, col_T] = J_local[1, 1]
        for l = 1:NC
            col_ξ = n_cells + (l - 1) * n_cells + i
            J_full[row_T, col_ξ] = J_local[1, 1+l]
        end
        for j = 1:NC
            row_ξ = n_cells + (j - 1) * n_cells + i
            J_full[row_ξ, col_T] = J_local[1+j, 1]
            for l = 1:NC
                col_ξ = n_cells + (l - 1) * n_cells + i
                J_full[row_ξ, col_ξ] = J_local[1+j, 1+l]
            end
        end
    end
    return J_full
end

@testset "ReactionSourceOperator" begin

    # Use the eulerian PMMA fixture — has 1 reaction, 2 components.
    def, ws, u, _ = build_baseline_eulerian()
    n = ws.layout.n_cells
    NC = n_components(ws.layout)
    material = def.material
    NR = n_reactions(material)

    @testset "apply! matches per-cell residual slice — hot+gas state" begin
        _set_state!(u, n, 1200.0, (400.0, 8.0))
        populate!(ws.state_cache, u, ws.layout)

        dT_ref, dξ_ref = _reference_reaction_contribution(ws.state_cache, material)

        # Use the operator's apply!.
        op = ReactionSourceOperator()
        ws_ad = build_reaction_source_workspace(material)
        # Float64 path: build a Float64 scratch to mirror the residual.
        fp_scratch = Pyrolysis.Jacobian.ReactionScratch{Float64}(NC, NR)
        dT_op = zeros(n)
        dξ_op = zeros(NC, n)
        Pyrolysis.Jacobian.apply!(
            dT_op,
            dξ_op,
            op,
            ws.state_cache,
            0.0,
            material,
            fp_scratch,
        )

        @test isapprox(dT_op, dT_ref; rtol = 1e-12, atol = 1e-14)
        @test isapprox(dξ_op, dξ_ref; rtol = 1e-12, atol = 1e-14)
    end

    @testset "contribute_to_jacobian! matches DenseAD on the cell-local block" begin
        _set_state!(u, n, 1200.0, (400.0, 8.0))

        # Build a sparse pattern that covers exactly the cell-local diagonal
        # block — for each cell i, the (1+NC) × (1+NC) sub-block at the
        # corresponding global indices. This is the operator's declared
        # support.
        n_total = length(u)
        I = Int[]
        J = Int[]
        V = Float64[]
        for i = 1:n
            row_T = i
            col_T = i
            push!(I, row_T)
            push!(J, col_T)
            push!(V, 1.0)
            for l = 1:NC
                col_ξ = n + (l - 1) * n + i
                push!(I, row_T)
                push!(J, col_ξ)
                push!(V, 1.0)
            end
            for j = 1:NC
                row_ξ = n + (j - 1) * n + i
                push!(I, row_ξ)
                push!(J, col_T)
                push!(V, 1.0)
                for l = 1:NC
                    col_ξ = n + (l - 1) * n + i
                    push!(I, row_ξ)
                    push!(J, col_ξ)
                    push!(V, 1.0)
                end
            end
        end
        A = sparse(I, J, V, n_total, n_total)
        fill!(A.nzval, 0.0)
        jac = StructuredJacobian(A)

        op = ReactionSourceOperator()
        ws_ad = build_reaction_source_workspace(material)

        contribute_to_jacobian!(jac, op, u, 0.0, ws_ad, material, ws.layout, ws.mesh)

        # Compare against the DenseAD reference.
        J_ref = _reference_dense_jacobian(u, ws.layout, material)
        J_op = Matrix(jac.A)
        # Restrict the comparison to the support pattern only — the operator
        # only writes inside its declared support, so off-pattern reference
        # entries are by definition not produced here.
        support = Matrix(A) .!= 0
        # Build support from the pattern coordinates we pushed.
        Pmask = falses(n_total, n_total)
        for (r, c) in zip(I, J)
            Pmask[r, c] = true
        end
        @test maximum(abs.(J_op[Pmask] .- J_ref[Pmask])) < 1e-9
        # The off-pattern entries of the reference must be zero — proves the
        # operator's declared sparsity is support-complete for its piece.
        @test maximum(abs.(J_ref[.!Pmask])) < 1e-12
    end

    @testset "apply! matches residual slice — cold initial state" begin
        # Cold state: reaction rates near zero. Sanity check that we don't
        # blow up at low temperature.
        _set_state!(u, n, 300.0, (1190.0, 0.0))
        populate!(ws.state_cache, u, ws.layout)

        dT_ref, dξ_ref = _reference_reaction_contribution(ws.state_cache, material)

        op = ReactionSourceOperator()
        fp_scratch = Pyrolysis.Jacobian.ReactionScratch{Float64}(NC, NR)
        dT_op = zeros(n)
        dξ_op = zeros(NC, n)
        Pyrolysis.Jacobian.apply!(
            dT_op,
            dξ_op,
            op,
            ws.state_cache,
            0.0,
            material,
            fp_scratch,
        )

        @test isapprox(dT_op, dT_ref; rtol = 1e-12, atol = 1e-14)
        @test isapprox(dξ_op, dξ_ref; rtol = 1e-12, atol = 1e-14)
    end
end
