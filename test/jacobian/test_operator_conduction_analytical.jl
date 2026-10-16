"""
Analytical override equivalence gate for `ConductionOperator`.

Asserts that `_contribute_to_jacobian_analytical!` (hand-derived chain
rule) and `_contribute_to_jacobian_ad!` (ForwardDiff over the 3-cell
stencil) produce the same Jacobian rows at rtol = 1e-12 on the
canonical configurations × probe states. Only floating-point summation
order differs between the two paths; the underlying math is the same.

If this test fails, the analytical override has drifted from the
residual's `apply_cell!` and must be fixed before it ships as the
public `contribute_to_jacobian!` default.
"""

using Test
using SparseArrays
using ForwardDiff
using Pyrolysis
using Pyrolysis.Materials: Material
using Pyrolysis.Jacobian:
    StructuredJacobian,
    ConductionOperator,
    ConductionADWorkspace,
    build_conduction_workspace,
    JacobianSpec,
    Structured,
    DirectSolve,
    build_sparsity,
    FastSparsityMap
using Pyrolysis.Jacobian: _contribute_to_jacobian_analytical!, _contribute_to_jacobian_ad!
using Pyrolysis.Problem:
    ProblemDef,
    Workspace,
    populate!,
    update_state_geometry!,
    update_cross_section_area!,
    mesh_mode,
    chi_mode

include(joinpath(@__DIR__, "..", "perf", "perf_helpers.jl"))

# --- Test helper: build a fresh StructuredJacobian + workspace, run one of
#     the two paths, return the assembled J (sparse) ---
function _assemble_conduction(def, ws, u, t, path::Symbol)
    S_phys, S_geom = build_sparsity(ws.layout, def.mode)
    S_union = S_phys .| S_geom
    A = SparseMatrixCSC{Float64,Int}(
        S_union.m,
        S_union.n,
        copy(S_union.colptr),
        copy(S_union.rowval),
        zeros(Float64, length(S_union.nzval)),
    )
    jac = StructuredJacobian(A)

    cond_ws = build_conduction_workspace(def.material)

    # Populate state cache & geometry (mirrors what
    # `Structured._compute_jacobian!` does before walking operators).
    populate!(ws.state_cache, u, ws.layout)
    if mesh_mode(def.mode) isa Pyrolysis.Problem.ALE
        update_cross_section_area!(
            ws.state_cache,
            def.material,
            ws.mesh,
            u,
            ws.layout,
            chi_mode(def.mode),
        )
        update_state_geometry!(ws.state_cache, ws.mesh, mesh_mode(def.mode))
    end

    if path === :analytical
        _contribute_to_jacobian_analytical!(
            jac,
            ConductionOperator(),
            u,
            Float64(t),
            cond_ws,
            ws.state_cache,
            ws.mesh,
            def.material,
            ws.layout,
        )
    elseif path === :ad
        _contribute_to_jacobian_ad!(
            jac,
            ConductionOperator(),
            u,
            Float64(t),
            cond_ws,
            ws.state_cache,
            ws.mesh,
            def.material,
            ws.layout,
        )
    else
        error("unknown path: $path")
    end
    return jac.A
end

# --- Probe states ---
function _ramp_state!(u, n_cells, NC, T_lo, T_hi, ξ_vals)
    @inbounds for i = 1:n_cells
        α = (i - 1) / max(1, n_cells - 1)
        u[i] = T_lo * (1 - α) + T_hi * α
    end
    @inbounds for (j, ξj) in enumerate(ξ_vals)
        off = n_cells + (j - 1) * n_cells
        for i = 1:n_cells
            u[off+i] = ξj
        end
    end
    return u
end

function _hot_uniform_state!(u, n_cells, NC, T_uniform, ξ_vals)
    fill!(view(u, 1:n_cells), T_uniform)
    @inbounds for (j, ξj) in enumerate(ξ_vals)
        off = n_cells + (j - 1) * n_cells
        for i = 1:n_cells
            u[off+i] = ξj
        end
    end
    return u
end

@testset "ConductionOperator analytical override == AD-default" begin

    probe_states = [
        ("default IC", (n, NC, u) -> u),  # whatever the builder gives
        ("uniform 600K", (n, NC, u) -> _hot_uniform_state!(u, n, NC, 600.0, (700.0, 1.5))),
        (
            "uniform 1100K",
            (n, NC, u) -> _hot_uniform_state!(u, n, NC, 1100.0, (200.0, 0.5)),
        ),
        (
            "ramp 600→1100",
            (n, NC, u) -> _ramp_state!(u, n, NC, 600.0, 1100.0, (800.0, 1.5)),
        ),
        (
            "ramp 300→1500",
            (n, NC, u) -> _ramp_state!(u, n, NC, 300.0, 1500.0, (500.0, 1.0)),
        ),
    ]

    configs = [
        ("Eulerian PMMA", build_baseline_eulerian),
        ("Eulerian Wood Darcy", build_baseline_eulerian_darcy),
        ("Eulerian SurfAbs", build_eulerian_surface_absorption),
        ("ALE PMMA", build_ale_no_radiation),
    ]

    for (cfg_name, builder) in configs
        def, ws, u_base, _ = builder()
        n_cells = ws.layout.n_cells
        NC = length(def.material.components)

        @testset "$cfg_name" begin
            for (state_name, set_state!) in probe_states
                u = copy(u_base)
                set_state!(n_cells, NC, u)

                J_ana = _assemble_conduction(def, ws, u, 0.0, :analytical)
                J_ad = _assemble_conduction(def, ws, u, 0.0, :ad)

                # Compare nzval directly — same sparsity pattern.
                diff_abs = maximum(abs, J_ana.nzval .- J_ad.nzval)
                scale = max(maximum(abs, J_ad.nzval), 1e-30)
                rel = diff_abs / scale

                @test rel ≤ 1e-12 || (
                    @info "Conduction ana vs ad" cfg_name state_name diff_abs scale rel;
                    false
                )
            end
        end
    end
end
