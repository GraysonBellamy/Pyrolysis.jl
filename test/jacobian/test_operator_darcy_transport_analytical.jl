"""
Analytical override equivalence gate for `DarcyTransportOperator`.

Asserts `_contribute_to_jacobian_analytical!` matches
`_contribute_to_jacobian_ad!` at rtol = 1e-12 on the canonical probe
states, for both `build_baseline_eulerian_darcy()` and
`build_ale_beer_lambert_darcy()` configs.

Probe ξ values are chosen so `Σⱼ ξⱼ_solid / ρⱼ_intrinsic ≤ 1`.
For the wood-like material the solid is component 1 with ρ_intrinsic = 500,
so ξ_solid ≤ 500 keeps porosity above the floor and `dP/dξ` analytical.
"""

using Test
using SparseArrays
using ForwardDiff
using Pyrolysis
using Pyrolysis.Materials: Material
using Pyrolysis.Jacobian:
    StructuredJacobian,
    DarcyTransportOperator,
    ResidualContext,
    build_residual_context,
    DarcyTransportADWorkspace,
    build_darcy_transport_workspace,
    build_sparsity
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

function _assemble_darcy(def, ws, u, t, path::Symbol)
    S_phys, S_geom = build_sparsity(ws.layout, def.mode)
    S_union = S_phys .| S_geom
    A = SparseMatrixCSC{Float64,Int}(
        S_union.m,
        S_union.n,
        copy(S_union.colptr),
        copy(S_union.rowval),
        zeros(Float64, length(S_union.nzval)),
    )
    mode_tag = mesh_mode(def.mode) isa Pyrolysis.Problem.ALE ? :ale : :eulerian
    jac = StructuredJacobian(mode_tag, A)

    dw = build_darcy_transport_workspace(def.material)
    ctx = build_residual_context([DarcyTransportOperator()])

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
            DarcyTransportOperator(),
            u,
            ctx,
            Float64(t),
            dw,
            ws.state_cache,
            ws.mesh,
            def.material,
            ws.layout,
        )
    elseif path === :ad
        _contribute_to_jacobian_ad!(
            jac,
            DarcyTransportOperator(),
            u,
            ctx,
            Float64(t),
            dw,
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

# A Darcy fixture whose pore-gas viscosity is NOT the built-in Sutherland form.
# The analytical override takes `dμ/dT` by `ForwardDiff.derivative` on
# `t -> effective_gas_viscosity(t, material)`, so it must track whatever curve a
# material carries. With only Sutherland fixtures, a hand-coded Sutherland
# derivative would pass; this one would not. Coefficients are the Cantera
# wood-pyrolyzate blend fit (ascending, `evalpoly` order), which is also a
# stiffer test than Sutherland: a quartic has a sign-varying second derivative.
function _build_eulerian_darcy_polyvisc()
    virgin = Pyrolysis.Materials.SolidComponent(
        :wood,
        ρ = 500.0,
        c = 1500.0,
        k = 0.15,
        ε = 0.9,
        α = 1000.0,
        κ = 1e-14,
        λ = 1e-5,
    )
    gas = Pyrolysis.Materials.GaseousComponent(
        :pyrgas,
        M = 0.06,
        c = 1800.0,
        k = 0.03,
        λ = 5e-5,
    )
    rxn = Pyrolysis.Materials.Reaction(:pyrolysis, 1 => 2, A = 1e10, E = 150e3, h = -300e3)
    mat = Material(
        name = :wood_polyvisc,
        components = (virgin, gas),
        reactions = (rxn,),
        gas_viscosity = Pyrolysis.Materials.PolynomialProperty((
            -2.590750377731154e-7,
            3.664256350094585e-8,
            -6.54432524120188e-14,
            -6.352113427713809e-15,
            2.3401486488384493e-18,
        )),
    )
    mesh = create_uniform_mesh(0.01, PERF_N_CELLS, Val(2))
    bc = BoundaryConditionSet(thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()))
    problem =
        PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
    return _finalize(problem; use_ale = false)
end

@testset "DarcyTransportOperator analytical override == AD-default" begin

    # Wood solid ρ_intrinsic = 500, so all ξ_solid values stay below 500 to
    # keep porosity above the floor.
    probe_states = [
        ("default IC", (n, NC, u) -> u),
        ("uniform 600K", (n, NC, u) -> _hot_uniform_state!(u, n, NC, 600.0, (300.0, 1.5))),
        (
            "uniform 1100K",
            (n, NC, u) -> _hot_uniform_state!(u, n, NC, 1100.0, (200.0, 0.5)),
        ),
        (
            "ramp 600→1100",
            (n, NC, u) -> _ramp_state!(u, n, NC, 600.0, 1100.0, (300.0, 1.5)),
        ),
        (
            "ramp 300→1500",
            (n, NC, u) -> _ramp_state!(u, n, NC, 300.0, 1500.0, (250.0, 1.0)),
        ),
    ]

    configs = [
        ("Eulerian Wood Darcy", build_baseline_eulerian_darcy),
        ("ALE Beer-Lambert Darcy", build_ale_beer_lambert_darcy),
        # Non-Sutherland μ(T): gates the mobility chain rule against the
        # MATERIAL's curve rather than the built-in correlation.
        ("Eulerian Wood Darcy, polynomial μ(T)", _build_eulerian_darcy_polyvisc),
    ]

    for (cfg_name, builder) in configs
        def, ws, u_base, _ = builder()
        n_cells = ws.layout.n_cells
        NC = length(def.material.components)

        @testset "$cfg_name" begin
            for (state_name, set_state!) in probe_states
                u = copy(u_base)
                set_state!(n_cells, NC, u)
                J_ana = _assemble_darcy(def, ws, u, 0.0, :analytical)
                J_ad = _assemble_darcy(def, ws, u, 0.0, :ad)
                diff_abs = maximum(abs, J_ana.nzval .- J_ad.nzval)
                scale = max(maximum(abs, J_ad.nzval), 1e-30)
                rel = diff_abs / scale
                @test rel ≤ 1e-12 ||
                      (@info "DT ana vs ad" cfg_name state_name diff_abs scale rel; false)
            end
        end
    end
end
