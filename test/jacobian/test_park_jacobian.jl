"""
Structured Jacobian on parked states. The specs are built on the fresh
problem (fixed sparsity pattern and caches) and must stay exact after cells
are parked: `verify_jacobian_at` against the dense ForwardDiff reference,
strict for `DirectSolve`, masked for `ApproxSolve`; ghost rows are zero in
both, and the parked-state Jacobian equals the compacted problem's on the
live block.
"""

using Test
using SparseArrays
using Pyrolysis
using Pyrolysis.Materials: BEER_LAMBERT, n_components
using Pyrolysis.Jacobian:
    JacobianSpec,
    Structured,
    DirectSolve,
    ApproxSolve,
    DenseAD_ForwardDiff,
    verify_jacobian_at,
    compute_jacobian!
using Pyrolysis.Solver: DepletionCallbackState, park!
using Pyrolysis.Problem:
    update_mesh_from_state!, get_node_positions, z_offset, chi_range, mesh_mode, ALE

include(joinpath(@__DIR__, "..", "helpers.jl"))

const NPJ = 8

function _pmma_like_j()
    virgin = SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86)
    gas = GaseousComponent(:MMA, M = 0.100, c = 1500, k = 0.02, λ = 1e-5)
    rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
    return Material(name = :PMMA, components = (virgin, gas), reactions = (rxn,))
end

function _wood_like_j()
    virgin = SolidComponent(
        :wood,
        ρ = 500.0,
        c = 1500.0,
        k = 0.15,
        ε = 0.9,
        α = 1000.0,
        κ = 1e-14,
        λ = 1e-5,
    )
    gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5)
    rxn = Reaction(:pyrolysis, 1 => 2, A = 1e10, E = 150e3, h = -300e3)
    return Material(name = :wood, components = (virgin, gas), reactions = (rxn,))
end

const JFIX = Dict(
    :eulerian_fickian =>
        () -> (
            _pmma_like_j(),
            BoundaryConditionSet(
                thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()),
            ),
            (; use_ale = false),
        ),
    :ale_fickian =>
        () -> (
            _pmma_like_j(),
            BoundaryConditionSet(
                thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()),
            ),
            (; use_ale = true),
        ),
    :eulerian_darcy =>
        () -> (
            _wood_like_j(),
            BoundaryConditionSet(
                thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()),
            ),
            (; use_ale = false),
        ),
    :ale_darcy_beer_lambert =>
        () -> (
            _wood_like_j(),
            BoundaryConditionSet(
                thermal = (
                    top = ConvectiveBC(h = 10.0, T_inf = 300.0),
                    bottom = AdiabaticBC(),
                ),
            ),
            (; use_ale = true, radiation_model = BEER_LAMBERT, radiation_intensity = 50e3),
        ),
    :ale_with_chi =>
        () -> (
            fixture_material(lateral = true),
            BoundaryConditionSet(
                thermal = (top = HeatFluxBC(1e4), bottom = AdiabaticBC()),
                mass = (top = DiffusiveMassBC(0.01), bottom = ImpermeableBC()),
            ),
            (; use_ale = true),
        ),
)

function _build_j(material, bc_set, z_nodes, kw)
    NC = n_components(material)
    mesh = create_mesh_from_nodes(z_nodes, Val(NC))
    problem = PyrolysisProblem(
        mesh = mesh,
        material = material,
        bc_set = bc_set,
        tspan = (0.0, 1.0),
    )
    return initialize_state(problem; kw...)
end

function _perturb_j!(u, layout, NC)
    n = layout.n_cells
    for i = 1:n
        u[i] = 350.0 + 60.0 * i + 5.0 * sin(0.7 * i)
        u[n+i] = u[n+i] * (1.0 - 0.08 * i)
        for j = 2:NC
            u[n+(j-1)*n+i] = 0.3 * i + 0.1 * j
        end
    end
    rng = chi_range(layout)
    isempty(rng) || (
        for i = 1:n

            u[first(rng)+i-1] = 0.04 * i
        end
    )
    return u
end

# Global indices of the live block of a parked state (T, ξ, live nodes, χ).
function _live_indices(layout, n_live, ale)
    n = layout.n_cells
    NC = n_components(layout)
    idx = collect(1:n_live)
    for j = 1:NC
        append!(idx, n .+ (j - 1) * n .+ (1:n_live))
    end
    ale && append!(idx, z_offset(layout) .+ (1:(n_live+1)))
    rng = chi_range(layout)
    isempty(rng) || append!(idx, first(rng) .+ (0:(n_live-1)))
    return idx
end

function _check(name, parks)
    material, bc_set, kw = JFIX[name]()
    NC = n_components(material)
    z0 = collect(range(0.0, 0.01; length = NPJ + 1))
    def, ws, u = _build_j(material, bc_set, z0, kw)
    layout = ws.layout
    ale = mesh_mode(layout.mode) isa ALE

    # Specs built on the fresh problem, as solve() does.
    ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
    spec_direct = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
    spec_approx =
        ale ? JacobianSpec(Structured(strategy = ApproxSolve()), def, ws) : nothing

    _perturb_j!(u, layout, NC)
    update_mesh_from_state!(ws.mesh, u, layout)
    cb = DepletionCallbackState(layout, ws.mesh, material, u)
    for p in parks
        park!(u, layout, ws.mesh, material, cb, p, 0.0)
        ws.top_face_ids[1] = ws.mesh.boundary_faces[1][1]
    end
    n_live = ws.mesh.n_active_cells
    n = layout.n_cells

    report = verify_jacobian_at(
        spec_direct,
        def,
        ws,
        u,
        0.0;
        reference_spec = ref_spec,
        rtol = 1e-6,
        atol = 1e-10,
    )
    @test report.passed
    report.passed ||
        @info "DirectSolve on parked state" name parks report.max_rel_dev report.worst_index

    # Ghost T/ξ/χ rows are identically zero in the reference and the backend.
    J_ref = report.J_ref
    J_act = Matrix(report.J_active)
    ghost_rows = Int[]
    for i = (n_live+1):n
        push!(ghost_rows, i)
        for j = 1:NC
            push!(ghost_rows, n + (j - 1) * n + i)
        end
        rng = chi_range(layout)
        isempty(rng) || push!(ghost_rows, first(rng) + i - 1)
    end
    @test all(iszero, J_ref[ghost_rows, :])
    @test all(iszero, J_act[ghost_rows, :])
    # Ghost node rows copy the surface node row (dz/dt = w_top for every ghost).
    if ale
        zoff = z_offset(layout)
        for k = (n_live+2):layout.n_nodes
            @test J_ref[zoff+k, :] == J_ref[zoff+n_live+1, :]
        end
    end

    # Live block equals the compacted problem's Jacobian.
    z_live = ale ? collect(get_node_positions(u, layout)[1:(n_live+1)]) : z0[1:(n_live+1)]
    def_c, ws_c, u_c = _build_j(material, bc_set, z_live, kw)
    layout_c = ws_c.layout
    for i = 1:n_live
        u_c[i] = u[i]
        for j = 1:NC
            u_c[n_live+(j-1)*n_live+i] = u[n+(j-1)*n+i]
        end
    end
    if ale
        for k = 1:(n_live+1)
            u_c[z_offset(layout_c)+k] = u[z_offset(layout)+k]
        end
    end
    rng, rng_c = chi_range(layout), chi_range(layout_c)
    if !isempty(rng)
        for i = 1:n_live
            u_c[first(rng_c)+i-1] = u[first(rng)+i-1]
        end
        ws_c.mesh.m_dry_solid_initial .= ws.mesh.m_dry_solid_initial[1:n_live]
    end
    J_c = zeros(length(u_c), length(u_c))
    compute_jacobian!(
        J_c,
        JacobianSpec(DenseAD_ForwardDiff(), def_c, ws_c),
        u_c,
        def_c,
        0.0,
        ws_c,
    )
    live = _live_indices(layout, n_live, ale)
    J_live = J_ref[live, live]
    scale = max(maximum(abs, J_c), 1.0)
    @test maximum(abs, J_live .- J_c) <= 1e-9 * scale

    # ApproxSolve drops couplings by design, so it is not compared against
    # the dense reference on a reacting state; parking must leave what it
    # assembles unchanged: parked live block ≡ compacted problem's block.
    if spec_approx !== nothing
        J_ap = Matrix(_assemble(spec_approx, u, def, ws))
        spec_approx_c = JacobianSpec(Structured(strategy = ApproxSolve()), def_c, ws_c)
        J_ac = Matrix(_assemble(spec_approx_c, u_c, def_c, ws_c))
        scale_a = max(maximum(abs, J_ac), 1.0)
        @test maximum(abs, J_ap[live, live] .- J_ac) <= 1e-9 * scale_a
        @test all(iszero, J_ap[ghost_rows, :])
    end
    return nothing
end

function _assemble(spec, u, def, ws)
    S = spec.sparsity_phys .| spec.sparsity_geom
    J = SparseMatrixCSC{Float64,Int}(
        S.m,
        S.n,
        copy(S.colptr),
        copy(S.rowval),
        zeros(length(S.nzval)),
    )
    compute_jacobian!(J, spec, u, def, 0.0, ws)
    return J
end

@testset "structured Jacobian on parked states" begin
    @testset "$name $parks" for (name, parks) in (
        (:eulerian_fickian, (NPJ,)),
        (:eulerian_darcy, (NPJ, NPJ - 1)),
        (:ale_fickian, (NPJ,)),
        (:ale_fickian, (5,)),
        (:ale_fickian, (NPJ, 4, 1)),
        (:ale_darcy_beer_lambert, (NPJ, 6)),
        (:ale_with_chi, (NPJ, 3)),
    )
        _check(name, parks)
    end
end
