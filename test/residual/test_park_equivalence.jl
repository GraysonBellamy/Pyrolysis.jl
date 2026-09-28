"""
Parked-state equivalence: the residual of the live cells of a state with
parked ghosts equals, to roundoff, the residual of the same physical
configuration built as an n_live-cell problem. Ghost rows are exactly zero
(T, ξ, χ) or carry the surface velocity (ALE ghost nodes). Pins every kernel's
live-range discipline at once, for every simulation mode.
"""

using Test
using Pyrolysis
using Pyrolysis.Materials: BEER_LAMBERT, NO_RADIATION, n_components
using Pyrolysis.Residual: residual!
using Pyrolysis.Solver: DepletionCallbackState, park!
using Pyrolysis.Problem:
    update_mesh_from_state!, get_node_positions, z_offset, chi_range, mesh_mode, ALE

include(joinpath(@__DIR__, "..", "helpers.jl"))

const NPE = 8

function _pmma_like()
    virgin = SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86)
    gas = GaseousComponent(:MMA, M = 0.100, c = 1500, k = 0.02, λ = 1e-5)
    rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
    return Material(name = :PMMA, components = (virgin, gas), reactions = (rxn,))
end

function _wood_like()
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

# (material, bc_set, kwargs for initialize_state)
const FIXTURES = Dict(
    :eulerian_fickian =>
        () -> (
            _pmma_like(),
            BoundaryConditionSet(
                thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()),
            ),
            (; use_ale = false),
        ),
    :ale_fickian =>
        () -> (
            _pmma_like(),
            BoundaryConditionSet(
                thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()),
            ),
            (; use_ale = true),
        ),
    :eulerian_darcy =>
        () -> (
            _wood_like(),
            BoundaryConditionSet(
                thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()),
            ),
            (; use_ale = false),
        ),
    :ale_darcy_beer_lambert =>
        () -> (
            _wood_like(),
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

function _build(material, bc_set, z_nodes, kw)
    NC = n_components(material)
    mesh = create_mesh_from_nodes(z_nodes, Val(NC))
    problem = PyrolysisProblem(
        mesh = mesh,
        material = material,
        bc_set = bc_set,
        tspan = (0.0, 1.0),
    )
    def, ws, u = initialize_state(problem; kw...)
    return def, ws, u
end

# A state with structure in every block so no kernel is trivially zero.
function _perturb!(u, layout, NC)
    n = layout.n_cells
    for i = 1:n
        u[i] = 350.0 + 60.0 * i + 5.0 * sin(0.7 * i)
        u[n+i] = u[n+i] * (1.0 - 0.08 * i)          # partial conversion of component 1
        for j = 2:NC
            u[n+(j-1)*n+i] = 0.3 * i + 0.1 * j       # some gas everywhere
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

function _blocks(du, layout, n_live, ale)
    n = layout.n_cells
    NC = n_components(layout)
    T = du[1:n_live]
    ξ = [du[n+(j-1)*n+i] for j = 1:NC, i = 1:n_live]
    z = ale ? du[z_offset(layout) .+ (1:(n_live+1))] : Float64[]
    rng = chi_range(layout)
    χ = isempty(rng) ? Float64[] : du[first(rng) .+ (0:(n_live-1))]
    return T, ξ, z, χ
end

_close(a, b) =
    isempty(a) ? isempty(b) :
    isapprox(a, b; rtol = 1e-11, atol = 1e-11 * max(maximum(abs, b), 1.0))

function _equivalence(name, parks)
    material, bc_set, kw = FIXTURES[name]()
    NC = n_components(material)
    z0 = collect(range(0.0, 0.01; length = NPE + 1))
    def, ws, u = _build(material, bc_set, z0, kw)
    layout = ws.layout
    ale = mesh_mode(layout.mode) isa ALE
    _perturb!(u, layout, NC)
    update_mesh_from_state!(ws.mesh, u, layout)
    cb = DepletionCallbackState(layout, ws.mesh, material, u)
    for p in parks
        park!(u, layout, ws.mesh, material, cb, p, 0.0)
        ws.top_face_ids[1] = ws.mesh.boundary_faces[1][1]
        ws.top_cell_ids[1] = ws.mesh.n_active_cells
    end
    n_live = ws.mesh.n_active_cells
    n = layout.n_cells

    du_p = similar(u)
    residual!(du_p, u, def, 0.0, ws)

    # Compacted reference on the live geometry with the live state.
    z_live = ale ? collect(get_node_positions(u, layout)[1:(n_live+1)]) : z0[1:(n_live+1)]
    def_c, ws_c, u_c = _build(material, bc_set, z_live, kw)
    layout_c = ws_c.layout
    for i = 1:n_live
        u_c[i] = u[i]
        for j = 1:NC
            u_c[n_live+(j-1)*n_live+i] = u[n+(j-1)*n+i]
        end
    end
    rng, rng_c = chi_range(layout), chi_range(layout_c)
    if !isempty(rng)
        for i = 1:n_live
            u_c[first(rng_c)+i-1] = u[first(rng)+i-1]
        end
        ws_c.mesh.m_dry_solid_initial .= ws.mesh.m_dry_solid_initial[1:n_live]
    end
    du_c = similar(u_c)
    residual!(du_c, u_c, def_c, 0.0, ws_c)

    Tp, ξp, zp, χp = _blocks(du_p, layout, n_live, ale)
    Tc, ξc, zc, χc = _blocks(du_c, layout_c, n_live, ale)
    @test _close(Tp, Tc)
    @test _close(ξp, ξc)
    @test _close(zp, zc)
    @test _close(χp, χc)

    # Ghost rows: zero for T/ξ/χ, surface velocity for ALE ghost nodes.
    for i = (n_live+1):n
        @test du_p[i] == 0.0
        for j = 1:NC
            @test du_p[n+(j-1)*n+i] == 0.0
        end
        isempty(rng) || @test du_p[first(rng)+i-1] == 0.0
    end
    if ale
        zoff = z_offset(layout)
        for k = (n_live+2):layout.n_nodes
            @test du_p[zoff+k] == du_p[zoff+n_live+1]
        end
    end
    return nothing
end

@testset "parked residual ≡ compacted residual" begin
    @testset "$name" for (name, parks) in (
        (:eulerian_fickian, (NPE,)),
        (:eulerian_darcy, (NPE, NPE - 1)),
        (:ale_fickian, (NPE,)),
        (:ale_fickian, (5,)),
        (:ale_fickian, (NPE, 4, 1)),
        (:ale_darcy_beer_lambert, (NPE, 6)),
        (:ale_with_chi, (NPE, 3)),
    )
        _equivalence(name, parks)
    end
end
