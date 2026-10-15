"""
Surface-radiation Jacobian with a temperature-dependent component emissivity
and a per-component surface absorptivity.

Under SURFACE_ABSORPTION the surface energy balance re-radiates at
ε_eff(T_n, ξ_·,n) and absorbs the incident flux at a_eff(T_n, ξ_·,n). With
ε_j = ε_j(T) on a Kirchhoff component both carry a ∂/∂T_n from the
emissivity slope alone (constant densities here). The surface Newton is
Dual-aware, so the chain enters the structured Jacobian through
`effective_emissivity` / `effective_absorptivity` evaluated on the seeded
(T, ξ) — no operator-level bookkeeping. Checked three ways:

1. structured vs DenseAD, strict per-entry;
2. the surface-temperature column of both against central finite
   differences of `residual!`;
3. load-bearing: a material with the *same* ε value frozen as a constant
   reproduces the residual at this state but not its Jacobian.
"""

using Test
using Pyrolysis
using Pyrolysis.Materials: SURFACE_ABSORPTION
using Pyrolysis.Residual: residual!
using Pyrolysis.Jacobian:
    verify_jacobian_at, JacobianSpec, DenseAD_ForwardDiff, Structured, DirectSolve

const _SR_SLOPE = -1.5e-4

function _build_surface_case(; ε_virgin = (0.95, _SR_SLOPE), use_ale::Bool = false)
    virgin =
        SolidComponent(:virgin, ρ = 600.0, c = 1500.0, k = 0.15, ε = ε_virgin, λ = 1e-5)
    # The char keeps its own heater absorptivity, so a_eff ≠ ε_eff while
    # the virgin component's Kirchhoff absorptivity still tracks ε(T).
    char = SolidComponent(
        :char,
        ρ = 150.0,
        c = 1700.0,
        k = 0.08,
        ε = 0.85,
        absorptivity = 0.6,
        λ = 2e-5,
    )
    gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5)
    rxn = Reaction(:pyro, 1 => (2, 3), yields = (0.3, 0.7), A = 1e10, E = 150e3, h = -300e3)
    mat = Material(name = :wood, components = (virgin, char, gas), reactions = (rxn,))
    mesh = create_uniform_mesh(0.01, 6, Val(3))
    bc = BoundaryConditionSet(
        thermal = (
            top = RadiativeFluxBC(flux = 50e3) +
                  RadiativeBC(T_inf = 300.0) +
                  ConvectiveBC(h = 10.0, T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
    )
    problem =
        PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
    def, ws, u =
        initialize_state(problem; use_ale = use_ale, radiation_model = SURFACE_ABSORPTION)
    n = ws.layout.n_cells
    for i = 1:n
        u[i] = 500.0 + 250.0 * (i - 1) / (n - 1)
        u[n+i] = 300.0
        u[2n+i] = 100.0
        u[3n+i] = 5.0
    end
    return def, ws, u
end

function _strict_report(def, ws, u)
    ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
    spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
    return verify_jacobian_at(spec, def, ws, u, 0.0; reference_spec = ref_spec)
end

function _fd_column(def, ws, u, col; h = 5e-2)
    up = copy(u)
    um = copy(u)
    up[col] += h
    um[col] -= h
    dp = similar(u)
    dm = similar(u)
    residual!(dp, up, def, 0.0, ws)
    residual!(dm, um, def, 0.0, ws)
    return (dp .- dm) ./ (2h)
end

@testset "surface radiation Jacobian with ε(T) and component absorptivity ($(case))" for (
    case,
    use_ale,
) in (
    ("Eulerian", false),
    ("ALE", true),
)
    def, ws, u = _build_surface_case(; use_ale = use_ale)
    n = ws.layout.n_cells
    report = _strict_report(def, ws, u)
    @test report.mode == :strict
    @test report.passed

    # Surface-temperature column: the boundary-cell energy row and its
    # conduction neighbour, against finite differences of the residual.
    fd = _fd_column(def, ws, u, n)
    J_s = Matrix(report.J_active)
    for row in (n, n - 1)
        @test J_s[row, n] ≈ fd[row] rtol = 1e-5
        @test report.J_ref[row, n] ≈ fd[row] rtol = 1e-5
    end
end

@testset "the ε(T) slope is load-bearing in the surface energy row" begin
    def, ws, u = _build_surface_case()
    n = ws.layout.n_cells
    # Same emissivity value at the surface cell, frozen as a constant.
    ε_frozen = 0.95 + _SR_SLOPE * u[n]
    def0, ws0, u0 = _build_surface_case(; ε_virgin = ε_frozen)

    du = similar(u)
    du0 = similar(u)
    residual!(du, u, def, 0.0, ws)
    residual!(du0, u0, def0, 0.0, ws0)
    # Only the surface cell's ε enters under SURFACE_ABSORPTION, so the two
    # materials agree on the residual at this state…
    @test du ≈ du0 rtol = 1e-12

    # …but not on ∂/∂T_n: the live material carries F·σ·(T_∞⁴ − T_s⁴)·dε/dT
    # (emission) and q_inc·da/dT (absorption, Kirchhoff virgin component).
    # The slope term is a small fraction of the diagonal (which the local
    # reaction heat dominates), so the load-bearing check is relative to the
    # finite-difference tolerance the live Jacobian must meet, not to the
    # diagonal itself.
    fd = _fd_column(def, ws, u, n)
    fd_tol = 1e-5 * abs(fd[n])
    J_live = Matrix(_strict_report(def, ws, u).J_active)
    J_frozen = Matrix(_strict_report(def0, ws0, u0).J_active)
    @test abs(J_live[n, n] - fd[n]) <= fd_tol
    @test abs(J_frozen[n, n] - fd[n]) > 20 * fd_tol
end
