"""
Jacobian-mirror test for the Kirchhoff ε_eff scaling of the Beer–Lambert
injection.

The residual scales the Beer–Lambert injection by the live surface
cell's ε_eff(T_n, ξ_·,n) on the Kirchhoff default path. That makes every
energy row depend on the surface cell's state through the injected
intensity — two Jacobian channels the structured backend must carry:

- T_row(i) × T_col(n), i < n: `∂a_eff/∂T_n` (active when component
  densities are T-dependent and component emissivities differ, or when a
  Kirchhoff component's emissivity is itself a function of temperature).
  The sparsity pattern carries these slots.
- T_row(i) × ξ_col(j, n): the `∂a_eff/∂ξ_{j,n}` chain folded into the
  existing suffix couplings as `v[n] -= (∂a_eff/∂x_n)/a_eff`.

The surface transmissivity is the mixture surface absorptivity
`a_eff = Σ v_j a_j` with `a_j = ε_j(T)` for components without an explicit
`absorptivity` — identical to ε_eff on the pure Kirchhoff path, distinct
when a component overrides it.

Verified structured-vs-DenseAD with the strict per-entry criterion at
small n (the regime where a single global normalization would hide
missing blocks), for the Kirchhoff, T-dependent-ε, component-absorptivity, and
explicit-absorptivity conventions, on Eulerian and ALE meshes; the
T-dependent-ε surface column is additionally checked against central
finite differences of the residual.
"""

using Test
using Pyrolysis
using Pyrolysis.Materials: BEER_LAMBERT
using Pyrolysis.Residual: residual!
using Pyrolysis.Jacobian:
    verify_jacobian_at, JacobianSpec, DenseAD_ForwardDiff, Structured, DirectSolve

# Component emissivities differ AND (by default) the virgin density is
# T-dependent, so a_eff = Σ v_j(ξ, T)·a_j(T) has non-zero ∂/∂T and ∂/∂ξ at
# mixed composition. `ε_virgin` may be any property-function form;
# `T_dependent_density = false` leaves ∂a_eff/∂T_n to the component
# emissivity slope alone; `absorptivity_char` decouples a_eff from ε_eff.
function _build_beer_lambert_case(;
    explicit_α::Union{Nothing,Float64} = nothing,
    use_ale::Bool = false,
    ε_virgin = 0.95,
    T_dependent_density::Bool = true,
    absorptivity_char::Union{Nothing,Float64} = nothing,
    α_virgin::Float64 = 2.0,
    α_char::Float64 = 50.0,
)
    virgin = SolidComponent(
        :virgin,
        ρ = T_dependent_density ? (T -> 600.0 - 0.05 * (T - 300.0)) : 600.0,
        c = 1500.0,
        k = 0.15,
        ε = ε_virgin,
        α = α_virgin,
        λ = 1e-5,
    )
    char = SolidComponent(
        :char,
        ρ = 150.0,
        c = 1700.0,
        k = 0.08,
        ε = 0.85,
        absorptivity = absorptivity_char,
        α = α_char,
        λ = 2e-5,
    )
    gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5)
    rxn = Reaction(:pyro, 1 => (2, 3), yields = (0.3, 0.7), A = 1e10, E = 150e3, h = -300e3)
    mat = Material(name = :wood, components = (virgin, char, gas), reactions = (rxn,))
    mesh = create_uniform_mesh(0.01, 6, Val(3))
    rad_bc =
        explicit_α === nothing ? RadiativeFluxBC(flux = 50e3) :
        RadiativeFluxBC(flux = 50e3, absorptivity = explicit_α)
    bc = BoundaryConditionSet(
        thermal = (
            top = rad_bc +
                  RadiativeBC(T_inf = 300.0) +
                  ConvectiveBC(h = 10.0, T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
    )
    problem =
        PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
    I_eff = explicit_α === nothing ? 50e3 : explicit_α * 50e3
    def, ws, u = initialize_state(
        problem;
        use_ale = use_ale,
        radiation_model = BEER_LAMBERT,
        radiation_intensity = I_eff,
        radiation_kirchhoff = explicit_α === nothing,
    )
    n = ws.layout.n_cells
    # Hot graded state, mixed composition (ε_eff away from its clamp and
    # strictly between the component values).
    for i = 1:n
        u[i] = 600.0 + 400.0 * (i - 1) / (n - 1)
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

# Central-difference column of the residual w.r.t. u[col].
function _fd_column(def, ws, u, col; h = 1e-2)
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

@testset "Beer-Lambert a_eff chain in the structured Jacobian ($(case))" for (
    case,
    kwargs,
) in (
    ("Eulerian, Kirchhoff", (;)),
    ("Eulerian, explicit α", (explicit_α = 0.8,)),
    ("ALE, Kirchhoff", (use_ale = true,)),
    (
        "Eulerian, Kirchhoff, ε(T) only",
        (ε_virgin = (0.95, -1.5e-4), T_dependent_density = false),
    ),
    ("Eulerian, component absorptivity", (absorptivity_char = 0.6,)),
    (
        "Eulerian, component absorptivity + ε(T)",
        (absorptivity_char = 0.6, ε_virgin = (0.95, -1.5e-4), T_dependent_density = false),
    ),
    (
        "ALE, component absorptivity + ε(T)",
        (absorptivity_char = 0.6, ε_virgin = (0.95, -1.5e-4), use_ale = true),
    ),
)
    def, ws, u = _build_beer_lambert_case(; kwargs...)
    n = ws.layout.n_cells

    report = _strict_report(def, ws, u)
    @test report.mode == :strict
    @test report.passed

    J_s = Matrix(report.J_active)
    J_r = report.J_ref

    if get(kwargs, :explicit_α, nothing) === nothing
        # Kirchhoff path: the ∂ε_eff/∂T_n surface column must be load-bearing
        # in the reference (T-dependent density + differing component ε) and
        # reproduced per entry by the structured backend, so the sparsity
        # pattern must carry these slots.
        col_T = J_r[1:(n-1), n]
        @test maximum(abs, col_T) > 0
        @test maximum(abs, J_s[1:(n-1), n] .- col_T) <=
              1e-9 * max(1.0, maximum(abs, J_r[1:(n-1), :]))
    else
        # Explicit-α path: no ε_eff chain — interior rows must NOT couple to
        # the surface temperature through radiation (conduction only reaches
        # neighbors, so test a row at least two cells below the surface).
        @test abs(J_r[1, n]) < 1e-12 * max(1.0, maximum(abs, J_r[1, :]))
    end

    # ξ-channel surface column (attenuation + ε_eff chain): per-entry match.
    for j = 1:3
        ξn_col = n + (j - 1) * n + n
        col = J_r[1:(n-1), ξn_col]
        @test maximum(abs, J_s[1:(n-1), ξn_col] .- col) <=
              1e-9 * max(1.0, maximum(abs, col), maximum(abs, J_r[1:(n-1), :]))
    end
end

@testset "surface T column comes from the component ε(T) slope" begin
    # Constant densities: the only ∂a_eff/∂T_n is the Kirchhoff virgin
    # component's dε/dT. Rows i ≤ n − 2 couple to T_n through the injected
    # a_eff·I alone (no conduction neighbour, reactions are cell-local), so
    # the finite-difference column is the pure radiation chain — linear in
    # T_n here, so the central difference is exact up to roundoff. The
    # default fixture is optically thick (τ ≈ 9 per cell) and the chain
    # vanishes below finite-difference resolution two cells down; use a thin
    # medium so every row sees the surface.
    thin = (α_virgin = 0.2, α_char = 2.0, T_dependent_density = false)
    def, ws, u = _build_beer_lambert_case(; ε_virgin = (0.95, -1.5e-4), thin...)
    n = ws.layout.n_cells
    report = _strict_report(def, ws, u)
    @test report.passed
    J_s = Matrix(report.J_active)
    fd = _fd_column(def, ws, u, n; h = 1.0)
    rows = 1:(n-2)
    @test maximum(abs, fd[rows]) > 0
    @test maximum(abs, J_s[rows, n] .- fd[rows]) <= 1e-6 * maximum(abs, fd[rows])
    @test maximum(abs, report.J_ref[rows, n] .- fd[rows]) <= 1e-6 * maximum(abs, fd[rows])

    # Constant ε and ρ: no surface-temperature chain at all.
    def0, ws0, u0 = _build_beer_lambert_case(; ε_virgin = 0.95, thin...)
    report0 = _strict_report(def0, ws0, u0)
    @test report0.passed
    fd0 = _fd_column(def0, ws0, u0, n; h = 1.0)
    scale0 = maximum(abs, report0.J_ref[rows, :])
    @test maximum(abs, fd0[rows]) <= 1e-9 * scale0
    @test maximum(abs, report0.J_ref[rows, n]) <= 1e-12 * scale0
end

@testset "Jacobians read the incident flux set on the workspace ($(case))" for (
    case,
    use_ale,
) in (
    ("Eulerian", false),
    ("ALE", true),
)
    # The specs are built at the construction flux; the flux then changes, as
    # `solve` changes it between evaluations for a time-varying flux. Both the
    # structured backend (geometry block included) and the dense reference
    # must see the new value.
    def, ws, u = _build_beer_lambert_case(; use_ale = use_ale)
    ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
    spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
    at_50 = verify_jacobian_at(spec, def, ws, u, 0.0; reference_spec = ref_spec)
    ws.I_radiation = 80e3
    at_80 = verify_jacobian_at(spec, def, ws, u, 0.0; reference_spec = ref_spec)
    fresh = verify_jacobian_at(
        JacobianSpec(Structured(strategy = DirectSolve()), def, ws),
        def,
        ws,
        u,
        0.0;
        reference_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws),
    )
    @test at_80.passed
    @test !isapprox(at_80.J_ref, at_50.J_ref; rtol = 1e-6)
    @test isapprox(at_80.J_ref, fresh.J_ref; rtol = 1e-12)
    @test isapprox(Matrix(at_80.J_active), Matrix(fresh.J_active); rtol = 1e-12)
end
