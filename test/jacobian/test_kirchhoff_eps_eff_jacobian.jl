"""
Jacobian-mirror regression test for the Kirchhoff ε_eff scaling of the
Beer–Lambert injection.

The residual scales the Beer–Lambert injection by the live surface
cell's ε_eff(T_n, ξ_·,n) on the Kirchhoff default path. That makes every
energy row depend on the surface cell's state through the injected
intensity — two Jacobian channels the structured backend must carry:

- T_row(i) × T_col(n), i < n: `∂ε_eff/∂T_n` (active when component
  densities are T-dependent and component emissivities differ). These
  slots did not exist in the sparsity pattern pre-fix.
- T_row(i) × ξ_col(j, n): the `∂ε_eff/∂ξ_{j,n}` chain folded into the
  existing suffix couplings as `v[n] -= (∂ε_eff/∂x_n)/ε_eff`.

Verified structured-vs-DenseAD with the strict per-entry criterion at
small n (the regime where the legacy global normalization hid missing
blocks), for both the Kirchhoff and explicit-absorptivity conventions, on
Eulerian and ALE meshes.
"""

using Test
using Pyrolysis
using Pyrolysis.Materials: BEER_LAMBERT
using Pyrolysis.Jacobian:
    verify_jacobian_at, JacobianSpec, DenseAD_ForwardDiff, Structured, DirectSolve

# Component emissivities differ AND the virgin density is T-dependent, so
# ε_eff = Σ v_j(ξ, T)·ε_j has non-zero ∂/∂T and ∂/∂ξ at mixed composition.
function _build_beer_lambert_case(;
    explicit_α::Union{Nothing,Float64} = nothing,
    use_ale::Bool = false,
)
    virgin = SolidComponent(
        :virgin,
        ρ = T -> 600.0 - 0.05 * (T - 300.0),
        c = 1500.0,
        k = 0.15,
        ε = 0.95,
        α = 2.0,
        λ = 1e-5,
    )
    char =
        SolidComponent(:char, ρ = 150.0, c = 1700.0, k = 0.08, ε = 0.85, α = 50.0, λ = 2e-5)
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

@testset "Beer-Lambert ε_eff chain in the structured Jacobian ($(case))" for (
    case,
    kwargs,
) in (
    ("Eulerian, Kirchhoff", (;)),
    ("Eulerian, explicit α", (explicit_α = 0.8,)),
    ("ALE, Kirchhoff", (use_ale = true,)),
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
        # reproduced per entry by the structured backend. Pre-fix these slots
        # did not exist in the sparsity pattern at all.
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
