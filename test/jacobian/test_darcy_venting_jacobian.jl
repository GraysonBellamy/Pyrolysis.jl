"""
Structured-vs-dense Jacobian agreement for the boundary Darcy venting
problem: the intrinsic-density φ_face chain, solid-basis κ mixing, the
zero-mobility rule at impermeable faces, and the boundary venting term.

Per-entry tolerance at small n — deliberately NOT a global-max
normalization, which hides missing O(small) blocks behind the largest
Arrhenius entries. The state exercises every part of the venting chain at
once:

- a char/virgin front (φ and κ vary cell-to-cell → the φ_face chain and the
  solid-basis κ mixing feed the analytical operator),
- a face with one impermeable side (the zero-mobility rule),
- a top `DiffusiveMassBC + DirichletPressureBC` with overpressure (the
  boundary venting term in `BCMassOperator`).

Only the ξ-rows are compared here: the boundary cells' T-row coupling to
the boundary mass flux is owned by `BCMassOperator` and is covered by the
boundary-flux energy-row tests in `test_jacobian_coupling_regressions.jl`.
"""

using Test
using SparseArrays
using Pyrolysis
using Pyrolysis.Materials
using Pyrolysis.BoundaryConditions
using Pyrolysis.Geometry
using Pyrolysis.Jacobian: JacobianSpec, DenseAD_ForwardDiff, Structured

@testset "Darcy venting: structured ξ-rows == DenseAD per-entry" begin
    mat = Material(
        name = :venting_jac,
        components = (
            SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2),       # κ = 0
            SolidComponent(:char, ρ = 500.0, c = 1000.0, k = 0.1, κ = 1e-12),
            GaseousComponent(:gas, M = 0.028, c = 1000.0, k = 0.03, λ = 1e-6),
        ),
        reactions = (),
    )
    n = 6
    mesh = create_uniform_mesh(0.01, n, Val(3))
    bc = BoundaryConditionSet(
        thermal = (top = HeatFluxBC(0.0), bottom = AdiabaticBC()),
        mass = (
            top = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0),
            bottom = ImpermeableBC(),
        ),
    )
    problem = PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 1.0))
    def, ws, u = initialize_state(problem; use_ale = false)

    # Char fraction ramps toward the surface: φ, κ, and P vary cell-to-cell,
    # the first faces have an impermeable side, the top cell has a large
    # overpressure against the Dirichlet BC.
    for i = 1:n
        α = (i - 1) / (n - 1)
        u[i] = 400.0 + 500.0 * α
        u[n+i] = 800.0 * (1 - α)          # virgin
        u[2n+i] = 300.0 * α               # char
        u[3n+i] = 0.5 + 1.5 * α           # gas
    end

    ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
    spec = JacobianSpec(Structured(), def, ws)

    J_ref = zeros(length(u), length(u))
    J_str = Pyrolysis.Jacobian._alloc_target(spec, length(u))
    Pyrolysis.Jacobian.compute_jacobian!(J_ref, ref_spec, u, def, 0.0, ws)
    Pyrolysis.Jacobian.compute_jacobian!(J_str, spec, u, def, 0.0, ws)

    Jr = Matrix(J_ref)
    Js = Matrix(J_str)

    n_bad = 0
    for r = (n+1):(4n), c = 1:size(Jr, 2)
        dev = abs(Js[r, c] - Jr[r, c])
        if dev > 1e-6 + 1e-6 * abs(Jr[r, c])
            n_bad += 1
        end
    end
    @test n_bad == 0

    # The venting term is actually present (not trivially zero): the top
    # cell's gas-row diagonal must include −(v_b/φ + film)·A/V on top of the
    # interior-face couplings, and the dense reference must see it too.
    gas_row_top = 3n + n
    @test Jr[gas_row_top, gas_row_top] < 0
    @test abs(Js[gas_row_top, gas_row_top] - Jr[gas_row_top, gas_row_top]) <=
          1e-6 * abs(Jr[gas_row_top, gas_row_top])
end
