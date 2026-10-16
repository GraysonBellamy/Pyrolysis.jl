"""
Regression tests for boundary Darcy venting and the impermeable-interface
mobility rule.

- The pressure-BC boundary Darcy velocity enters the species equations as
  an advective venting flux J = v_b·ξ_int/φ_int (outflow-upwinded against a
  zero-concentration ambient). The failing regime before the fix: venting was
  independent of permeability and ∇P (film law only).
- The venting flux carries the 1/φ intrinsic-density factor.
- A char(κ>0)/virgin(κ=0) front has exactly zero face mobility — no flow
  leaks into the impermeable side.
- NeumannPressureBC's ∂P/∂n is an outward-normal gradient — the bottom
  boundary maps it to ∂P/∂z with a sign flip.
- ConvectivePressureBC cannot vent through an impermeable surface cell
  (κ_int = 0 ⇒ v_b = 0).
"""

using Test
using Pyrolysis
using Pyrolysis.Residual: residual!
using Pyrolysis.Materials
using Pyrolysis.BoundaryConditions
using Pyrolysis.Geometry
using Pyrolysis.Properties
using Pyrolysis.Physics
using Pyrolysis.Discretization: boundary_darcy_velocity, update_darcy_properties!
using Pyrolysis.Properties: PropertyCache, update_properties!

# Two-component venting fixture: hot uniform column with a top
# DiffusiveMassBC + DirichletPressureBC (the documented Darcy venting idiom).
function _venting_problem(κ_char; n = 6, T0 = 600.0, ξ_solid = 250.0, ξ_gas = 2.0)
    mat = Material(
        name = :venting,
        components = (
            SolidComponent(:char, ρ = 500.0, c = 1000.0, k = 0.1, κ = κ_char),
            GaseousComponent(:gas, M = 0.028, c = 1000.0, k = 0.03, λ = 1e-6),
        ),
        reactions = (),
    )
    mesh = create_uniform_mesh(0.01, n, Val(2))
    bc = BoundaryConditionSet(
        thermal = (top = HeatFluxBC(0.0), bottom = AdiabaticBC()),
        mass = (
            top = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0),
            bottom = ImpermeableBC(),
        ),
    )
    problem = PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 1.0))
    def, ws, u = initialize_state(problem; use_ale = false)
    for i = 1:n
        u[i] = T0
        u[n+i] = ξ_solid
        u[2n+i] = ξ_gas
    end
    return def, ws, u, similar(u), n
end

@testset "Boundary Darcy venting" begin

    @testset "venting rate depends on permeability and ∇P" begin
        rate = Dict{Float64,Float64}()
        for κ in (1e-14, 1e-12, 1e-11)
            def, ws, u, du, n = _venting_problem(κ)
            residual!(du, u, def, 0.0, ws)
            rate[κ] = du[2n+n]   # gas ξ tendency in the surface cell
        end
        # Gas leaves the surface cell, faster with higher κ — the defining
        # regression: before the fix the rate was the κ-independent film law.
        @test rate[1e-12] < 0
        @test abs(rate[1e-12]) > 10 * abs(rate[1e-14])
        @test abs(rate[1e-11]) > 5 * abs(rate[1e-12])
    end

    @testset "boundary flux = v_b·ξ/φ + film law (uniform state, exact)" begin
        κ = 1e-12
        def, ws, u, du, n = _venting_problem(κ)
        residual!(du, u, def, 0.0, ws)

        mat = def.material
        T = 600.0
        ξ = (250.0, 2.0)
        φ = porosity(ξ, T, mat)              # = 0.5 for this fixture
        @test φ ≈ 0.5 rtol = 1e-12
        P_int = ideal_gas_pressure(ξ, T, φ, mat)
        μ = SutherlandProperty()(T)
        κ_eff = effective_permeability(ξ, T, mat)
        Δz = 0.01 / n
        d = Δz / 2
        v_b = κ_eff / μ * (P_int - 101325.0) / d      # outward (+z) at the top
        @test v_b > 0                                  # overpressure vents out

        J_adv = v_b * ξ[2] / φ                         # intrinsic density
        λ = effective_gas_transfer(ξ, T, mat)
        h_m = 0.05
        J_diff = λ * h_m * ξ[2] / (λ + h_m * d)

        # Uniform state ⇒ all interior fluxes vanish; the surface-cell gas
        # tendency is exactly the boundary outflux divided by the cell height.
        expected = -(J_adv + J_diff) / Δz
        @test du[2n+n] ≈ expected rtol = 1e-10

        # The 1/φ factor is present: the advective part alone is 2× the
        # bulk-concentration flux v_b·ξ at φ = 0.5.
        @test J_adv ≈ 2.0 * v_b * ξ[2] rtol = 1e-12
    end

    @testset "no pressure BC ⇒ no advective venting (film law only)" begin
        # Same fixture but mass BC without the pressure component.
        mat = Material(
            name = :no_vent,
            components = (
                SolidComponent(:char, ρ = 500.0, c = 1000.0, k = 0.1, κ = 1e-12),
                GaseousComponent(:gas, M = 0.028, c = 1000.0, k = 0.03, λ = 1e-6),
            ),
            reactions = (),
        )
        n = 6
        mesh = create_uniform_mesh(0.01, n, Val(2))
        bc = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(0.0), bottom = AdiabaticBC()),
            mass = (top = DiffusiveMassBC(0.05), bottom = ImpermeableBC()),
        )
        problem =
            PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 1.0))
        def, ws, u = initialize_state(problem; use_ale = false)
        for i = 1:n
            u[i] = 600.0
            u[n+i] = 250.0
            u[2n+i] = 2.0
        end
        du = similar(u)
        residual!(du, u, def, 0.0, ws)

        T = 600.0
        ξ = (250.0, 2.0)
        λ = effective_gas_transfer(ξ, T, mat)
        Δz = 0.01 / n
        d = Δz / 2
        h_m = 0.05
        J_diff = λ * h_m * ξ[2] / (λ + h_m * d)
        @test du[2n+n] ≈ -J_diff / Δz rtol = 1e-10
    end
end

@testset "Impermeable-interface mobility" begin
    # char(κ>0) cell next to virgin(κ=0) cell: the shared face must carry
    # exactly zero Darcy velocity — infinite series resistance.
    mat = Material(
        name = :front,
        components = (
            SolidComponent(:virgin, ρ = 1000.0, c = 1500.0, k = 0.2),       # κ = 0
            SolidComponent(:char, ρ = 500.0, c = 1000.0, k = 0.1, κ = 1e-12),
            GaseousComponent(:gas, M = 0.028, c = 1000.0, k = 0.03, λ = 1e-6),
        ),
        reactions = (),
    )
    n = 4
    mesh = create_uniform_mesh(0.01, n, Val(3))
    bc = BoundaryConditionSet(
        thermal = (top = HeatFluxBC(0.0), bottom = AdiabaticBC()),
        mass = (top = DiffusiveMassBC(0.05), bottom = ImpermeableBC()),
    )
    problem = PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 1.0))
    def, ws, u = initialize_state(problem; use_ale = false)
    state = ws.state_cache

    # Cells 1–2 virgin (impermeable), cells 3–4 char with a huge gas
    # overpressure differential across the front.
    gas_profile = (0.1, 0.1, 5.0, 2.0)     # pressure jump across the front AND
    for i = 1:n                             # within the char region
        u[i] = 600.0
        u[n+i] = i <= 2 ? 800.0 : 0.0      # virgin
        u[2n+i] = i <= 2 ? 0.0 : 400.0     # char
        u[3n+i] = gas_profile[i]
    end
    Pyrolysis.Problem.populate!(state, u, ws.layout)
    Pyrolysis.Problem.update_state_geometry!(
        state,
        ws.mesh,
        Pyrolysis.Problem.mesh_mode(def.mode),
    )
    update_properties!(ws.prop_cache, state, ws.mesh, def.material)
    update_darcy_properties!(
        ws.prop_cache,
        state,
        ws.mesh,
        def.material,
        def.bcs,
        ws.top_face_ids,
        ws.bottom_face_ids,
    )
    darcy = ws.prop_cache.darcy

    # Face 3 separates cell 2 (virgin, κ=0) from cell 3 (char, κ>0): the
    # one-sided fallback used to assign the char side's full mobility here.
    front_face = 3
    @test ws.mesh.faces[front_face].cell_ids == (2, 3)
    @test darcy.v_gas[front_face] == 0.0
    # Char-side interior face still flows.
    @test darcy.v_gas[4] != 0.0
end

@testset "Boundary Darcy velocity kernel" begin
    P_int = 150_000.0
    κ = 1e-12
    μ = 2e-5
    d = 1e-3

    @testset "Neumann ∂P/∂n is outward-normal (bottom sign flip)" begin
        # Pressure DECREASING outward on both boundaries: ∂P/∂n < 0.
        bc = NeumannPressureBC(-1e6)
        v_top = boundary_darcy_velocity(bc, P_int, κ, μ, d, Val(:top))
        v_bot = boundary_darcy_velocity(bc, P_int, κ, μ, d, Val(:bottom))
        # Outflow at the top is +z; outflow at the bottom is −z.
        @test v_top > 0
        @test v_bot < 0
        @test v_top ≈ -v_bot rtol = 1e-12
        @test v_top ≈ κ / μ * 1e6 rtol = 1e-12
    end

    @testset "Convective BC gated on interior permeability" begin
        bc = ConvectivePressureBC(1e-6, 101325.0)
        v = boundary_darcy_velocity(bc, P_int, κ, μ, d, Val(:top))
        @test v ≈ 1e-6 * (P_int - 101325.0) rtol = 1e-12
        # Impermeable surface cell cannot vent regardless of h_P.
        @test boundary_darcy_velocity(bc, P_int, 0.0, μ, d, Val(:top)) == 0.0
        # Bottom orientation: outflow (P_int > P_ext) is −z.
        @test boundary_darcy_velocity(bc, P_int, κ, μ, d, Val(:bottom)) < 0
    end

    @testset "Dirichlet overpressure vents outward on both boundaries" begin
        bc = DirichletPressureBC(101325.0)
        v_top = boundary_darcy_velocity(bc, P_int, κ, μ, d, Val(:top))
        v_bot = boundary_darcy_velocity(bc, P_int, κ, μ, d, Val(:bottom))
        @test v_top > 0       # +z out the top
        @test v_bot < 0       # −z out the bottom
        @test v_top ≈ -v_bot rtol = 1e-12
        # Impermeable interior cell: no flow.
        @test boundary_darcy_velocity(bc, P_int, 0.0, μ, d, Val(:top)) == 0.0
    end

    @testset "Combined BC sums children; non-pressure BCs contribute zero" begin
        comb = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0)
        v_comb = boundary_darcy_velocity(comb, P_int, κ, μ, d, Val(:top))
        v_dir = boundary_darcy_velocity(
            DirichletPressureBC(101325.0),
            P_int,
            κ,
            μ,
            d,
            Val(:top),
        )
        @test v_comb ≈ v_dir rtol = 1e-14
        @test boundary_darcy_velocity(DiffusiveMassBC(0.05), P_int, κ, μ, d, Val(:top)) ==
              0.0
    end
end
