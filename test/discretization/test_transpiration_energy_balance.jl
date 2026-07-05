"""
Surface-transpiration energy-balance regression tests.

Guards against the gas-convective-enthalpy DOUBLE-COUNT documented in
`docs/surface_transpiration_energy_review.md`. The volumetric advection source
`S_conv = -Σ_g cp_g J_g ∂T/∂x` is always active; the surface transpiration term
`Σ_g ṁ_g·Δh_g(T_s,T_∞)`, folded into the boundary energy balance, overlaps it by
`∫S_conv`. Applying BOTH (the old default) over-applies the gas sink — the
supported convention (ThermaKin Eq.17) is S_conv only, with transpiration OFF.

These tests:
  1. Assert the construction-time guard rejects the double-count combination.
  2. Build an INDEPENDENT global energy ledger from one `compute_rhs!` evaluation
     and prove the both-on configuration removes an extra `Σ_g ṁ_g·Δh_g·A` of
     energy on top of the (unchanged) volumetric `S_conv` — exactly the spurious
     sink — while the default (transpiration-off) does not.
"""

using Test
using Pyrolysis
using Pyrolysis.Materials
using Pyrolysis.BoundaryConditions
using Pyrolysis.Geometry
using Pyrolysis.Properties
using Pyrolysis.Physics
using Pyrolysis.Discretization
using Pyrolysis.Problem: StateCache

include(joinpath(@__DIR__, "..", "helpers.jl"))

@testset "Surface transpiration energy balance" begin

    # ------------------------------------------------------------------
    # 1. Construction-time guard rejects the double-count combination.
    # ------------------------------------------------------------------
    @testset "Guard rejects use_transpiration_bc=true" begin
        problem = fixture_problem()

        # Default (transpiration off) constructs fine.
        @test (to_problem_def(problem); true)
        @test (to_problem_def(problem; use_transpiration_bc = false); true)

        # Enabling transpiration alongside the always-on S_conv is rejected.
        @test_throws ArgumentError to_problem_def(problem; use_transpiration_bc = true)
        @test_throws ArgumentError to_problem_def(
            problem;
            energy_form = :with_generation_sink,
            use_transpiration_bc = true,
        )

        # `solve` fails fast on the same combination.
        @test_throws ArgumentError solve(problem; use_transpiration_bc = true)
    end

    # ------------------------------------------------------------------
    # 2. Independent energy ledger: both-on removes an extra surface
    #    outflow on top of the unchanged volumetric S_conv (the double-count).
    # ------------------------------------------------------------------
    @testset "Double-count magnitude via independent ledger" begin
        # Constant-cp gas so the surface enthalpy is exactly cp_g·(T_s−T_∞)
        # (no midpoint-rule ambiguity in the analytic cross-check).
        cp_gas = 1500.0
        solid = SolidComponent(:solid, ρ = 1000.0, c = 1000.0, k = 0.2, ε = 0.9)
        gas = GaseousComponent(:vol, M = 0.030, c = cp_gas, k = 0.03, λ = 1e-5)
        mat = Material(name = :transp_test, components = (solid, gas))

        L = 0.01
        n = 20
        mesh = create_uniform_mesh(L, n, Val(2))

        # Hot surface (top), gas present throughout so it flows out the top and
        # S_conv is non-trivial. T_s will solve well above the 298 K ambient
        # datum, so the outflow enthalpy is a genuine sink.
        for i = 1:n
            z = mesh.cell_centers[i]
            T = 400.0 + 600.0 * z / L           # 400 K (back) → 1000 K (surface)
            vol = mesh.cell_states[i].volume
            mesh.cell_states[i] = Pyrolysis.Geometry.CellState{2}(T, (900.0, 50.0), vol)
        end

        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()),
            mass = (top = DiffusiveMassBC(0.05), bottom = ImpermeableBC()),
        )

        top_face_ids = mesh.boundary_faces[1]
        bottom_face_ids = mesh.boundary_faces[2]

        prop_cache = PropertyCache(mesh)
        update_properties!(prop_cache, StateCache(mesh, mat), mesh, mat)
        rhs_cache = RHSCache(n, n_faces(mesh), 2, 0)

        # Helper: one RHS evaluation; returns the global energy-rate ledger and
        # the boundary diagnostics. The energy rate into the slab is
        # Σ_i ρcp_i · dT_i · V_i = ∫(-∇·q + S_conv + S_gen + Q_rad + Q_rxn) dV.
        function ledger(use_transp::Bool)
            state = StateCache(mesh, mat)
            dT = zeros(n)
            dξ = zeros(2, n)
            compute_rhs!(
                dT,
                dξ,
                state,
                mesh,
                prop_cache,
                rhs_cache,
                mat,
                bc_set,
                0.0,
                top_face_ids,
                bottom_face_ids;
                use_transpiration_bc = use_transp,
            )
            E_rate = sum(prop_cache.ρcp[i] * dT[i] * state.volume[i] for i = 1:n)
            S_conv_int = sum(rhs_cache.S_conv[i] * state.volume[i] for i = 1:n)
            bstate = mesh.boundary_states[1]
            A_top = state.A_face[first(top_face_ids)]
            return (
                E_rate = E_rate,
                S_conv_int = S_conv_int,
                q_transpired = bstate.q_transpired,
                T_s = bstate.T_surface,
                A_top = A_top,
            )
        end

        off = ledger(false)   # default supported convention (S_conv only)
        both = ledger(true)   # buggy double-count (S_conv + surface outflow)

        # The volumetric advection sink is active in both runs. Since the
        # Minor-5 fix the boundary cell's S_conv reads the surface slope from
        # the thermal-BC pass (s_face = −q_total/k), so S_conv now responds to
        # the transpiration-shifted T_s — the two runs are no longer equal.
        @test abs(off.S_conv_int) > 0.0
        @test off.S_conv_int != both.S_conv_int

        # Surface is hot ⇒ outflow carries a positive enthalpy sink only in the
        # both-on run; the default reports zero transpiration.
        @test off.q_transpired == 0.0
        @test both.q_transpired > 0.0

        # Cross-check the transpiration diagnostic against the closed form
        # Σ_g ṁ_g·cp_g·(T_s−T_∞): recompute ṁ_top from the gas flux at the
        # surface face used by the boundary solve.
        ṁ_top = rhs_cache.J_gas[2, first(top_face_ids)]   # +z outward at the top
        T_inf = Float64(Pyrolysis.Materials.T_REF)
        @test ṁ_top > 0.0
        @test both.q_transpired ≈ ṁ_top * cp_gas * (both.T_s - T_inf) rtol = 1e-6

        # CORE ASSERTION: the exact flag-difference ledger. The divergence
        # term telescopes to the boundary q_total, and the surface Newton
        # gives q_cond = q_BC − q_transpire with a constant q_BC here, so
        #   E_off − E_both = q_transpired·A_top + (∫S_conv_off − ∫S_conv_both).
        # The first term is the spurious surface double-count; the second is
        # the boundary cell's S_conv responding to the transpiration-shifted
        # T_s (Minor-5 coupling). With the pre-Minor-5 flag-independent
        # S_conv this reduces to the original identity
        # E_off − E_both = q_transpired·A_top.
        @test (off.E_rate - both.E_rate) ≈
              both.q_transpired * both.A_top + (off.S_conv_int - both.S_conv_int) rtol =
            1e-8

        # Magnitude context: the spurious extra sink is the same order as the
        # advective enthalpy the slab already carries via S_conv — i.e. roughly
        # doubling the gas-energy sink (the review's factor-of-2 result). Not an
        # exact 2× for a general profile, but confirmed comparable here.
        extra_sink = both.q_transpired * both.A_top
        @test extra_sink > 0.2 * abs(off.S_conv_int)
    end
end
