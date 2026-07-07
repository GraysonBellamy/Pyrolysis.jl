"""
Regression tests for Beer–Lambert surface absorptivity: the default
(Kirchhoff) Beer–Lambert path never applied a surface absorptivity — the
medium absorbed 100 % of the incident flux while emitting with ε_eff, so
switching SURFACE_ABSORPTION → BEER_LAMBERT changed the absorbed power by
1/ε_eff (≈ 11 % at ε = 0.9).

Post-fix convention:
- Kirchhoff path (`apply_kirchhoff_absorptivity = true`, no explicit BC
  absorptivity): the flux entering the medium is ε_eff·I, with ε_eff the
  *live* surface cell's effective emissivity — absorption and emission use
  the same coefficient (detailed balance restored).
- Explicit-absorptivity path (`false`): α(t)·I is folded into the intensity
  upstream (`extract_radiative_flux_from_bcs`) and the kernel does not
  rescale.
"""

using Test
using Pyrolysis
using Pyrolysis.Materials: NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT, P1_QUASI_STEADY
using Pyrolysis.Physics: compute_beer_lambert_source!, compute_radiation_source!
using Pyrolysis.Properties: effective_emissivity
using Pyrolysis.BoundaryConditions: extract_radiative_flux_from_bcs
using Pyrolysis.Problem: populate!

# Two solids with different ε so ε_eff is composition-dependent; non-zero
# absorption coefficients so Beer-Lambert is non-trivial.
function _radiation_fixture()
    virgin = SolidComponent(
        :virgin,
        ρ = 600.0,
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
    mesh = create_uniform_mesh(0.01, 8, Val(3))
    bc = BoundaryConditionSet(
        thermal = (
            top = RadiativeFluxBC(flux = 50e3) + RadiativeBC(T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
    )
    problem =
        PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
    def, ws, u = initialize_state(problem; radiation_model = BEER_LAMBERT)
    n = ws.layout.n_cells
    # Graded hot state with mixed composition (ε_eff strictly between the
    # component values, away from the clamp).
    for i = 1:n
        u[i] = 600.0 + 300.0 * (i - 1) / (n - 1)
        u[n+i] = 350.0
        u[2n+i] = 120.0
        u[3n+i] = 4.0
    end
    populate!(ws.state_cache, u, ws.layout)
    Pyrolysis.Properties.update_properties!(ws.prop_cache, ws.state_cache, ws.mesh, mat)
    return problem, mat, ws
end

@testset "Beer-Lambert surface absorptivity" begin
    problem, mat, ws = _radiation_fixture()
    n = ws.layout.n_cells
    A = ws.state_cache.A_section[]
    I0 = 50e3

    Δz = [ws.state_cache.volume[i] / A for i = 1:n]
    τ_tot = sum(ws.prop_cache.α[i] * Δz[i] for i = 1:n)
    ε_surf = ws.prop_cache.ε[n]
    @test 0.85 < ε_surf < 0.95   # mixed composition: strictly between components

    q_bare = zeros(n)
    compute_beer_lambert_source!(q_bare, ws.state_cache, ws.prop_cache, I0, false)
    q_kirchhoff = zeros(n)
    compute_beer_lambert_source!(q_kirchhoff, ws.state_cache, ws.prop_cache, I0, true)

    @testset "Kirchhoff path scales by the live surface ε_eff" begin
        @test q_kirchhoff ≈ ε_surf .* q_bare rtol = 1e-13
        # ε_eff matches the property used for surface *emission* (Kirchhoff).
        ξ_surf = ntuple(k -> ws.state_cache.ξ[k, n], 3)
        @test ε_surf ≈ effective_emissivity(ξ_surf, ws.state_cache.T[n], mat) rtol = 1e-13
    end

    @testset "absorbed power: ε_eff·I·(1 − e^{−τ_tot}), conservation exact" begin
        absorbed = sum(q_kirchhoff[i] * Δz[i] for i = 1:n)
        @test absorbed ≈ ε_surf * I0 * (1 - exp(-τ_tot)) rtol = 1e-12
        # Pre-fix the absorbed power was I0·(1 − e^{−τ_tot}) — the 1/ε_eff
        # model-switch jump. Guard against regression:
        @test !isapprox(absorbed, I0 * (1 - exp(-τ_tot)), rtol = 1e-3)
    end

    @testset "model consistency: BEER_LAMBERT → SURFACE_ABSORPTION as τ → ∞" begin
        # SURFACE_ABSORPTION applies ε_eff·I at the face (boundary Newton).
        # The Beer-Lambert absorbed power must approach the same ε_eff·I for
        # an optically thick sample — NOT 1/ε_eff times it (the old defect).
        absorbed = sum(q_kirchhoff[i] * Δz[i] for i = 1:n)
        @test absorbed / (ε_surf * I0) ≈ 1 - exp(-τ_tot) rtol = 1e-12
        @test absorbed <= ε_surf * I0 + 1e-9
    end

    @testset "explicit-absorptivity path is not double-scaled" begin
        # extract folds α into the intensity and reports it did so.
        bc_explicit = BoundaryConditionSet(
            thermal = (
                top = RadiativeFluxBC(flux = 50e3, absorptivity = 0.8),
                bottom = AdiabaticBC(),
            ),
        )
        I, has_α = extract_radiative_flux_from_bcs(bc_explicit, Val(:top), 0.0)
        @test has_α
        @test I ≈ 0.8 * 50e3
        # Kernel with the flag off leaves the (already-scaled) intensity alone.
        q = zeros(n)
        compute_beer_lambert_source!(q, ws.state_cache, ws.prop_cache, I, false)
        @test q ≈ (0.8 * 50e3 / I0) .* q_bare rtol = 1e-13

        # Kirchhoff default reports no explicit α and returns the bare flux.
        I_k, has_α_k = extract_radiative_flux_from_bcs(problem.bc_set, Val(:top), 0.0)
        @test !has_α_k
        @test I_k ≈ I0
    end

    @testset "dispatcher zeroes non-Beer-Lambert models" begin
        q = fill(1.0, n)
        for model in (NO_RADIATION, SURFACE_ABSORPTION, P1_QUASI_STEADY)
            fill!(q, 1.0)
            compute_radiation_source!(
                q,
                ws.state_cache,
                ws.mesh,
                ws.prop_cache,
                I0,
                model,
                true,
            )
            @test all(iszero, q)
        end
    end
end
