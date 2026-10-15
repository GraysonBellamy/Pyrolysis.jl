"""
Donor-side and limiter regression tests for the ALE advection scheme.

The ALE correction is ADDED to the Eulerian RHS: dT/dτ = RHS + w·∂T/∂z, i.e.
in standard form dT/dτ + c·∂T/∂z = RHS with apparent transport velocity
c = −w. The donor (upwind) side must therefore follow the sign of −w:
w > 0 → donor on the right (forward difference), w < 0 → donor on the left
(backward difference). Selecting the donor on the sign of w (the historical
defect) is anti-diffusive for both signs — stable-looking on smooth fields
but exponentially unstable on step profiles, which is why the translation
tests below use a step.

The limited gradient must use the near/far slope ratio r = grad_1/grad_0:
with the inverted ratio r = grad_0/grad_1 the correction scales with the
far slope and sign-flips on monotone data (van Leer flips beyond r = 1+√2 —
the leading edge of every steep front).
"""

using Test
using Pyrolysis
using Pyrolysis.Problem: StateCache
using Pyrolysis.Discretization:
    compute_ale_advection_term!, _limited_gradient, MINMOD, VANLEER, SUPERBEE, NONE
using Pyrolysis.Geometry: CellState

@testset "ALE advection donor side" begin

    NC = 2
    L = 0.05
    n_cells = 50
    Δz = L / n_cells

    function make_mesh(T_of_z)
        mesh = create_uniform_mesh(
            L,
            n_cells,
            Val(NC);
            T_initial = 300.0,
            ξ_initial = (1000.0, 0.0),
        )
        for i = 1:n_cells
            cs = mesh.cell_states[i]
            mesh.cell_states[i] =
                CellState{NC}(T_of_z(mesh.cell_centers[i]), cs.ξ, cs.volume, cs.volume_ref)
        end
        return mesh
    end

    set_T!(mesh, T) =
        for i in eachindex(T)
            cs = mesh.cell_states[i]
            mesh.cell_states[i] = CellState{NC}(T[i], cs.ξ, cs.volume, cs.volume_ref)
        end

    @testset "first-order donor selection follows c = -w" begin
        # Non-symmetric smooth field so forward ≠ backward difference.
        mesh = make_mesh(z -> 300.0 + 1e4 * z + 4e5 * z^2)
        n_nodes = length(mesh.nodes)
        state = StateCache(mesh)
        T = [state.T[i] for i = 1:n_cells]

        # w > 0: apparent velocity c = -w < 0 → donor is the RIGHT neighbor.
        w = 2e-5
        advection = zeros(n_cells)
        compute_ale_advection_term!(
            advection,
            state,
            mesh,
            fill(w, n_nodes);
            limiter = NONE,
        )
        for i = 2:(n_cells-1)
            @test advection[i] ≈ w * (T[i+1] - T[i]) / Δz rtol = 1e-12
        end
        # Top cell's donor is the boundary ghost (zero-gradient) → zero.
        @test advection[n_cells] == 0.0
        @test advection[1] ≈ w * (T[2] - T[1]) / Δz rtol = 1e-12

        # w < 0 (standard shrinking): c = -w > 0 → donor is the LEFT neighbor.
        w = -2e-5
        compute_ale_advection_term!(
            advection,
            state,
            mesh,
            fill(w, n_nodes);
            limiter = NONE,
        )
        for i = 2:(n_cells-1)
            @test advection[i] ≈ w * (T[i] - T[i-1]) / Δz rtol = 1e-12
        end
        # Surface corollary: in the shrinking case the surface cell upwinds
        # into the INTERIOR (where ∇T is steepest), not the boundary ghost.
        @test advection[n_cells] ≈ w * (T[n_cells] - T[n_cells-1]) / Δz rtol = 1e-12
        @test advection[1] == 0.0
    end

    @testset "step-profile translation is stable and bounded" begin
        # Pure mesh-motion translation: dT/dτ = w·∂T/∂z with uniform w < 0
        # (standard shrinking sign) translates the profile at c = -w in +z.
        # The anti-diffusive (wrong-donor) scheme blows up exponentially here;
        # a monotone donor-side scheme obeys the discrete maximum principle
        # and does not increase total variation.
        T_hi, T_lo = 400.0, 300.0
        z_front0 = 0.02
        w = -1e-3
        c = -w
        dt = 0.5 * Δz / abs(w)   # CFL = 0.5
        nsteps = 20
        shift = c * dt * nsteps  # = 0.01 m (10 cells)

        total_variation(T) = sum(abs(T[i+1] - T[i]) for i = 1:(length(T)-1))

        front_z(mesh, T) = begin
            i = findfirst(j -> T[j] < 0.5 * (T_hi + T_lo), 1:length(T))
            mesh.cell_centers[i]
        end

        for (limiter, overshoot_tol, tv_tol) in ((NONE, 1e-9, 1e-9), (VANLEER, 1e-9, 1e-9))
            mesh = make_mesh(z -> z < z_front0 ? T_hi : T_lo)
            n_nodes = length(mesh.nodes)
            w_nodes = fill(w, n_nodes)
            T = [mesh.cell_states[i].T for i = 1:n_cells]
            tv0 = total_variation(T)
            advection = zeros(n_cells)

            for _ = 1:nsteps
                set_T!(mesh, T)
                state = StateCache(mesh)
                compute_ale_advection_term!(
                    advection,
                    state,
                    mesh,
                    w_nodes;
                    limiter = limiter,
                )
                @. T += dt * advection
            end

            @test maximum(T) <= T_hi + overshoot_tol
            @test minimum(T) >= T_lo - overshoot_tol
            @test total_variation(T) <= tv0 + tv_tol
            # The front actually translated by ≈ c·t (not zero, not backwards).
            @test abs(front_z(mesh, T) - (z_front0 + shift)) <= 3 * Δz
        end
    end
end

@testset "TVD-limited gradient slope ratio" begin
    # Reference case: monotone data T = {0, 10, 11}, uniform unit
    # spacing, donor on the left → grad_1 = 1, grad_0 = 10. With the inverted
    # ratio (r = grad_0/grad_1 = 10) minmod returned −3.5 and van Leer −7.18:
    # sign-flipped and far outside the [1, 10] hull of the one-sided slopes.
    # The near/far ratio r = grad_1/grad_0 = 0.1 keeps the result bounded and
    # positive.
    @testset "monotone data stays bounded, no sign flip" begin
        g_minmod = _limited_gradient(11.0, 10.0, 0.0, 1.0, 1.0, true, MINMOD, true)
        @test g_minmod ≈ 0.55 atol = 1e-12
        g_vanleer = _limited_gradient(11.0, 10.0, 0.0, 1.0, 1.0, true, VANLEER, true)
        @test g_vanleer ≈ 1.0 + 0.5 * (0.2 / 1.1) * (1.0 - 10.0) atol = 1e-12
        @test g_vanleer > 0.0

        for limiter in (MINMOD, VANLEER, SUPERBEE)
            g = _limited_gradient(11.0, 10.0, 0.0, 1.0, 1.0, true, limiter, true)
            @test 0.0 < g <= 10.0
        end

        # Mirrored right-donor case (monotone decreasing in the donor frame).
        g_right = _limited_gradient(11.0, 10.0, 0.0, 1.0, 1.0, true, MINMOD, false)
        @test g_right ≈ -0.55 atol = 1e-12
    end

    @testset "smooth and degenerate limits" begin
        # Linear field: grad_1 = grad_0 → ψ(1) = 1, correction vanishes.
        for limiter in (MINMOD, VANLEER, SUPERBEE, NONE)
            @test _limited_gradient(2.0, 1.0, 0.0, 1.0, 1.0, true, limiter, true) ≈ 1.0
        end
        # Smooth quadratic T = z² at z = {0, 1, 2}: ψ acts on r = 3 and the
        # van Leer correction stays on the bounded near-slope side.
        g = _limited_gradient(4.0, 1.0, 0.0, 1.0, 1.0, true, VANLEER, true)
        @test g ≈ 3.0 + 0.5 * (6.0 / 4.0) * (3.0 - 1.0) atol = 1e-12  # = 4.5
        # No donor-donor cell → first-order fallback.
        @test _limited_gradient(4.0, 1.0, 0.0, 1.0, 1.0, false, VANLEER, true) ≈ 3.0
        # Flat far slope (grad_0 = 0): ratio undefined → first-order fallback.
        @test _limited_gradient(1.0, 0.0, 0.0, 1.0, 1.0, true, VANLEER, true) ≈ 1.0
        # Extremum (slopes of opposite sign): ψ = 0 → first-order.
        @test _limited_gradient(1.0, 0.0, 1.0, 1.0, 1.0, true, VANLEER, true) ≈ 1.0
    end
end
