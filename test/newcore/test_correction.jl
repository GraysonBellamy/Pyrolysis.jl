using Test
using LinearAlgebra
using ForwardDiff
using SparseArrays
using Pyrolysis
using Pyrolysis.Materials: DepletionLimiter, BEER_LAMBERT
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.Kernels: RelativeFluxKernel, ReactionKernel, SWITCH_WIDTH, find_kernel
using Pyrolysis.Closures
using Pyrolysis.Assembly: residual!
using Pyrolysis.ColoredJacobian
using Pyrolysis.Deactivation
using Pyrolysis.Integration
using Pyrolysis.Integration: setup, integrate, DomainCheck, CORRECTION_BAND
using LinearSolve: UMFPACKFactorization

# ============================================================================
# Decks
# ============================================================================

# An inert solid on a mesh whose cells grow by `ratio` toward the back.
function _inert(n; L = 0.007, ratio = 2.0, t_end = 3.0)
    a = SolidComponent(:a, ρ = 1000.0, c = 1000.0, k = 0.5)
    return PyrolysisProblem(
        mesh = create_stretched_mesh(L, n, Val(1), ratio),
        material = Material(name = :inert, components = (a,), reactions = ()),
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
        ),
        T_initial = z -> 500.0,
        ξ_initial = [z -> 1000.0],
        tspan = (0.0, t_end),
        output_times = collect(0.0:0.5:t_end),
    )
end

# Two solids of equal density and different heat capacity, a smooth front
# between them (or both uniform), adiabatic and impermeable faces, on a mesh
# whose cells grow by `ratio` toward the back. Each solid holds a quarter of
# the mass or more everywhere: the relative fluxes' face state beside a
# boundary cell is the linear interpolation toward its neighbour, which
# drains a cell holding almost none of a component under a steep front.
function _composition(n; uniform = false, T0 = 800.0, ratio = 1.1)
    L = 0.01
    a = SolidComponent(:a, ρ = 1000.0, c = 1000.0, k = 0.5)
    b = SolidComponent(:b, ρ = 1000.0, c = 3000.0, k = 0.5)
    mat = Material(name = :composition, components = (a, b), reactions = ())
    front(z) = uniform ? 0.5 : 0.5 + 0.25 * tanh((z - L / 2) / (L / 10))
    return PyrolysisProblem(
        mesh = create_stretched_mesh(L, n, Val(2), ratio),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
        ),
        T_initial = z -> T0,
        ξ_initial = [z -> 1000.0 * (1 - front(z)), z -> 1000.0 * front(z)],
        tspan = (0.0, 10.0),
    )
end

# A uniform reaction at a uniform temperature: a uniform strain rate.
function _uniform_slab(n)
    a = SolidComponent(:a, ρ = 1000.0, c = 1000.0, k = 0.5)
    b = SolidComponent(:b, ρ = 400.0, c = 1000.0, k = 0.5)
    g = GaseousComponent(:g, M = 0.03, c = 1000.0, k = 0.03, λ = 0.0)
    mat = Material(
        name = :uniform,
        components = (a, b, g),
        reactions = (
            Reaction(:r, 1 => (2, 3); A = 0.5, E = 0.0, h = 0.0, yields = (0.5, 0.5)),
        ),
        depletion_limiter = DepletionLimiter(enabled = false),
    )
    return PyrolysisProblem(
        mesh = create_stretched_mesh(0.02, n, Val(3), 1.1),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = DirichletThermalBC(400.0), bottom = DirichletThermalBC(400.0)),
        ),
        T_initial = z -> 400.0,
        ξ_initial = [z -> 1000.0, z -> 0.0, z -> 0.0],
        tspan = (0.0, 4.0),
    )
end

# A gasifying layer of `n_top` cells absorbing in depth on an inert substrate
# of `n_bottom` cells, the interface at node `n_bottom + 1` (V-BO-01's
# configuration: three components, in-depth radiation, a moving mesh).
function _two_layer(n_top, n_bottom; t_end = 100.0)
    L_top, L_bottom = 2e-3, 3e-3
    top = SolidComponent(:top, ρ = 1200.0, c = 1500.0, k = 0.2, ε = 0.95, α = 2.0, λ = 1e-6)
    sub = SolidComponent(:sub, ρ = 250.0, c = 1100.0, k = 0.06, ε = 1.0, α = 30.0, λ = 1e-6)
    gas = GaseousComponent(:gas, M = 0.1, c = 1500.0, k = 0.02, λ = 1e-4)
    mat = Material(
        name = :two_layer,
        components = (top, sub, gas),
        reactions = (Reaction(:gasify, 1 => 3; A = 1e11, E = 1.6e5, h = 8e5),),
    )
    z = vcat(
        range(0.0, L_bottom, length = n_bottom + 1),
        range(L_bottom, L_bottom + L_top, length = n_top + 1)[2:end],
    )
    return PyrolysisProblem(
        mesh = create_mesh_from_nodes(z, Val(3)),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = RadiativeFluxBC(flux = 5e4) + RadiativeBC(T_inf = 300.0),
                bottom = AdiabaticBC(),
            ),
            mass = (top = DiffusiveMassBC(), bottom = ImpermeableBC()),
        ),
        T_initial = z -> 600.0 + 5e4 * (z - L_bottom) * (z > L_bottom),
        ξ_initial = [
            z -> z > L_bottom ? 1200.0 : 0.0,
            z -> z > L_bottom ? 0.0 : 250.0,
            z -> 0.0,
        ],
        tspan = (0.0, t_end),
    )
end

# A uniform zero-order burn at a uniform temperature (V-BO-04's
# configuration).
function _zero_order(n; t_end = 480.0)
    T0, ρ_s, r, E = 600.0, 1000.0, 10.0, 5.0e4
    A = r * exp(E / (8.314462618 * T0))
    solid = SolidComponent(:solid, ρ = ρ_s, c = 1000.0, k = 0.2, λ = 1e-6)
    gas = GaseousComponent(:gas, M = 0.03, c = 1000.0, k = 0.02, λ = 1e-5)
    mat = Material(
        name = :zero_order,
        components = (solid, gas),
        reactions = (Reaction(:burn, 1 => 2; A = A, E = E, h = 0.0, n = 0.0),),
    )
    return PyrolysisProblem(
        mesh = create_stretched_mesh(10e-3, n, Val(2), 1.05),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (
                top = ConvectiveMassBC(h_m = 1.0, ξ_inf = 0.0),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> T0,
        ξ_initial = [z -> ρ_s, z -> 0.0],
        tspan = (0.0, t_end),
        output_times = collect(0.0:2.0:t_end),
    )
end

# A heterogeneous layer: two gasifying solids of different densities mixed in
# every cell, the second's fraction graded through the layer, so the density
# varies within the layer and the monitor's mass term acts.
function _mixed(n)
    L = 4e-3
    a = SolidComponent(:a, ρ = 1200.0, c = 1500.0, k = 0.2, ε = 0.95, λ = 1e-6)
    b = SolidComponent(:b, ρ = 600.0, c = 1500.0, k = 0.2, ε = 0.95, λ = 1e-6)
    gas = GaseousComponent(:gas, M = 0.1, c = 1500.0, k = 0.02, λ = 1e-4)
    mat = Material(
        name = :mixed,
        components = (a, b, gas),
        reactions = (
            Reaction(:ra, 1 => 3; A = 1e11, E = 1.6e5, h = 8e5),
            Reaction(:rb, 2 => 3; A = 1e10, E = 1.6e5, h = 8e5),
        ),
    )
    f(z) = 0.1 + 0.4 * z / L
    return PyrolysisProblem(
        mesh = create_stretched_mesh(L, n, Val(3), 1.1),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(5e4), bottom = AdiabaticBC()),
            mass = (top = DiffusiveMassBC(), bottom = ImpermeableBC()),
        ),
        T_initial = z -> 650.0 + 2e4 * z,
        ξ_initial = [z -> 1200.0 * (1 - f(z)), z -> 600.0 * f(z), z -> 0.0],
        tspan = (0.0, 100.0),
    )
end

# An ablating solid heated at the top: its cells thin from the exposed side.
function _ablator(n; L = 3e-3, q = 2.8e6, t_end = 0.8)
    ρ, c, k = 1000.0, 1000.0, 1.0
    R, E, Tp, T0 = 8.314462618, 1.0e5, 700.0, 300.0
    v = 2.0e-3
    A = v^2 * E * (Tp - T0) * exp(E / (R * Tp)) / (k / (ρ * c) * R * Tp^2)
    solid = SolidComponent(:solid, ρ = ρ, c = c, k = k, λ = 1e-6)
    gas = GaseousComponent(:gas, M = 0.05, c = c, k = 0.02, λ = 1e-4)
    mat = Material(
        name = :ablator,
        components = (solid, gas),
        reactions = (Reaction(:ablate, 1 => 2; A = A, E = E, h = 1.0e6),),
    )
    return PyrolysisProblem(
        mesh = create_uniform_mesh(L, n, Val(2)),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(q), bottom = AdiabaticBC()),
            mass = (
                top = ConvectiveMassBC(h_m = 10.0, ξ_inf = 0.0),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> T0 + 1e5 * z,
        ξ_initial = [z -> ρ, z -> 0.0],
        tspan = (0.0, t_end),
        output_times = collect(0.0:0.05:t_end),
    )
end

const _ablation =
    (use_ale = true, deactivation = true, correction = EquidistributionCorrection(τ = 2.0))

# The state of `cp` with the cell thicknesses `h`, from the pinned node.
function _with_thicknesses(cp, h)
    u = copy(cp.u0)
    layout = cp.def.layout
    x = cumsum([0.0; h])
    for nd = 1:n_nodes(layout)
        u[position_index(layout, nd)] = x[nd]
    end
    return u
end

# The node corrections at the state `u`.
function _corrections(cp, u; t = 0.0)
    du = similar(u)
    residual!(du, u, cp.def, t, cp.ws)
    return copy(cp.ws.derived.c)
end

# ============================================================================
# The interface
# ============================================================================

@testset "setup binds the correction to the layers and refuses what it cannot carry" begin
    p = _two_layer(4, 3)
    τ = 100.0
    cp = setup(
        p;
        use_ale = true,
        correction = EquidistributionCorrection(; τ),
        lagrangian = [4],
    )
    corr = cp.def.closure.correction
    @test corr isa LayerCorrection
    @test find_kernel(RelativeFluxKernel, cp.def.kernels) !== nothing
    @test cp.def.setup.topology.lagrangian == [4]
    # The substrate (cells 1 to 3) and the layer (cells 4 to 7).
    @test corr.gain == [0, 9 / τ, 9 / τ, 0, 16 / τ, 16 / τ, 16 / τ, 0]
    @test corr.monitor_scale == ones(7)
    @test cp.def.setup.reference.thickness ≈ diff(cp.def.setup.reference.positions)

    # The interface between the two compositions must be Lagrangian.
    @test_throws ArgumentError setup(
        p;
        use_ale = true,
        correction = EquidistributionCorrection(),
    )
    # Without a correction nothing crosses the faces.
    plain = setup(p; use_ale = true)
    @test plain.def.closure.correction isa NoCorrection
    @test find_kernel(RelativeFluxKernel, plain.def.kernels) === nothing
    @test_throws ArgumentError setup(
        p;
        correction = EquidistributionCorrection(),
        lagrangian = [4],
    )
    @test_throws ArgumentError setup(
        _composition(6);
        use_ale = true,
        w_override = (x, t) -> 0.0,
        correction = EquidistributionCorrection(),
    )
    @test_throws ArgumentError EquidistributionCorrection(τ = 0.0)
    @test_throws ArgumentError EquidistributionCorrection(activation = (0.1, 0.3))
    @test_throws ArgumentError EquidistributionCorrection(monitor = (0.0, 1.0))
    @test_throws ArgumentError EquidistributionCorrection(monitor = (1.0, -1.0))

    # The band defaults to CORRECTION_BAND under the correction, to the
    # threshold without; a zero switch width removes the switch.
    kw = (use_ale = true, deactivation = true, lagrangian = [4])
    @test setup(p; kw..., correction = EquidistributionCorrection()).deactivation.band ==
          CORRECTION_BAND
    @test setup(p; kw...).deactivation.band == 1e-2
    @test setup(p; kw..., threshold = 0.2, correction = EquidistributionCorrection()).deactivation.band ==
          0.2
    @test find_kernel(ReactionKernel, setup(p; kw...).def.kernels).switch.width ==
          SWITCH_WIDTH
    @test find_kernel(ReactionKernel, setup(p; kw..., switch_width = 0).def.kernels).switch ===
          nothing
    @test_throws ArgumentError setup(p; kw..., switch_width = -1.0)
end

# ============================================================================
# §6.1 item 1: the zeros, the sign and the activation's smoothness
# ============================================================================

@testset "the correction is zero at the bounding nodes and moves a node into the thicker cell" begin
    τ = 10.0
    cp = setup(
        _inert(6; ratio = 1.0);
        use_ale = true,
        correction = EquidistributionCorrection(; τ),
        lagrangian = [4],
    )
    layout = cp.def.layout
    h0 = cp.def.setup.reference.thickness
    h = h0 .* [1.0, 1.0, 0.2, 0.05, 1.0, 0.2]
    c = _corrections(cp, _with_thicknesses(cp, h))
    # The pinned, the Lagrangian and the top node; two healthy cells.
    @test c[1] == 0 && c[4] == 0 && c[7] == 0
    @test c[2] == 0
    # A thin cell is fed from its thicker neighbour.
    for nd in (3, 5, 6)
        @test c[nd] != 0
        @test sign(c[nd]) == sign(h[nd] - h[nd-1])
    end
    # Cell 4 is below r_lo: the node beside it is fully active.
    @test c[5] ≈ 9 / τ * (h[5] - h[4]) rtol = 1e-14
    # The outermost live node holds still once the cell beyond it is dead.
    layout.mask[collect(cell_rows(layout, 6))] .= false
    layout.mask[collect(node_rows(layout, 7))] .= false
    @test _corrections(cp, _with_thicknesses(cp, h))[6] == 0
    fill!(layout.mask, true)
    # A deck that never thins below r_hi has c ≡ 0.
    @test all(iszero, _corrections(cp, copy(cp.u0)))
end

@testset "the activation is C¹ in ln h, one at and below zero thickness" begin
    cp = setup(_inert(3); use_ale = true, correction = EquidistributionCorrection())
    corr = cp.def.closure.correction
    h0 = 2e-4
    a(s) = activation(corr, h0 * exp(s), h0)
    da(s) = ForwardDiff.derivative(a, s)
    s_hi, s_lo = log(0.3), log(0.1)
    @test a(0.0) == 0
    @test activation(corr, 0.3h0, h0) == 0
    @test activation(corr, 0.1h0, h0) ≈ 1
    @test a(s_lo - 0.5) == 1
    # Both joins: the value continuous and the slope vanishing on each side.
    for δ in (1e-4, 1e-6)
        @test abs(a(s_hi - δ) - a(s_hi + δ)) <= 10δ^2
        @test abs(da(s_hi - δ)) <= 10δ && da(s_hi + δ) == 0
        @test abs(a(s_lo + δ) - a(s_lo - δ)) <= 10δ^2
        @test abs(da(s_lo + δ)) <= 10δ && da(s_lo - δ) == 0
    end
    # Between the joins the derivative is the finite difference's.
    s = log(0.2)
    @test da(s) ≈ (a(s + 1e-6) - a(s - 1e-6)) / 2e-6 rtol = 1e-6
    @test da(s) < 0
    # At and below zero thickness: one, with a finite derivative.
    for h in (0.0, -h0, 1e-12 * h0)
        @test activation(corr, h, h0) == 1
        @test isfinite(ForwardDiff.derivative(x -> activation(corr, x, h0), h))
    end
end

# ============================================================================
# §6.1 item 2: relaxed equidistribution
# ============================================================================

@testset "three cells relax to equal thickness at the closed-form rate" begin
    τ = 9.0
    # Every node fully active: the switch is one below 50 times the initial
    # thickness, so cₙ = (N²/τ)(hₙ − hₙ₋₁) throughout.
    cp = setup(
        _inert(3);
        use_ale = true,
        correction = EquidistributionCorrection(; τ, activation = (100.0, 50.0)),
        abstol = 1e-12,
        reltol = 1e-10,
    )
    sol = integrate(cp)
    @test sol.retcode === :Success
    k = 9 / τ
    A = k * [-1 1 0; 1 -2 1; 0 1 -1]
    h0 = cp.def.setup.reference.thickness
    @test h0[1] > h0[2] > h0[3]
    for (j, t) in enumerate(sol.t)
        @test diff(sol.z[:, j]) ≈ exp(A * t) * h0 rtol = 1e-7
    end
    # The layer mode decays at k: the spread falls by e^{−k t}.
    spread(j) = (diff(sol.z[:, j])[3] - diff(sol.z[:, j])[1])
    @test spread(length(sol.t)) / spread(1) ≈ exp(-k * sol.t[end]) rtol = 1e-6
    # The solid is at rest and isothermal; its mass moved across the faces.
    @test maximum(abs, sol.T .- 500) <= 1e-9
    @test sum(sol.m[:, :, end]) ≈ sum(sol.m[:, :, 1]) rtol = 1e-12
    @test maximum(abs, sol.ξ[1, :, end] .- 1000.0) <= 1e-9 * 1000
end

# ============================================================================
# §6.1 item 3: the balance rows carry the sweep
# ============================================================================

@testset "w − c is the solid velocity: θ x, θ r/2, θ r/3 with c ≠ 0" begin
    for (coords, d) in ((Planar(), 1), (Cylindrical(), 2), (Spherical(), 3))
        cp = setup(
            _uniform_slab(7);
            use_ale = true,
            coordinates = coords,
            correction = EquidistributionCorrection(τ = 1.0, activation = (2.0, 0.5)),
        )
        layout = cp.def.layout
        u = cp.u0
        du = similar(u)
        residual!(du, u, cp.def, 0.0, cp.ws)
        ds, c = cp.ws.derived, cp.ws.contributions
        θ = c.θ
        @test all(≈(θ[1]; rtol = 1e-13), θ)
        @test θ[1] > 0
        @test maximum(abs, ds.c) > 0.1 * maximum(abs, ds.w)
        for nd = 1:n_nodes(layout)
            x = u[position_index(layout, nd)]
            w = u[velocity_index(layout, nd)]
            @test w - ds.c[nd] ≈ θ[1] * x / d rtol = 1e-12 atol = 1e-18
            # The velocity rows are at round-off after the consistent sweep.
            @test abs(du[velocity_index(layout, nd)]) <= 1e-12 * maximum(abs, ds.w)
        end
        # Each cell's volume follows its strain rate and the sweep of c.
        g = ds.geometry
        for k = 1:n_cells(layout)
            f1, f2 = cp.def.setup.topology.cells[k].faces
            rate(v) = g.face_area[f2] * v[f2] - g.face_area[f1] * v[f1]
            @test rate(ds.w) ≈ θ[k] * ds.V[k] + rate(ds.c) rtol = 1e-11
        end
        # The scaled constraint residual, the sweep in its scale, is at round-off.
        check = DomainCheck(cp)
        @test !check(u, nothing, 0.0)
        @test check.monitor.last <= 1e-6
    end
end

# ============================================================================
# §6.1 item 4: a prescribed velocity is the correction c = w
# ============================================================================

@testset "a prescribed mesh velocity is the correction c = w" begin
    slide = (x, t) -> 2.5e-4 * sinpi(x / 0.01)
    p = _composition(20; ratio = 1.0)
    tol = (abstol = 1e-10, reltol = 1e-8)
    cp = setup(p; use_ale = true, w_override = slide, tol...)
    @test cp.def.closure.correction isa PrescribedCorrection
    @test find_kernel(RelativeFluxKernel, cp.def.kernels) !== nothing
    u = copy(cp.u0)
    du = similar(u)
    residual!(du, u, cp.def, 0.5, cp.ws)
    @test cp.ws.derived.c == cp.ws.derived.w
    # The cell rows read their nodes' velocities through c.
    layout = cp.def.layout
    @test cp.jacobian.pattern[energy_index(layout, 3), velocity_index(layout, 3)]
    a = integrate(cp; saveat = Float64[])
    b = integrate(
        setup(p; use_ale = true, correction = PrescribedCorrection(slide), tol...);
        saveat = Float64[],
    )
    @test a.retcode === :Success
    @test a.raw_solution.t == b.raw_solution.t
    @test a.raw_solution.u == b.raw_solution.u
end

# ============================================================================
# §6.1 item 5: isothermal advection and free stream under the correction
# ============================================================================

@testset "isothermal advection and free stream under the equidistribution correction: $(energy) form" for energy in
                                                                                                          (
    :temperature,
    :enthalpy,
)
    corr = EquidistributionCorrection(τ = 40.0, activation = (2.0, 0.5))
    kw = (energy_state = energy, use_ale = true, correction = corr)
    cp = setup(_composition(20); kw..., abstol = 1e-10, reltol = 1e-8)
    sol = integrate(cp; saveat = Float64[])
    @test sol.retcode === :Success
    @test maximum(abs, sol.T .- 800.0) <= 1e-9
    @test sum(sol.m[:, :, end]) ≈ sum(sol.m[:, :, 1]) rtol = 1e-12
    # The mesh relaxed toward equal cells and the front crossed it.
    h_first, h_last = diff(sol.z[:, 1]), diff(sol.z[:, end])
    @test maximum(h_last) / minimum(h_last) < 0.9 * maximum(h_first) / minimum(h_first)
    @test maximum(abs, sol.ξ[1, :, end] .- sol.ξ[1, :, 1]) > 20
    @test sol.z[1, end] == 0.0
    @test sol.z[end, end] ≈ 0.01 rtol = 1e-12

    # A uniform composition stays uniform.
    cp = setup(_composition(20; uniform = true); kw..., abstol = 1e-10, reltol = 1e-8)
    sol = integrate(cp; saveat = Float64[])
    @test sol.retcode === :Success
    @test maximum(abs, sol.ξ .- 500.0) <= 1e-9 * 500
    @test maximum(abs, sol.T .- 800.0) <= 1e-9
end

# ============================================================================
# §6.1 item 11: containment with the correction active
# ============================================================================

# Every node partially active at the initial state (r_hi > 1 > r_lo), so the
# activation's derivatives enter, and the monitor's mass term on. The upwind
# side of the relative fluxes switches where a node's correction is zero; at
# an active node with equal M̃h on both sides it is exactly zero with nonzero
# partials, and ForwardDiff resolves that tie by the partials, so the colored
# and the dense Jacobian take the kink's derivative along different seeds.
# The states move every interior node by up to 5 %, which breaks every tie.
const _active =
    EquidistributionCorrection(τ = 50.0, activation = (1.5, 0.5), monitor = (1.0, 1.0))
const _correction_configurations = (
    two_layer = (problem = n -> _two_layer(n - 3, 3), lagrangian = [4]),
    zero_order = (problem = _zero_order, lagrangian = Int[]),
    mixed = (problem = _mixed, lagrangian = Int[]),
)

function _contained(cp, u; t = 1.0)
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    du = zeros(length(u))
    residual!(du, u, def, t, ws)
    @test all(isfinite, du)
    Jd = dense_jacobian(def, ws, u, t)
    @test all(isfinite, Jd)
    @test isempty(missing_entries(Jd, cj.pattern))
    J = jacobian!(jac_prototype(cj), cj, u, t)
    @test same_structure(J, cj.pattern)
    @test compare_entries(Matrix(J), Jd; rtol = 1e-12).passed
    return Jd
end

@testset "containment with the correction active: $(name), $(energy) form, $(n) cells" for name in
                                                                                           keys(
        _correction_configurations,
    ),
    energy in (:temperature, :enthalpy),
    n in (6, 8)

    config = _correction_configurations[name]
    cp = setup(
        config.problem(n);
        energy_state = energy,
        use_ale = true,
        deactivation = true,
        radiation_model = name === :two_layer ? BEER_LAMBERT : Pyrolysis.NO_RADIATION,
        correction = _active,
        config.lagrangian,
    )
    layout = cp.def.layout
    u0 = copy(cp.u0)
    for nd = 2:(n_nodes(layout)-1)
        u0[position_index(layout, nd)] *= 1 + 0.05 * sin(nd)
    end
    L = length(u0)
    states = (
        skewed = u0,
        heated = begin
            u = copy(u0)
            for i = 1:n_cells(layout)
                u[energy_index(layout, i)] = energy_state_value(
                    energy_state(layout),
                    cp.def.setup,
                    700.0 + 10i,
                    ntuple(j -> u[mass_index(layout, i, j)], n_components(layout)),
                    i,
                )
            end
            u
        end,
        perturbed = u0 .* (1 .+ 0.2 .* sin.(1:L)),
        zero_state = zeros(L),
        negative = -abs.(u0),
    )
    for (state, u) in pairs(states)
        @testset "$state" begin
            Jd = _contained(cp, u)
            if state === :skewed
                # The correction is active: a balance row reads the positions
                # of the nodes beyond its cell and the masses of its
                # neighbours through the monitor (cell 4's row node is
                # interior to the two-layer deck's upper layer).
                k = name === :two_layer ? 4 : 3
                r = velocity_index(layout, cp.def.closure.row_node[k])
                @test Jd[r, position_index(layout, k + 2)] != 0
                @test Jd[r, mass_index(layout, k + 1, 1)] != 0
            end
        end
    end
end

# ============================================================================
# §6.1 item 14: outputs
# ============================================================================

@testset "outputs: the correction velocity, the thinnest cell and the invocations" begin
    sol = Integration.solve(_ablator(12); _ablation...)
    @test sol.retcode in (:Success, :Terminated)
    @test !isempty(sol.events)
    @test sol.statistics.invocations == last(sol.events).invocation >= 1
    @test all(e -> 1 <= e.invocation <= sol.statistics.invocations, sol.events)
    @test size(sol.correction_velocity) == size(sol.z)
    @test any(x -> isfinite(x) && x != 0, sol.correction_velocity)
    n = size(sol.T, 1)
    for k in eachindex(sol.t)
        live = findall(isfinite, view(sol.T, :, k))
        @test live == 1:length(live)
        top = length(live) + 1
        # Nodes beyond the outermost live one, and dead cells, read NaN.
        @test all(isnan, sol.correction_velocity[(top+1):end, k])
        @test all(isfinite, sol.correction_velocity[1:top, k])
        @test sol.correction_velocity[1, k] == 0
        @test sol.correction_velocity[top, k] == 0
        @test all(isnan, sol.ξ[:, (top):n, k])
        @test all(iszero, sol.m[:, (top):n, k])
        @test sol.extras.min_cell_thickness[k] == minimum(diff(sol.z[1:top, k]))
    end
    # Without a correction the series is zero on a moving mesh and absent on
    # a fixed one.
    plain = Integration.solve(_ablator(6; t_end = 0.1); use_ale = true)
    @test all(iszero, plain.correction_velocity)
    @test plain.statistics.invocations == 0
    fixed = Integration.solve(_ablator(6; t_end = 0.1))
    @test size(fixed.correction_velocity) == (0, 0)
    @test fixed.extras.min_cell_thickness ≈ fill(3e-3 / 6, length(fixed.t)) rtol = 1e-12
end

# ============================================================================
# §6.1 item 16: determinism
# ============================================================================

@testset "repeated solves with the correction and events are bit-identical under each linear solver" begin
    problem = _ablator(12)
    for linear_solver in (StructuredKLU(), UMFPACKFactorization())
        sols = [
            Integration.solve(problem; _ablation..., linear_solver, saveat = Float64[])
            for _ = 1:2
        ]
        @test !isempty(sols[1].events)
        @test sols[2].raw_solution.t == sols[1].raw_solution.t
        @test sols[2].raw_solution.u == sols[1].raw_solution.u
        @test [e.t for e in sols[2].events] == [e.t for e in sols[1].events]
    end
end
