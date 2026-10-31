# =============================================================================
# V-RX-06 — conversion-basis burnout on the moving mesh
# =============================================================================
#
# An adiabatic uniform slab at T₀ with one reaction solid → gas of order n and
# no heat of reaction, so T stays at T₀, k = A·exp(−E/RT₀) is a constant and the
# conversion x = m/m₀ obeys dx/dt = −k xⁿ. The order acts on the conversion
# whatever the mesh does: on the moving mesh, where the concentration stays at
# its bulk density and a concentration-based order would behave as first order,
# the mass and the thickness follow
#
#   m/m₀ = L/L₀ = [1 − (1−n) k t]^(1/(1−n)).
#
# Sub-cases at n = 0, 0.5 and 2, on the moving mesh and on a fixed one, the
# orders below one compared up to 0.9·t_b. The limiter's onset sits at 1e-6, so
# over the window it is exactly one and the finite base rounds away.

function _vrx06_setup(n_cells, order; t_end = _vrx06_window(order), x_th = 1e-6)
    ρ_s, T0, L0, k, E = 500.0, 600.0, 10e-3, 0.02, 5.0e4
    A = k * exp(E / (8.314462618 * T0))
    solid = SolidComponent(:solid, ρ = ρ_s, φ = 0.3, c = 1000.0, k = 0.2, λ = 1e-6)
    gas = GaseousComponent(:gas, M = 0.03, c = 1000.0, k = 0.02)
    material = Material(
        name = :VRX06,
        components = (solid, gas),
        reactions = (Reaction(:burn, 1 => 2; A = A, E = E, h = 0.0, n = order),),
        depletion_limiter = DepletionLimiter(conversion_threshold = x_th),
    )
    bc_set = BoundaryConditionSet(
        thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
        mass = (
            top = AmbientGasBC(GasReservoir(P = 0.0, T = T0)),
            bottom = ImpermeableBC(),
        ),
    )
    problem = PyrolysisProblem(
        mesh = create_uniform_mesh(L0, n_cells),
        material = material,
        bc_set = bc_set,
        T_initial = z -> T0,
        ξ_initial = [z -> ρ_s, z -> 0.0],
        tspan = (0.0, t_end),
        output_times = collect(range(0.0, t_end, length = 16)),
    )
    return (; problem, T0, L0, k, order)
end

_vrx06_window(order) = order < 1 ? 0.9 * nth_order_burnout_time(k = 0.02, n = order) : 150.0

_vrx06_exact(S) = t -> nth_order_xi(t; ξ0 = 1.0, k = S.k, n = S.order)

# The condensed mass of a solution against its value at the start.
_vrx06_mass(sol) =
    [sum(view(sol.m, 1, :, j)) / sum(view(sol.m, 1, :, 1)) for j in eachindex(sol.t)]

function _vrx06_case(id, order, mass_tol, fixed_tol, thickness_tol, difference_tol)
    n_default = 20
    keywords = (abstol = 1e-10, reltol = 1e-8)
    build(n) = _vrx06_setup(n, order).problem
    function qois(sol)
        S = _vrx06_setup(n_default, order)
        x = _vrx06_exact(S)
        fixed = suite_solve(S.problem; keywords...)
        moving = _vrx06_mass(sol)
        at_rest = _vrx06_mass(fixed)
        mass = maximum(abs(moving[j] - x(sol.t[j])) / x(sol.t[j]) for j in eachindex(sol.t))
        rest =
            maximum(abs(at_rest[j] - x(sol.t[j])) / x(sol.t[j]) for j in eachindex(sol.t))
        L = sol.extras.thickness_history
        thickness =
            maximum(abs(L[j] / S.L0 - x(sol.t[j])) / x(sol.t[j]) for j in eachindex(sol.t))
        difference = maximum(abs(moving[j] - at_rest[j]) for j in eachindex(sol.t))
        drift = maximum(abs(T - S.T0) for T in sol.T)
        return [
            QoIResult(
                name = "m/m₀ on the moving mesh",
                value = moving[end],
                exact = x(sol.t[end]),
                error = mass,
                tol = mass_tol,
                comment = "max relative error over the saved times",
            ),
            QoIResult(
                name = "m/m₀ on the fixed mesh",
                value = at_rest[end],
                exact = x(sol.t[end]),
                error = rest,
                tol = fixed_tol,
                comment = "max relative error over the saved times",
            ),
            QoIResult(
                name = "L/L₀ on the moving mesh",
                value = L[end] / S.L0,
                exact = x(sol.t[end]),
                error = thickness,
                tol = thickness_tol,
                comment = "max relative error over the saved times",
            ),
            QoIResult(
                name = "the two meshes agree",
                value = difference,
                exact = 0.0,
                error = difference,
                tol = difference_tol,
                comment = "largest difference of m/m₀, moving against fixed",
            ),
            QoIResult(
                name = "isothermal: T stays at T₀",
                value = S.T0,
                exact = S.T0,
                error = drift,
                tol = 1.4e-11,  # 3× observed 4.7e-12
                comment = "K; no heat of reaction, no gradient",
            ),
        ]
    end
    return VerificationCase(
        id = id,
        title = "conversion-basis burnout on the moving mesh (n = $(order))",
        tags = [:ale, :kinetics, :reaction_order, :exact],
        references = ["closed form m/m₀ = L/L₀ = [1 − (1−n) k t]^(1/(1−n))"],
        build = build,
        solve_kwargs = (use_ale = true, keywords...),
        qois = qois,
        # The fixed-mesh solve in the QoIs compiles inside the timed window.
        ci = CIConfig(n_cells = n_default, budget_s = 20.0),
        report = (
            curves = sol -> begin
                S = _vrx06_setup(n_default, order)
                [
                    ReportCurve(
                        name = :mass_history,
                        title = "$(id): mass and thickness at order $(order)",
                        xlabel = "t (s)",
                        ylabel = "fraction of the initial value",
                        series = [
                            CurveSeries(
                                label = "m/m₀, moving mesh",
                                x = sol.t,
                                y = _vrx06_mass(sol),
                                exact = _vrx06_exact(S),
                            ),
                            CurveSeries(
                                label = "L/L₀, moving mesh",
                                x = sol.t,
                                y = sol.extras.thickness_history ./ S.L0,
                                exact = _vrx06_exact(S),
                            ),
                        ],
                    ),
                ]
            end,
        ),
    )
end

# Tolerances at 3× the larger observed of the two energy forms: moving mass,
# fixed mass, thickness, the meshes' difference.
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _vrx06_case("V-RX-06a", 0.0, 5e-9, 5e-9, 5e-9, 2e-14),  # 1.7e-9, 1.7e-9, 1.7e-9, 6.4e-15
)
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _vrx06_case("V-RX-06b", 0.5, 1e-8, 2.3e-7, 1e-8, 2.2e-9),  # 3.3e-9, 7.6e-8, 3.3e-9, 7.3e-10
)
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _vrx06_case("V-RX-06c", 2.0, 4.2e-10, 2.3e-9, 4.2e-10, 5e-10),  # 1.4e-10, 7.6e-10, 1.4e-10, 1.6e-10
)
