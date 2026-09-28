# =============================================================================
# V-BO-04 — uniform zero-order burnout (exact L(t) through the park cascade)
# =============================================================================
#
# An adiabatic slab at uniform T₀ with one zero-order reaction solid → gas,
# γ_gas = 0 and no heat of reaction. The temperature stays at T₀ exactly (no
# heat source, no gradient), so the volumetric rate r = A·exp(−E/RT₀) is a
# constant and every cell loses the same fraction r/ρ_s of its current
# volume per unit time. Under ALE the solid concentration stays at ρ_s
# (dilation removes volume, not concentration) and the exact thickness is
#
#   L(t) = L₀ · exp(−r t / ρ_s).
#
# Every cell reaches the depletion threshold at the same instant, so the
# parking machinery meets its worst case: a cascade of park events, one per
# step, surface-first, halving the live cell count at t_k = ln(2^k·20) ρ_s/r.
# L(t) and the solid mass must follow the exact law through every cascade,
# and ξ_s must stay at ρ_s to roundoff.

function _vbo04_setup(n_cells::Int)
    ρ_s, c, k = 1000.0, 1000.0, 0.2
    T0 = 600.0
    L0 = 10e-3
    r = 10.0                                     # kg m⁻³ s⁻¹ volumetric rate at T0
    E = 5.0e4
    R = 8.314462618
    A = r * exp(E / (R * T0))
    t_end = 480.0                                # three cascades (40 → 20 → 10 → 5 cells)
    solid = SolidComponent(:solid, ρ = ρ_s, c = c, k = k, λ = 1e-6)  # non-degenerate gas transport
    gas = GaseousComponent(:gas, M = 0.03, c = 1000.0, k = 0.02, λ = 1e-5)
    rxn = Reaction(:burn, 1 => 2; A = A, E = E, h = 0.0, n = 0.0)
    mat = Material(name = :VBO04, components = (solid, gas), reactions = (rxn,))
    mesh = create_uniform_mesh(L0, n_cells, Val(2); T_initial = T0, ξ_initial = (ρ_s, 0.0))
    bcs = BoundaryConditionSet(
        thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
        mass = (top = ConvectiveMassBC(h_m = 1.0, ξ_inf = 0.0), bottom = ImpermeableBC()),
    )
    problem = PyrolysisProblem(
        mesh = mesh,
        material = mat,
        bc_set = bcs,
        T_initial = z -> T0,
        ξ_initial = [z -> ρ_s, z -> 0.0],
        tspan = (0.0, t_end),
        output_times = collect(0.0:2.0:t_end),
    )
    return (; problem, ρ_s, T0, L0, r, t_end)
end

function _verification_case_vbo04()
    n_default = 40
    function build(n)
        return _vbo04_setup(n).problem
    end
    function qois(sol)
        S = _vbo04_setup(n_default)
        L_exact = t -> S.L0 * exp(-S.r * t / S.ρ_s)
        th = sol.extras.thickness_history
        err_L = maximum(
            abs(th[i] - L_exact(sol.t[i])) / L_exact(sol.t[i]) for i in eachindex(sol.t)
        )
        # Solid concentration stays at ρ_s in every live cell.
        ξs = sol.ξ[1, :, :]
        err_ξ = maximum(abs(x - S.ρ_s) / S.ρ_s for x in ξs if !isnan(x))
        # Temperature stays at T0.
        err_T = maximum(abs(x - S.T0) for x in sol.T if !isnan(x))
        n_parks = length(sol.extras.park_times)
        live_end = sol.extras.live_cells[end]
        return [
            QoIResult(
                name = "thickness follows L₀ exp(−r t/ρ) through the park cascades",
                value = th[end],
                exact = L_exact(sol.t[end]),
                error = err_L,
                tol = 1e-4,
                comment = "max relative error over all output times",
            ),
            QoIResult(
                name = "solid concentration stays at ρ_s (dilation removes volume only)",
                value = S.ρ_s,
                exact = S.ρ_s,
                error = err_ξ,
                tol = 1e-8,
                comment = "max |ξ_s − ρ_s|/ρ_s over live cells and times",
            ),
            QoIResult(
                name = "isothermal: T stays at T₀",
                value = S.T0,
                exact = S.T0,
                error = err_T,
                tol = 1e-6,
                comment = "K; no heat of reaction, no gradient",
            ),
            QoIResult(
                name = "park cascade: 40 → 5 live cells",
                value = Float64(live_end),
                exact = 5.0,
                error = abs(live_end - 5.0),
                tol = 0.5,
                comment = "$(n_parks) park events; every cell hits the threshold at once",
            ),
        ]
    end
    return VerificationCase(
        id = "V-BO-04",
        title = "Uniform zero-order burnout: exact thickness through park cascades",
        tags = [:ale, :burnout, :depletion, :kinetics, :exact],
        references = ["closed form L(t) = L₀ exp(−r t/ρ_s)"],
        build = build,
        solve_kwargs = (
            use_ale = true,
            handle_depletion = true,
            depletion_threshold = 0.05,
            min_cells = 2,
            min_thickness = 1e-6,
            abstol = 1e-9,
            reltol = 1e-7,
        ),
        qois = qois,
        ci = CIConfig(n_cells = n_default, budget_s = 10.0),
        # Single rung: the solution is uniform in space, so a resolution
        # ladder has no spatial error to converge.
        report = (
            curves = sol -> begin
                S = _vbo04_setup(n_default)
                [
                    ReportCurve(
                        name = :thickness_history,
                        title = "V-BO-04: thickness through three park cascades",
                        xlabel = "t (s)",
                        ylabel = "thickness (m)",
                        series = [
                            CurveSeries(
                                label = "live thickness",
                                x = sol.t,
                                y = sol.extras.thickness_history,
                                exact = t -> S.L0 * exp(-S.r * t / S.ρ_s),
                            ),
                        ],
                    ),
                    ReportCurve(
                        name = :live_cells_history,
                        title = "V-BO-04: live cells",
                        xlabel = "t (s)",
                        ylabel = "live cells",
                        series = [
                            CurveSeries(
                                label = "live cells",
                                x = sol.t,
                                y = Float64.(sol.extras.live_cells),
                            ),
                        ],
                    ),
                ]
            end,
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vbo04)
