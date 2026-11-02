# =============================================================================
# V-RX-07 — series scheme with an order other than one on the intermediate
# =============================================================================
#
# One cell held at 400 K by Dirichlet faces, A → 0.9 B + 0.1 gas at first order,
# k₁ = 0.05 s⁻¹, and B → gas at order 1.8, k₂ = 0.1 s⁻¹, with E = 0 and no limiter.
# The intermediate's reference is its lineage inventory, 0.9 of the mass the
# cell started with, so its conversion e_B = m_B/(0.9 m_A,0) obeys the Kinetics
# Neo balance on conversions,
#
#   de_A/dt = −k₁ e_A,    de_B/dt = k₁ e_A − k₂ e_B^1.8,
#
# integrated by `series_chain_conversion`. With the limiter off only the finite
# base of the conversion remains, a sink of about 2e-15 k₂ y m_A,0 at e_B = 0.

function _vrx07_problem(n_cells)
    ξ0, T_iso, y = 1000.0, 400.0, 0.9
    A = SolidComponent(:A, ρ = ξ0, φ = 0.3, c = 1000.0, k = 0.5, λ = 1e-6)
    B = SolidComponent(:B, ρ = ξ0, φ = 0.3, c = 1000.0, k = 0.5, λ = 1e-6)
    gas = GaseousComponent(:gas, M = 0.03, c = 1000.0, k = 0.02)
    material = Material(
        name = :VRX07,
        components = (A, B, gas),
        reactions = (
            Reaction(:first, 1 => (2, 3); A = 0.05, E = 0.0, h = 0.0, yields = (y, 1 - y)),
            Reaction(:second, 2 => 3; A = 0.1, E = 0.0, h = 0.0, n = 1.8),
        ),
        depletion_limiter = DepletionLimiter(enabled = false),
    )
    return PyrolysisProblem(
        mesh = uniform_mesh(1e-3, n_cells),
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (top = DirichletThermalBC(T_iso), bottom = DirichletThermalBC(T_iso)),
            mass = (
                top = AmbientGasBC(GasReservoir(P = 0.0, T = T_iso)),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> T_iso,
        ξ_initial = [z -> ξ0, z -> 0.0, z -> 0.0],
        tspan = (0.0, 80.0),
        output_times = [0.0, 5.0, 10.0, 20.0, 40.0, 80.0],
    )
end

function _verification_case_vrx07()
    y = 0.9
    exact(t) = series_chain_conversion(t; k1 = 0.05, k2 = 0.1, n2 = 1.8)
    function qois(sol)
        m0 = sol.m[1, 1, 1]
        times = eachindex(sol.t)[2:end]
        A = maximum(abs(sol.m[1, 1, j] / m0 - exact(sol.t[j])[1]) for j in times)
        B = maximum(abs(sol.m[2, 1, j] / (y * m0) - exact(sol.t[j])[2]) for j in times)
        total = maximum(
            abs(
                (sol.m[1, 1, j] + sol.m[2, 1, j]) / m0 - exact(sol.t[j])[1] -
                y * exact(sol.t[j])[2],
            ) for j in times
        )
        return [
            QoIResult(
                name = "conversion of A",
                value = sol.m[1, 1, end] / m0,
                exact = exact(sol.t[end])[1],
                error = A,
                tol = 2.4e-12,  # 3× observed 7.7e-13
                comment = "max absolute error over 5, 10, 20, 40, 80 s",
            ),
            QoIResult(
                name = "conversion of B on its lineage",
                value = sol.m[2, 1, end] / (y * m0),
                exact = exact(sol.t[end])[2],
                error = B,
                tol = 3.7e-11,  # 3× observed 1.2e-11
                comment = "max absolute error over 5, 10, 20, 40, 80 s",
            ),
            QoIResult(
                name = "condensed mass",
                value = (sol.m[1, 1, end] + sol.m[2, 1, end]) / m0,
                exact = exact(sol.t[end])[1] + y * exact(sol.t[end])[2],
                error = total,
                tol = 3.3e-11,  # 3× observed 1.1e-11
                comment = "max absolute error over 5, 10, 20, 40, 80 s",
            ),
        ]
    end
    return VerificationCase(
        id = "V-RX-07",
        title = "series scheme with an order other than one on the intermediate",
        tags = [:kinetics, :reaction_order, :isothermal],
        references = [
            "NETZSCH Kinetics Neo FAQ, What Equations Are Used in My Kinetic Model?",
        ],
        build = _vrx07_problem,
        solve_kwargs = (abstol = 1e-11, reltol = 1e-9),
        qois = qois,
        ci = CIConfig(n_cells = 1, budget_s = 10.0),
        report = (
            curves = sol -> begin
                m0 = sol.m[1, 1, 1]
                [
                    ReportCurve(
                        name = :conversions,
                        title = "V-RX-07: the series scheme on conversions",
                        xlabel = "t (s)",
                        ylabel = "conversion",
                        series = [
                            CurveSeries(
                                label = "A",
                                x = sol.t,
                                y = sol.m[1, 1, :] ./ m0,
                                exact = t -> t == 0 ? 1.0 : exact(t)[1],
                            ),
                            CurveSeries(
                                label = "B on its lineage",
                                x = sol.t,
                                y = sol.m[2, 1, :] ./ (y * m0),
                                exact = t -> t == 0 ? 0.0 : exact(t)[2],
                            ),
                        ],
                    ),
                ]
            end,
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vrx07)
