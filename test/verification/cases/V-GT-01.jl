# =============================================================================
# V-GT-01 — Gaussian diffusion of a gas pulse (ThermaKin plane-source case)
# =============================================================================
#
# Pure Fickian transport of gas released inside the solid (Crank/ThermaKin
# plane source). Instead of approximating a delta function, the initial
# condition is a finite-width Gaussian — itself an exact solution — so the
# comparison is IC-discretization-free:
#
#   ξ(z,t) = M_tot/√(2πσ²(t)) · exp(−(z−z₀)²/(2σ²(t))),  σ²(t) = σ₀² + 2dt
#
# Setup notes:
# - All permeabilities are 0 ⇒ `has_darcy_flow` is false ⇒ the residual takes
#   the pure-diffusion path (`compute_gas_fluxes!`), whose kernel
#   −(λ/T)·Δ(ξT)/Δz reduces exactly to Fick's law with diffusivity λ at
#   uniform T (adiabatic faces + no sources keep T uniform to machine).
# - λ = d is set on BOTH components: PARALLEL mixing is volume-fraction
#   weighted, so equal per-component λ makes λ_eff ≡ d independent of
#   composition (otherwise the pulse would see a composition-dependent d).
# - Default ImpermeableBC on both faces: mass is conserved to machine
#   precision; the sealed-boundary image-source correction to the infinite-
#   domain reference is ≤ e^(−11.7) ≈ 8e-6 of the peak at t_end (domain
#   half-width = 4.9σ(t_end)) — below the discretization error.
#
# QoIs: profile L∞ (relative to the evolving peak), total-mass conservation,
# and the parameter-free variance growth law Var(t) = σ₀² + 2dt.

function _verification_case_vgt01()
    d = 1e-6            # gas-transfer coefficient λ [m²/s]
    L, z0, σ0 = 0.04, 0.02, 1e-3
    A_amp = 0.1         # peak bulk concentration [kg/m³]
    T_iso = 400.0
    M_tot = A_amp * sqrt(2π) * σ0
    ts = [2.0, 4.0, 8.0]

    σ2(t) = σ0^2 + 2 * d * t
    ξ_exact(z, t) = M_tot / sqrt(2π * σ2(t)) * exp(-(z - z0)^2 / (2 * σ2(t)))

    function build(n_cells)
        solid = SolidComponent(:solid, ρ = 1000.0, c = 1000.0, k = 0.5, λ = d)
        gas = GaseousComponent(:gas, M = 0.028, c = 1000.0, k = 0.03, λ = d)
        mat = Material(name = :VGT01, components = (solid, gas), reactions = ())
        mesh = create_uniform_mesh(
            L,
            n_cells,
            Val(2);
            T_initial = T_iso,
            ξ_initial = (1000.0, 0.0),
        )
        bc_set =
            BoundaryConditionSet(thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()))
        return PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T_iso,
            ξ_initial = [z -> 1000.0, z -> A_amp * exp(-(z - z0)^2 / (2 * σ0^2))],
            tspan = (0.0, ts[end]),
            output_times = vcat(0.0, ts),
        )
    end

    function qois(sol)
        out = QoIResult[]
        z = cell_positions(sol)
        Δz = L / length(z)
        for (j, t) in enumerate(ts)
            peak = M_tot / sqrt(2π * σ2(t))
            err =
                maximum(abs(sol.ξ[2, i, j+1] - ξ_exact(z[i], t)) for i in eachindex(z)) /
                peak
            push!(
                out,
                QoIResult(
                    name = "L∞ ξ-profile error / peak, t=$(t) s",
                    value = err,
                    error = err,
                    tol = 4e-3,  # 3× observed max: 1.25e-3 at t=2 s (n=160)
                    comment = "Gaussian stays Gaussian; n=160 uniform, σ₀ = 4Δz",
                ),
            )
        end

        # Total gas mass (sealed system): machine conservation.
        m0 = sum(sol.ξ[2, i, 1] for i in eachindex(z)) * Δz
        m_err = maximum(
            abs(sum(sol.ξ[2, i, j] for i in eachindex(z)) * Δz - m0) / m0 for
            j = 2:length(sol.t)
        )
        push!(
            out,
            QoIResult(
                name = "gas-mass conservation (sealed)",
                value = m_err,
                error = m_err,
                tol = 1e-13,  # observed 1.7e-15
                comment = "antisymmetric face fluxes telescope exactly",
            ),
        )

        # Variance growth: discrete second moment vs σ₀² + 2dt (slope 2d is
        # parameter-free). Sealed-boundary truncation affects both identically.
        for (j, t) in enumerate(ts)
            m = sum(sol.ξ[2, i, j+1] for i in eachindex(z)) * Δz
            zbar = sum(sol.ξ[2, i, j+1] * z[i] for i in eachindex(z)) * Δz / m
            var_num = sum(sol.ξ[2, i, j+1] * (z[i] - zbar)^2 for i in eachindex(z)) * Δz / m
            push!(
                out,
                QoIResult(
                    name = "variance growth, t=$(t) s",
                    value = var_num,
                    exact = σ2(t),
                    error = abs(var_num - σ2(t)) / σ2(t),
                    tol = 2e-5,  # 3× observed max: 5.0e-6 at t=8 s (grows w/ boundary truncation)
                    comment = "Var(t) = σ₀² + 2dt; relative error reported",
                ),
            )
        end
        return out
    end

    function l2_error_t4(sol)
        z = cell_positions(sol)
        Δz = L / length(z)
        return sqrt(
            sum((sol.ξ[2, i, 3] - ξ_exact(z[i], 4.0))^2 * Δz for i in eachindex(z)) / L,
        )
    end

    return VerificationCase(
        id = "V-GT-01",
        title = "Gaussian gas-pulse diffusion (plane source)",
        tags = [:gas_transport, :fickian, :mass_conservation],
        references = [
            "J. Crank, The Mathematics of Diffusion (plane-source solution)",
            "Stoliarov & Lyon, FSS 9:1141 (2008), gas-transport verification",
        ],
        build = build,
        solve_kwargs = (abstol = 1e-12, reltol = 1e-10),
        qois = qois,
        ci = CIConfig(
            n_cells = 160,
            budget_s = 20.0,
            convergence = ConvergenceSpec(
                # observed rates 2.009, 2.002
                n_cells = [80, 160, 320],
                error_fn = l2_error_t4,
                expected_order = 2.0,
                min_order = 1.5,
                min_mean_order = 1.8,
                solve_kwargs = (abstol = 1e-12, reltol = 1e-10),
            ),
        ),
        report = (
            n_cells_ladder = [40, 80, 160, 320, 640],
            profile_times = ts,
            qoi_curves = (:xi_profiles, :variance_history),
            curves = sol -> begin
                z = cell_positions(sol)
                Δz = L / length(z)
                function var_num(j)
                    m = sum(sol.ξ[2, i, j] for i in eachindex(z)) * Δz
                    zbar = sum(sol.ξ[2, i, j] * z[i] for i in eachindex(z)) * Δz / m
                    return sum(sol.ξ[2, i, j] * (z[i] - zbar)^2 for i in eachindex(z)) *
                           Δz / m
                end
                [
                    ReportCurve(
                        name = :xi_profiles,
                        title = "V-GT-01: Gaussian gas-pulse diffusion",
                        xlabel = "z (m)",
                        ylabel = "ξ_gas (kg/m³)",
                        series = [
                            CurveSeries(
                                label = "t = $(t) s",
                                x = z,
                                y = sol.ξ[2, :, j+1],
                                exact = zz -> ξ_exact(zz, t),
                            ) for (j, t) in enumerate(ts)
                        ],
                    ),
                    ReportCurve(
                        name = :variance_history,
                        title = "V-GT-01: variance growth Var(t) = σ₀² + 2dt",
                        xlabel = "t (s)",
                        ylabel = "Var (m²)",
                        series = [
                            CurveSeries(
                                label = "discrete second moment",
                                x = sol.t,
                                y = [var_num(j) for j in eachindex(sol.t)],
                                exact = σ2,
                            ),
                        ],
                    ),
                ]
            end,
            error_norms = sol -> begin
                z = cell_positions(sol)
                peak = M_tot / sqrt(2π * σ2(4.0))
                e = [sol.ξ[2, i, 3] - ξ_exact(z[i], 4.0) for i in eachindex(z)]
                (L2 = sqrt(sum(abs2, e) / length(e)) / peak, Linf = maximum(abs, e) / peak)
            end,
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vgt01)
