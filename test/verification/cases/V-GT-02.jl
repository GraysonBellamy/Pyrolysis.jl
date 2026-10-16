# =============================================================================
# V-GT-02 — Steady compressible Darcy outflow with uniform generation (P²)
# =============================================================================
#
# The canonical porous-flow exact solution (gas-reservoir form): isothermal
# ideal gas, uniform volumetric generation ġ (zero-order gas-producing
# reaction at clamped T), sealed back at z = 0, Dirichlet pressure P∞ at
# z = L. The steady mass balance d/dz[(κM/μRT)·P·dP/dz] = −ġ integrates to
#
#   P²(z) = P∞² + (μRT·ġ/(κM))·(L² − z²),    exit mass flux = ġ·L  (exact)
#
# Nonlinear in P; exercises the algebraic ideal-gas pressure closure, face
# mobility, Darcy upwinding, and the DirichletPressureBC outflow — all at
# once. The existing Darcy test only checks v = −(κ/μ)∇P pointwise.
#
# Setup notes (the four pitfalls the survey flagged):
# - P is NOT persisted in PyrolysisSolution: reconstructed per cell from
#   sol.ξ/sol.T via porosity + ideal_gas_pressure with the same material.
# - Initial gas is seeded with background_gas_concentration so P(z,0) = P∞
#   (starting from ξ_gas = 0 would start from vacuum).
# - μ in the exact profile is the MATERIAL's own viscosity curve,
#   effective_gas_viscosity(T_iso, mat), not a hand-picked constant — so the
#   case stays exact for whatever `gas_viscosity` the material carries, which
#   is what lets the ×2-viscosity check in test_darcy.jl exercise that field.
# - λ = 0 on every component kills the Fickian term: pure Darcy outflow.
# - ρ_intrinsic = 1e5 (φ = 0.99, verification-only value): solid→gas
#   conversion frees pore volume at dφ/dt = ġ/ρ_intrinsic, and the gas
#   stored in that new volume steals (P̄M/RT)/ρ_intrinsic of the exit flux —
#   a mesh-independent bias (~3.6e-4 at ρ_intrinsic = 2000) that would floor
#   the convergence study; at 1e5 it is ~7e-6, below discretization error.
#   Permeability mixes on the renormalized solid basis, so κ_eff = κ exactly
#   regardless of φ.
#
# Timescales: pressure equilibrates in ≪ 0.1 s (gas inventory / production
# ≈ 0.07 s); compared at t = 1 s ≈ steady, with an explicit steadiness QoI.

function _verification_case_vgt02()
    κ_perm = 1e-12
    L = 0.01
    T_iso = 500.0
    M_gas = 0.028
    P∞ = 101325.0
    ġ = 5.0                 # zero-order generation [kg/(m³·s)]
    R = 8.314462618
    t_end = 1.0

    function build(n_cells)
        solid = SolidComponent(
            :solid,
            ρ = 1000.0,
            c = 1000.0,
            k = 0.5,
            κ = κ_perm,
            λ = 0.0,
            ρ_intrinsic = 1e5,
        )
        gas = GaseousComponent(:gas, M = M_gas, c = 1000.0, k = 0.03, λ = 0.0)
        rxn = Reaction(:generate, 1 => 2; A = ġ, E = 0.0, h = 0.0, n = 0.0)
        mat = Material(
            name = :VGT02,
            components = (solid, gas),
            reactions = (rxn,),
            depletion_limiter = DepletionLimiter(enabled = false),
        )
        ξ_g0 = background_gas_concentration(mat, 2, (1000.0, 0.0), T_iso; P = P∞)
        mesh = create_uniform_mesh(
            L,
            n_cells,
            Val(2);
            T_initial = T_iso,
            ξ_initial = (1000.0, ξ_g0),
        )
        bc_set = BoundaryConditionSet(
            thermal = (top = DirichletThermalBC(T_iso), bottom = DirichletThermalBC(T_iso)),
            mass = (top = DirichletPressureBC(P∞), bottom = ImpermeableBC()),
        )
        return PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T_iso,
            ξ_initial = [z -> 1000.0, z -> ξ_g0],
            tspan = (0.0, t_end),
            output_times = [0.0, t_end / 2, t_end],
        )
    end

    function qois(sol)
        out = QoIResult[]
        mat = sol.problem.material
        z = get_cell_positions(sol)
        μ = effective_gas_viscosity(T_iso, mat)
        ΔP²_coeff = μ * R * T_iso * ġ / (κ_perm * M_gas)
        P²_exact(zz) = P∞^2 + ΔP²_coeff * (L^2 - zz^2)
        ΔP²_scale = ΔP²_coeff * L^2  # back-face overpressure² scale

        P_num = map(eachindex(z)) do i
            ξi = (sol.ξ[1, i, end], sol.ξ[2, i, end])
            Ti = sol.T[i, end]
            φi = porosity(ξi, Ti, mat)
            ideal_gas_pressure(ξi, Ti, φi, mat)
        end

        err_P² = maximum(abs(P_num[i]^2 - P²_exact(z[i])) / ΔP²_scale for i in eachindex(z))
        push!(
            out,
            QoIResult(
                name = "P² profile error / overpressure² scale",
                value = err_P²,
                error = err_P²,
                tol = 2.5e-3,  # 3× observed 7.1e-4 (n=80)
                comment = "steady gas-reservoir parabola in P²; n=80 uniform",
            ),
        )

        P_back_exact = sqrt(P²_exact(0.0))
        push!(
            out,
            QoIResult(
                name = "back-face pressure",
                value = P_num[1],
                exact = P_back_exact,
                error = abs(P_num[1] - P_back_exact) / (P_back_exact - P∞),
                tol = 2.5e-3,  # 3× observed 7.2e-4 relative to overpressure
                comment = "overpressure ≈ $(round(P_back_exact - P∞, digits = 0)) Pa " *
                          "(wood-pyrolysis magnitude)",
            ),
        )

        flux_num = sol.extras.surface_mass_flux[end]
        flux_exact = ġ * L
        push!(
            out,
            QoIResult(
                name = "exit mass flux vs ġL",
                value = flux_num,
                exact = flux_exact,
                error = abs(flux_num - flux_exact) / flux_exact,
                tol = 5e-5,  # 3× observed 7.2e-6 (residual φ-growth storage bias)
                comment = "parameter-free exact steady outflow",
            ),
        )

        steadiness =
            abs(sol.extras.surface_mass_flux[end] - sol.extras.surface_mass_flux[end-1]) /
            flux_exact
        push!(
            out,
            QoIResult(
                name = "steadiness of exit flux",
                value = steadiness,
                error = steadiness,
                tol = 1e-6,  # observed 2.6e-8
                comment = "fill timescale ≈ 0.07 s ≪ t_end = 1 s",
            ),
        )
        return out
    end

    # L∞ of the P² profile (normalized) for the order check, at t_end.
    function p2_error(sol)
        mat = sol.problem.material
        z = get_cell_positions(sol)
        μ = effective_gas_viscosity(T_iso, mat)
        ΔP²_coeff = μ * R * T_iso * ġ / (κ_perm * M_gas)
        maximum(eachindex(z)) do i
            ξi = (sol.ξ[1, i, end], sol.ξ[2, i, end])
            φi = porosity(ξi, T_iso, mat)
            P = ideal_gas_pressure(ξi, T_iso, φi, mat)
            abs(P^2 - (P∞^2 + ΔP²_coeff * (L^2 - z[i]^2))) / (ΔP²_coeff * L^2)
        end
    end

    return VerificationCase(
        id = "V-GT-02",
        title = "steady compressible Darcy outflow (P² profile)",
        tags = [:gas_transport, :darcy, :pressure_closure, :steady],
        references = [
            "gas-reservoir P² solution (Muskat-form compressible Darcy, standard porous-media texts)",
            "technical reference Ch. 9 (wood-pyrolysis overpressure magnitudes)",
        ],
        build = build,
        solve_kwargs = (abstol = 1e-10, reltol = 1e-8),
        qois = qois,
        ci = CIConfig(
            n_cells = 80,
            budget_s = 30.0,
            convergence = ConvergenceSpec(
                # first-order upwind advection; observed rates 0.89, 0.94
                n_cells = [40, 80, 160],
                error_fn = p2_error,
                expected_order = 1.0,
                min_order = 0.7,
                min_mean_order = 0.8,
                solve_kwargs = (abstol = 1e-10, reltol = 1e-8),
            ),
        ),
        report = (
            n_cells_ladder = [20, 40, 80, 160, 320],
            profile_times = [t_end],
            qoi_curves = (:P2_profile, :exit_flux_history),
            curves = sol -> begin
                mat = sol.problem.material
                z = get_cell_positions(sol)
                μ = effective_gas_viscosity(T_iso, mat)
                ΔP²_coeff = μ * R * T_iso * ġ / (κ_perm * M_gas)
                P_num = map(eachindex(z)) do i
                    ξi = (sol.ξ[1, i, end], sol.ξ[2, i, end])
                    φi = porosity(ξi, sol.T[i, end], mat)
                    ideal_gas_pressure(ξi, sol.T[i, end], φi, mat)
                end
                [
                    ReportCurve(
                        name = :P2_profile,
                        title = "V-GT-02: steady gas-reservoir parabola in P²",
                        xlabel = "z (m)",
                        ylabel = "P² (Pa²)",
                        series = [
                            CurveSeries(
                                label = "reconstructed cell pressure²",
                                x = z,
                                y = P_num .^ 2,
                                exact = zz -> P∞^2 + ΔP²_coeff * (L^2 - zz^2),
                            ),
                        ],
                    ),
                    ReportCurve(
                        name = :exit_flux_history,
                        title = "V-GT-02: exit mass flux vs exact ġL",
                        xlabel = "t (s)",
                        ylabel = "ṁ″ (kg/m²s)",
                        series = [
                            CurveSeries(
                                label = "surface mass flux",
                                x = sol.t,
                                y = sol.extras.surface_mass_flux,
                                exact = t -> ġ * L,
                            ),
                        ],
                    ),
                ]
            end,
            error_norms = sol -> begin
                mat = sol.problem.material
                z = get_cell_positions(sol)
                μ = effective_gas_viscosity(T_iso, mat)
                ΔP²_coeff = μ * R * T_iso * ġ / (κ_perm * M_gas)
                e = map(eachindex(z)) do i
                    ξi = (sol.ξ[1, i, end], sol.ξ[2, i, end])
                    φi = porosity(ξi, T_iso, mat)
                    P = ideal_gas_pressure(ξi, T_iso, φi, mat)
                    (P^2 - (P∞^2 + ΔP²_coeff * (L^2 - z[i]^2))) / (ΔP²_coeff * L^2)
                end
                (L2 = sqrt(sum(abs2, e) / length(e)), Linf = maximum(abs, e))
            end,
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vgt02)
