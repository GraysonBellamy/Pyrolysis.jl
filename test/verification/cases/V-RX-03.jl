# =============================================================================
# V-RX-03 — Competing parallel reactions A→νB·B + G, A→νC·C + G
# =============================================================================
#
# Two reactions consume the same reactant A with different product splits and
# different heats of reaction (FDS `matl_e_cons` analog). Exact solution:
#
#   ξ_A = ξ₀·e^(−(k₁+k₂)t);  converted mass splits k₁:k₂, scaled by yields:
#   ξ_B = ν_B·(k₁/k)·(ξ₀−ξ_A),  ξ_C = ν_C·(k₂/k)·(ξ₀−ξ_A),
#   ξ_G = (γ_B·k₁ + γ_C·k₂)/k·(ξ₀−ξ_A),   k = k₁+k₂, γ = gas yields
#
# Verifies stoichiometric bookkeeping (two-product yields) and heat-release
# superposition. E = 0 keeps rates exactly T-independent, so the reaction
# heats (h₁ endothermic, h₂ exothermic — h > 0 endothermic in THIS code, as
# in Gpyro and FDS, opposite to ThermaKin) perturb only the temperature, not
# the species solution.
#
# Heat-release superposition target: a single cell pinned by Dirichlet faces
# satisfies the DISCRETE quasi-steady balance exactly:
#
#   2·k_cell·(T_bc − T_c)/(Δz/2) = −Q̇·Δz  ⇒  T_c − T_bc = Q̇·Δz²/(4·k_cell)
#
# with Q̇ = −(h₁k₁ + h₂k₂)·ξ_A(t) the net volumetric source. Quasi-steady is
# valid because the cell diffusion time Δz²/α = 0.02 s ≪ 1/(k₁+k₂) = 2.2 s.
# The gas is retained (default ImpermeableBC, no Darcy/Fickian transport), so
# the system is closed and Σξ = ξ₀ throughout.

function _verification_case_vrx03()
    ξ0 = 1000.0
    k1, k2 = 0.3, 0.15
    ktot = k1 + k2
    νB, γB = 0.7, 0.3   # r1 solid/gas split
    νC, γC = 0.4, 0.6   # r2 solid/gas split
    h1, h2 = +2.0e4, -1.0e4  # J/kg: endo + exo superposition
    T_iso = 400.0
    L = 1e-4
    k_cond, c_sol = 0.5, 1000.0
    ts = [1.0, 2.0, 4.0, 7.0, 10.0, 15.0]

    conv(t) = ξ0 * (1 - exp(-ktot * t))
    ξ_exact(t) = (
        ξ0 * exp(-ktot * t),
        νB * (k1 / ktot) * conv(t),
        νC * (k2 / ktot) * conv(t),
        (γB * k1 + γC * k2) / ktot * conv(t),
    )
    # Net volumetric heat source (h > 0 endothermic ⇒ source −Σhᵢrᵢ):
    Q̇(t) = -(h1 * k1 + h2 * k2) * ξ0 * exp(-ktot * t)

    function build(n_cells)
        compA = SolidComponent(:A, ρ = ξ0, c = c_sol, k = k_cond)
        compB = SolidComponent(:B, ρ = ξ0, c = c_sol, k = k_cond)
        compC = SolidComponent(:C, ρ = ξ0, c = c_sol, k = k_cond)
        gas = GaseousComponent(:G, M = 0.028, c = 1000.0, k = 0.03, λ = 0.0)
        r1 = Reaction(:r1, 1 => (2, 4); A = k1, E = 0.0, h = h1, yields = (νB, γB))
        r2 = Reaction(:r2, 1 => (3, 4); A = k2, E = 0.0, h = h2, yields = (νC, γC))
        mat = Material(
            name = :VRX03,
            components = (compA, compB, compC, gas),
            reactions = (r1, r2),
            depletion_limiter = DepletionLimiter(enabled = false),
        )
        mesh = create_uniform_mesh(L, n_cells)
        bc_set = BoundaryConditionSet(
            thermal = (top = DirichletThermalBC(T_iso), bottom = DirichletThermalBC(T_iso)),
        )
        return PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T_iso,
            ξ_initial = [z -> ξ0, z -> 0.0, z -> 0.0, z -> 0.0],
            tspan = (0.0, ts[end]),
            output_times = vcat(0.0, ts),
        )
    end

    function qois(sol)
        out = QoIResult[]
        for (species, idx) in ((:A, 1), (:B, 2), (:C, 3), (:G, 4))
            err = maximum(
                abs(sol.ξ[idx, 1, j+1] - ξ_exact(t)[idx]) / ξ0 for (j, t) in enumerate(ts)
            )
            push!(
                out,
                QoIResult(
                    name = "max relative ξ_$(species) error",
                    value = err,
                    error = err,
                    tol = 1e-8,  # 3× observed max: 2.3e-9 (ξ_A; abstol=1e-12)
                    comment = "competing-parallel exact split k₁:k₂ × yields",
                ),
            )
        end
        mass_err = maximum(abs(sum(sol.ξ[:, 1, j]) - ξ0) / ξ0 for j = 1:length(sol.t))
        push!(
            out,
            QoIResult(
                name = "total-mass conservation (incl. retained gas)",
                value = mass_err,
                error = mass_err,
                tol = 1e-13,  # observed 8.0e-16
                comment = "yields sum to 1 per reaction; gas sealed in",
            ),
        )
        # Heat-release superposition: quasi-steady discrete drift.
        drift_err = maximum(
            abs((sol.T[1, j+1] - T_iso) - Q̇(t) * L^2 / (4 * k_cond)) for
            (j, t) in enumerate(ts)
        )
        push!(
            out,
            QoIResult(
                name = "reaction-heat superposition (quasi-steady cell drift)",
                value = drift_err,
                error = drift_err,
                tol = 1e-4,  # 3× observed 2.8e-5 K of ~0.02 K amplitude
                comment = "T_c − T_bc = Q̇Δz²/(4k) exactly in the discrete steady balance; " *
                          "drift amplitude ~0.02 K at t=0 (h₁k₁+h₂k₂ superposed)",
            ),
        )
        return out
    end

    return VerificationCase(
        id = "V-RX-03",
        title = "competing parallel reactions with yield splits",
        tags = [:kinetics, :parallel_reactions, :stoichiometry, :energy],
        references = [
            "FDS Verification Guide `matl_e_cons` (energy-superposition analog)",
            "standard parallel first-order kinetics (exact)",
        ],
        build = build,
        solve_kwargs = (abstol = 1e-12, reltol = 1e-10),
        qois = qois,
        ci = CIConfig(n_cells = 1, budget_s = 10.0),
        report = (
            qoi_curves = (:xi_histories, :T_drift_history),
            heating = :isothermal,
            curves = sol -> [
                ReportCurve(
                    name = :xi_histories,
                    title = "V-RX-03: competing parallel reactions, yield splits",
                    xlabel = "t (s)",
                    ylabel = "ξ (kg/m³)",
                    series = [
                        CurveSeries(
                            label = "ξ_$(species)",
                            x = sol.t,
                            y = sol.ξ[idx, 1, :],
                            exact = t -> ξ_exact(t)[idx],
                        ) for (species, idx) in ((:A, 1), (:B, 2), (:C, 3), (:G, 4))
                    ],
                ),
                ReportCurve(
                    name = :T_drift_history,
                    title = "V-RX-03: reaction-heat superposition (quasi-steady drift)",
                    xlabel = "t (s)",
                    ylabel = "T − T_iso (K)",
                    series = [
                        CurveSeries(
                            label = "cell drift",
                            x = sol.t,
                            y = sol.T[1, :] .- T_iso,
                            exact = t -> Q̇(t) * L^2 / (4 * k_cond),
                        ),
                    ],
                ),
            ],
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vrx03)
