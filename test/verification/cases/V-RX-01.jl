# =============================================================================
# V-RX-01 — Two-step series chain A→B→C (the FDS/Gpyro shared benchmark)
# =============================================================================
#
# Identical problem in FDS (`two_step_solid_reaction`) and Gpyro §3.4.4:
# constant rates k₁ = 0.389 s⁻¹, k₂ = 0.262 s⁻¹ (published values, so numbers
# are cross-checkable against both codes), isothermal, exact Bateman solution
# in exact/kinetics.jl.
#
# Exercises multi-step reaction coupling: intermediate B is produced by
# reaction 1 and consumed by reaction 2 SIMULTANEOUSLY (the existing
# isothermal test is single-step).
#
# Setup notes (solver conventions that would corrupt the case):
# - constant rates realized as E = 0 (rate = A exactly, independent of T, so
#   the Dirichlet-pinned temperature field cannot perturb the kinetics);
# - depletion limiter DISABLED (tanh roll-off below ~1 kg/m³ would distort
#   the late decay: ξ_A < 1 kg/m³ for t > ln(ξ₀)/k₁ ≈ 17.8 s);
# - temperature gates left at defaults (fully open);
# - h = 0 on both reactions keeps the single cell exactly isothermal.

function _verification_case_vrx01()
    ξ0 = 1000.0
    k1, k2 = 0.389, 0.262
    T_iso = 400.0
    L = 1e-3
    ts = [1.0, 2.0, 4.0, 8.0, 12.0, 20.0]

    function build(n_cells)
        compA = SolidComponent(:A, ρ = ξ0, c = 1000.0, k = 0.5)
        compB = SolidComponent(:B, ρ = ξ0, c = 1000.0, k = 0.5)
        compC = SolidComponent(:C, ρ = ξ0, c = 1000.0, k = 0.5)
        r1 = Reaction(:r1, 1 => 2; A = k1, E = 0.0, h = 0.0)
        r2 = Reaction(:r2, 2 => 3; A = k2, E = 0.0, h = 0.0)
        mat = Material(
            name = :VRX01,
            components = (compA, compB, compC),
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
            ξ_initial = [z -> ξ0, z -> 0.0, z -> 0.0],
            tspan = (0.0, ts[end]),
            output_times = vcat(0.0, ts),
        )
    end

    function qois(sol)
        out = QoIResult[]
        for (species, idx) in ((:A, 1), (:B, 2), (:C, 3))
            err = maximum(
                abs(
                    sol.ξ[idx, 1, j+1] - series_chain_xi(t; ξ0 = ξ0, k1 = k1, k2 = k2)[idx],
                ) / ξ0 for (j, t) in enumerate(ts)
            )
            push!(
                out,
                QoIResult(
                    name = "max relative ξ_$(species) error",
                    value = err,
                    error = err,
                    tol = 1e-8,  # 3× observed max: 2.8e-9 (ξ_B; abstol=1e-12)
                    comment = "Bateman solution, k₁=0.389, k₂=0.262 (FDS/Gpyro values)",
                ),
            )
        end
        mass_err = maximum(abs(sum(sol.ξ[:, 1, j]) - ξ0) / ξ0 for j = 1:length(sol.t))
        push!(
            out,
            QoIResult(
                name = "total-mass conservation (closed chain)",
                value = mass_err,
                error = mass_err,
                tol = 1e-13,  # observed 2.3e-15
                comment = "Σξ must equal ξ₀; unit yields",
            ),
        )
        # Isothermality guard: with E = 0 nothing depends on T, but a drift
        # would signal an unexpected coupling (h = 0 ⇒ no reaction heat).
        T_err = maximum(abs(sol.T[1, j] - T_iso) for j = 1:length(sol.t))
        push!(
            out,
            QoIResult(
                name = "isothermality drift",
                value = T_err,
                error = T_err,
                tol = 1e-12,  # observed 5.7e-14
                comment = "h = 0 both reactions; Dirichlet-pinned single cell",
            ),
        )
        return out
    end

    return VerificationCase(
        id = "V-RX-01",
        title = "two-step series chain A→B→C, isothermal",
        tags = [:kinetics, :multi_step, :isothermal],
        references = [
            "FDS Verification Guide `two_step_solid_reaction`",
            "Lautenberger dissertation (2007) §3.4.4",
        ],
        build = build,
        solve_kwargs = (abstol = 1e-12, reltol = 1e-10),
        qois = qois,
        ci = CIConfig(n_cells = 1, budget_s = 10.0),
        report = (
            qoi_curves = (:xi_histories,),
            heating = :isothermal,
            curves = sol -> begin
                # FDS 6.11.0 overlay (reference/fds/custom/): the span-matched
                # re-run of two_step_solid_reaction (T_END = 20 s = our
                # window; identical kinetics 0.389/0.262 s⁻¹, E = 0,
                # isothermal — the same constants as Gpyro §3.4.4). FDS devc
                # values are time-averages over the output interval, so
                # samples plot at t − Δt_out/2; their Y (ρ_s with
                # ρ₀ = 1 kg/m³) scales to ξ by ξ0.
                fds = reference_columns(
                    "fds/custom/two_step_solid_reaction_matched_devc.csv",
                    "Time",
                    "Ya",
                    "Yb",
                    "Yc",
                )
                t_fds = fds[1][2:end] .- (fds[1][2] - fds[1][1]) / 2
                # Gpyro 0.8200 overlay (reference/gpyro/): dissertation §3.4.4
                # 0D isothermal deck, identical constants, dt = 1 ms (Gpyro's
                # species update is first-order in time).
                gp = reference_columns(
                    "gpyro/vrx01_result.csv",
                    "t_s",
                    "Y_A",
                    "Y_B",
                    "Y_C";
                    header_rows = 1,
                )
                [
                    ReportCurve(
                        name = :xi_histories,
                        title = "V-RX-01: series chain A→B→C (Bateman solution)",
                        xlabel = "t (s)",
                        ylabel = "ξ (kg/m³)",
                        series = vcat(
                            [
                                CurveSeries(
                                    label = "ξ_$(species)",
                                    x = sol.t,
                                    y = sol.ξ[idx, 1, :],
                                    exact = t ->
                                        series_chain_xi(t; ξ0 = ξ0, k1 = k1, k2 = k2)[idx],
                                ) for (species, idx) in ((:A, 1), (:B, 2), (:C, 3))
                            ],
                            [
                                CurveSeries(
                                    label = "FDS 6.11.0 ξ_$(species)",
                                    crosscode = true,
                                    x = t_fds,
                                    y = ξ0 .* fds[1+idx][2:end],
                                    exact = t ->
                                        series_chain_xi(t; ξ0 = ξ0, k1 = k1, k2 = k2)[idx],
                                ) for (species, idx) in ((:A, 1), (:B, 2), (:C, 3))
                            ],
                            [
                                CurveSeries(
                                    label = "Gpyro 0.8200 ξ_$(species)",
                                    crosscode = true,
                                    x = gp[1],
                                    y = ξ0 .* gp[1+idx],
                                    exact = t ->
                                        series_chain_xi(t; ξ0 = ξ0, k1 = k1, k2 = k2)[idx],
                                ) for (species, idx) in ((:A, 1), (:B, 2), (:C, 3))
                            ],
                        ),
                    ),
                ]
            end,
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vrx01)
