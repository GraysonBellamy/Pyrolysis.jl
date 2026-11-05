# =============================================================================
# V-RX-05a–c — Insulated reaction with unequal species heat capacities
# =============================================================================
#
# Inspired by Gpyro Verification and Validation Manual v1.0 (2026), §3.6,
# p. 23, Eq. (3.11): endothermic (+1e6 J/kg), exothermic (-1e5 J/kg),
# and zero-reaction-heat cases, the latter producing 60% gas by mass.
#
# Adaptation: a sealed, insulated single cell starts at T₀ = 1000 K. Set
# E = 0 and A = 0.2 s⁻¹ so the first-order reactant mass fraction is exactly
# a(t) = exp(-0.2t), independently of the evolving temperature. Unlike the
# guide's source-term reconstruction from numerical reaction rates, both
# the species histories and temperature here have independent references.
#
# c_A = 1200 and c_B = 1900 J/(kg K); retained gas has c_G = 1000 J/(kg K).
# For solid yield ν, c_products = ν c_B + (1-ν)c_G. The temperature obeys
# [c_A a + c_products(1-a)] dT/dt = -h A a, with the closed-form logarithmic
# solution in exact/gpyro_verification.jl. The sensible-enthalpy reference is
# `H = m₀[c_A a + c_products(1-a)](T-T_REF)`, evaluated without solver rates.
#
# Cases a/b convert A entirely to B; c produces 0.4 B + 0.6 G with h = 0.
# Case c must keep T constant even though its sensible enthalpy changes:
# the H equation includes formation terms Σ h_j dm_j/dt. All products stay
# inside the cell, the depletion limiter is disabled, and no mesh moves.

function _vrx05_case(id, heat; solid_yield = 1.0)
    T0, rate, ξ0, L = 1000.0, 0.2, 1000.0, 0.001
    ca, cb, cg = 1200.0, 1900.0, 1000.0
    cp = solid_yield*cb + (1-solid_yield)*cg
    a(t) = exp(-rate*t)
    temperature(t) = adiabatic_reaction_temperature(t; T0, rate, heat, ca, cb = cp)
    ξexact(t) = (ξ0*a(t), ξ0*solid_yield*(1-a(t)), ξ0*(1-solid_yield)*(1-a(t)))
    H(t) = ξ0*L*(ca*a(t)+cp*(1-a(t)))*(temperature(t)-Pyrolysis.Materials.T_REF)
    function build(n)
        mat = Material(
            name = Symbol(id),
            depletion_limiter = DepletionLimiter(enabled = false),
            components = (
                SolidComponent(:A; ρ = ξ0, φ = 0.3, c = ca, k = 0.2, λ = 0.0),
                SolidComponent(:B; ρ = ξ0, φ = 0.3, c = cb, k = 0.2, λ = 0.0),
                GaseousComponent(:G; M = 0.032, c = cg, k = 0.03),
            ),
            reactions = (
                Reaction(
                    :convert,
                    1 => (2, 3);
                    A = rate,
                    E = 0.0,
                    h = heat,
                    yields = (solid_yield, 1-solid_yield),
                ),
            ),
        )
        PyrolysisProblem(
            mesh = uniform_mesh(L, n),
            material = mat,
            bc_set = BoundaryConditionSet(
                thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            ),
            T_initial = z -> T0,
            ξ_initial = [z -> ξ0, z -> 0.0, z -> 0.0],
            tspan = (0.0, 20.0),
            output_times = collect(0.0:0.5:20.0),
        )
    end
    temperature_error(sol) = maximum(
        abs(sol.T[i, j]-temperature(t)) for (j, t) in enumerate(sol.t) for
        i in axes(sol.T, 1)
    )
    function qois(sol)
        species_error = maximum(
            abs(sol.ξ[k, i, j]-ξexact(t)[k])/ξ0 for (j, t) in enumerate(sol.t) for
            i in axes(sol.T, 1) for k = 1:3
        )
        mass_error = maximum(abs(sum(sol.m[:, :, j])/(ξ0*L)-1) for j in eachindex(sol.t))
        enthalpy_error = maximum(
            abs(sol.extras.live_energy[j]-H(t))/abs(H(0.0)) for (j, t) in enumerate(sol.t)
        )
        [
            QoIResult(
                name = "independent species histories",
                value = species_error,
                error = species_error,
                tol = 1e-10,
                comment = "exact first-order conversion; observed ≤3e-12 in temperature form, integration margin",
            ),
            QoIResult(
                name = "adiabatic temperature history",
                value = temperature_error(sol),
                error = temperature_error(sol),
                tol = 1e-7,
                comment = "logarithmic exact solution; observed ≤3.17e-8 K in temperature form, ~3× margin",
            ),
            QoIResult(
                name = "sensible enthalpy history / initial energy",
                value = enthalpy_error,
                error = enthalpy_error,
                tol = 2e-10,
                comment = "independent Σ m_j c_j (T-T_REF); observed ≤6.19e-11 relative, ~3× margin",
            ),
            QoIResult(
                name = "closed total mass",
                value = mass_error,
                error = mass_error,
                tol = 1e-12,
            ),
            QoIResult(
                name = "sealed boundary mass flux",
                value = maximum(abs, sol.mass_loss_rate),
                error = maximum(abs, sol.mass_loss_rate),
                tol = 1e-14,
            ),
        ]
    end
    VerificationCase(
        id = id,
        title = "insulated reaction ($(heat > 0 ? "endothermic" : heat < 0 ? "exothermic" : "zero heat, retained gas"))",
        tags = [:kinetics, :reaction_enthalpy, :adiabatic, :analytical],
        references = [
            "[Gpyro Verification and Validation Manual v1.0 (2026), §3.6, p. 23, Eq. (3.11)]($(GPYRO_VERIFICATION_GUIDE)#page=25)",
            "Independent closed-form adaptation: E = 0, insulated and sealed cell, unequal heat capacities; see problem statement.",
        ],
        build = build,
        solve_kwargs = (reltol = 1e-10, abstol = 1e-12),
        qois = qois,
        ci = CIConfig(n_cells = 1, budget_s = 15.0),
        report = (
            comparison_note = "Gpyro V2 (`bff0126c68f2`) uses two boundary nodes on a " *
                              "uniform insulated slab, Δt = 0.001 s, E = 0, and the same " *
                              "reaction heat, yields and heat capacities. A 1e-6 Pa initial " *
                              "gas pressure avoids a singular vacuum diffusivity; its initial " *
                              "mass is negligible (about 1e-15 kg/m²). Sensible enthalpy is " *
                              "reconstructed from native concentrations and temperature at " *
                              "the common 298.15 K datum. " *
                              (
                                  solid_yield < 1 ?
                                  "The retained-gas case has a Gpyro mass discrepancy of about " *
                                  "3.27e-5 of initial mass that persists when Δt is halved; its " *
                                  "density and pressure routines use different gas constants. " :
                                  ""
                              ) *
                              "ThermaKin2Ds 7 uses one insulated element at dt = 0.001 s; " *
                              "its positive reaction heat is exothermic, so the sign is " *
                              "reversed" *
                              (
                                  solid_yield < 1 ?
                                  ", and its gas product stays in the element without " *
                                  "occupying volume." : "."
                              ),
            curves = sol -> begin
                [
                    ReportCurve(
                        name = :temperature,
                        title = "$(id): adiabatic temperature",
                        xlabel = "t (s)",
                        ylabel = "T (K)",
                        series = [
                            CurveSeries(
                                label = "temperature",
                                x = sol.t,
                                y = vec(sol.T[1, :]),
                                exact = temperature,
                            ),
                        ],
                    ),
                    ReportCurve(
                        name = :species,
                        title = "$(id): species histories",
                        xlabel = "t (s)",
                        ylabel = "mass concentration (kg/m³)",
                        series = [
                            CurveSeries(
                                label = string(name),
                                x = sol.t,
                                y = vec(sol.ξ[k, 1, :]),
                                exact = t -> ξexact(t)[k],
                            ) for (k, name) in enumerate((:A, :B, :G))
                        ],
                    ),
                    ReportCurve(
                        name = :enthalpy,
                        title = "$(id): sensible enthalpy",
                        xlabel = "t (s)",
                        ylabel = "H (J)",
                        series = [
                            CurveSeries(
                                label = "enthalpy",
                                x = sol.t,
                                y = sol.extras.live_energy,
                                exact = H,
                            ),
                        ],
                    ),
                ]
            end,
        ),
    )
end
push!(VERIFICATION_CASE_BUILDERS, () -> _vrx05_case("V-RX-05a", 1e6))
push!(VERIFICATION_CASE_BUILDERS, () -> _vrx05_case("V-RX-05b", -1e5))
push!(VERIFICATION_CASE_BUILDERS, () -> _vrx05_case("V-RX-05c", 0.0; solid_yield = 0.4))
