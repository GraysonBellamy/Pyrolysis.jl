# =============================================================================
# V-GT-08a/b — Imposed gas inflow and a fixed-pressure Darcy outlet
# =============================================================================
#
# Adapted from Gpyro Verification and Validation Manual v1.0 (2026), §4.3,
# p. 29, Eqs. (4.5)–(4.6) and Appendix A.2: a 0.1 m isothermal porous
# slab with permeability 5e-10 m², T = 300 K, an outlet at 1 bar and an
# imposed incoming pure-gas flux of 10 kg/(m² s). Use constant dynamic
# viscosity μ = 2e-5 Pa s instead of Gpyro's Chapman–Enskog correlation.
#
# Darcy's law and the ideal-gas EOS give J = -(κM/μRT) P dP/dz, so
# P² = Pout² + 2μRT J x/(κM), where x is distance inward from the outlet.
# There is no reaction or diffusive wall transport. The initial pressure is
# uniform at 1 bar; the solve must relax to this independently derived profile.
# Member a feeds the bottom and vents the top; b mirrors both boundaries.
# MassFluxBC is outward-positive, so the inlet flux is -J in both cases.
# An adiabatic matrix, gas and incoming gas all at 300 K preserve T exactly.

function _vgt08_case(id; reverse = false)
    L, T, Pout, M, φ = 0.1, 300.0, 1e5, 0.029, 0.6
    κ, μ, J = 5e-10, 2e-5, 10.0
    depth(z) = reverse ? z : L-z
    pressure(z) = darcy_imposed_flux_pressure(depth(z); Pout, flux = J, μ, T, κ, M)
    ΔP² = pressure(reverse ? L : 0.0)^2-Pout^2
    inlet, outlet = reverse ? (:top, :bottom) : (:bottom, :top)
    function build(n)
        mat = Material(
            name = Symbol(id),
            gas_viscosity = μ,
            components = (
                SolidComponent(:matrix; ρ = 400.0, φ, c = 1200.0, k = 0.2, κ, λ = 0.0),
                GaseousComponent(:gas; M, c = 1000.0, k = 0.03),
            ),
        )
        feed = MassFluxBC(fluxes = (gas = -J,), inflow_temperature = (gas = T,))
        vent = AmbientGasBC(GasReservoir(P = Pout, T = T, mole_fractions = (gas = 1.0,)))
        PyrolysisProblem(
            mesh = uniform_mesh(L, n),
            material = mat,
            bc_set = BoundaryConditionSet(
                thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
                mass = reverse ? (top = feed, bottom = vent) : (top = vent, bottom = feed),
            ),
            T_initial = z -> T,
            ξ_initial = [z -> 400.0, z -> φ*Pout*M/(Pyrolysis.Materials.R_GAS*T)],
            tspan = (0.0, 1.0),
            output_times = [0.0, 0.5, 1.0],
        )
    end
    profile_error(sol) =
        maximum(
            abs.(sol.gas.pressure[:, end] .^ 2 .- pressure.(sol.cell_centers[:, end]) .^ 2),
        )/ΔP²
    function qois(sol)
        inflow = get_boundary_exchange(sol, inlet).fluxes[1, :]
        outflow = get_boundary_exchange(sol, outlet).fluxes[1, :]
        inventory_drift =
            abs(sum(sol.m[2, :, end])-sum(sol.m[2, :, end-1]))/sum(sol.m[2, :, end])
        [
            QoIResult(
                name = "linear P² profile / pressure-squared rise",
                value = profile_error(sol),
                error = profile_error(sol),
                tol = 2.5e-3,
                comment = "Gpyro Eq. (4.5), μRT/M in place of α; n=80 observed 8.13e-4, ~3× margin",
            ),
            QoIResult(
                name = "prescribed signed inlet flux",
                value = inflow[end],
                exact = -J,
                error = maximum(abs.(inflow .+ J))/J,
                tol = 1e-12,
            ),
            QoIResult(
                name = "outlet flux vs imposed inflow",
                value = outflow[end],
                exact = J,
                error = abs(outflow[end]/J-1),
                tol = 1e-8,
                comment = "independent steady mass balance after pressure relaxation",
            ),
            QoIResult(
                name = "steady signed boundary balance",
                value = outflow[end]+inflow[end],
                error = abs(outflow[end]+inflow[end])/J,
                tol = 1e-8,
            ),
            QoIResult(
                name = "gas inventory steadiness",
                value = inventory_drift,
                error = inventory_drift,
                tol = 1e-8,
                comment = "relative inventory change between 0.5 and 1 s",
            ),
            QoIResult(
                name = "isothermal invariant",
                value = maximum(abs.(sol.T .- T)),
                error = maximum(abs.(sol.T .- T)),
                tol = 1e-8,
            ),
        ]
    end
    kwargs = (reltol = 1e-9, abstol = 1e-11, pressure_tolerance = 0.001)
    VerificationCase(
        id = id,
        title = "imposed-inflow Darcy slab ($(reverse ? "top to bottom" : "bottom to top"))",
        tags = [:gas_transport, :darcy, :mass_flux_bc, :analytical],
        references = [
            "[Gpyro Verification and Validation Manual v1.0 (2026), §4.3, p. 29, Eqs. (4.5)–(4.6), Appendix A.2]($(GPYRO_VERIFICATION_GUIDE)#page=31)",
        ],
        build = build,
        solve_kwargs = kwargs,
        qois = qois,
        ci = CIConfig(
            n_cells = 80,
            budget_s = 60.0,
            convergence = ConvergenceSpec(
                n_cells = [20, 40, 80],
                error_fn = profile_error,
                expected_order = 1.0,
                min_order = 0.8,
                min_mean_order = 0.9,
                solve_kwargs = kwargs,
            ),
        ),
        report = (
            n_cells_ladder = [20, 40, 80],
            error_norms = sol -> (P2 = profile_error(sol),),
            note = "The steady mass-balance plot shows the saved samples at 0.5 and 1 s. " *
                   "Markers should lie on the reference lines: inlet −10, outlet +10, " *
                   "and their sum zero kg/(m² s). Flux is outward-positive, so negative " *
                   "means gas enters the slab. The initial sample at t = 0 is omitted " *
                   "from this steady comparison: its outlet flux is zero while pressure builds.",
            comparison_note = "Gpyro V2 (`bff0126c68f2`) uses 81 boundary-centered nodes " *
                              "(80 intervals) and Δt = 0.005 s, with temperature held at 300 K. " *
                              "Its molecular collision parameters match the Darcy coefficient " *
                              "νP = μRT/M used in this problem. Native gauge pressure is converted " *
                              "to absolute pressure, depth to z = L − depth, and boundary fluxes " *
                              "to outward-positive kg/(m² s). Both codes are compared with the " *
                              "same steady solution.",
            curves = sol -> begin
                [
                    ReportCurve(
                        name = :pressure,
                        title = "$(id): prescribed-inflow pressure",
                        xlabel = "z (m)",
                        ylabel = "P (Pa)",
                        series = [
                            CurveSeries(
                                label = "steady pressure",
                                x = sol.cell_centers[:, end],
                                y = sol.gas.pressure[:, end],
                                exact = pressure,
                            ),
                        ],
                    ),
                    ReportCurve(
                        name = :boundary_flux,
                        title = "$(id): steady mass balance",
                        xlabel = "t (s)",
                        ylabel = "outward mass flux (kg/m²/s)",
                        series = [
                            CurveSeries(
                                label = "inlet ($(inlet))",
                                x = sol.t[2:end],
                                y = vec(get_boundary_exchange(sol, inlet).fluxes[1, 2:end]),
                                exact = t -> -J,
                            ),
                            CurveSeries(
                                label = "outlet ($(outlet))",
                                x = sol.t[2:end],
                                y = vec(
                                    get_boundary_exchange(sol, outlet).fluxes[1, 2:end],
                                ),
                                exact = t -> J,
                            ),
                            CurveSeries(
                                label = "net (inlet + outlet)",
                                x = sol.t[2:end],
                                y = vec(
                                    get_boundary_exchange(sol, inlet).fluxes[1, 2:end] .+
                                    get_boundary_exchange(sol, outlet).fluxes[1, 2:end],
                                ),
                                exact = t -> 0.0,
                            ),
                        ],
                    ),
                ]
            end,
        ),
    )
end
push!(VERIFICATION_CASE_BUILDERS, () -> _vgt08_case("V-GT-08a"))
push!(VERIFICATION_CASE_BUILDERS, () -> _vgt08_case("V-GT-08b"; reverse = true))
