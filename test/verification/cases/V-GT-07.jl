# =============================================================================
# V-GT-07a/b — Transient binary diffusion through a mass-transfer film
# =============================================================================
#
# Adapted from Gpyro Verification and Validation Manual v1.0 (2026), §4.4,
# pp. 30–31, Eqs. (4.7)–(4.11): a 1 m porous slab initially containing
# carrier gas, exposed at its top to a reservoir with tracer fraction 0.22;
# the bottom is sealed. T = 300 K, P = 1 bar, φ = 0.6, D = 2e-5 m²/s.
# Members a/b use h_m = 1e-4 / 1e-2 kg/(m² s), respectively.
#
# Like the upstream input decks, both tracer gases have M = 0.032 kg/mol;
# equal molecular masses, heat capacities and wall coefficients preserve
# constant total pressure, density and temperature during counterdiffusion.
# LumpedTransfer with Λ = φD and κ = 0 realizes the constant-D equation;
# Chapman–Enskog diffusivity is replaced by the guide's specified constant.
# The zero-total-flow MassTransferFilm law reduces exactly to the Robin BC.
#
# At depth x = L-z, Y = Y∞ Φ(x,t), where Φ is the semi-infinite convective
# erfc solution with h/k = b = h_m/(φρ_gD). We stop at 400 s so the back
# face is over five diffusion lengths away; the reference tail is checked.
# Surface uptake is h_m Y∞ erfcx(b√(Dt)); its time integral independently
# predicts the stored tracer mass. All references use original inputs.
#
# Cell-center and cell-average profile errors are reported separately. The
# latter compares each finite-volume inventory with the analytical integral
# over its actual cell. Gpyro's boundary-node values retain point references.
# Resolution figures use the maximum over 50, 150 and 400 s; the existing
# spatial-order check uses the profile at 400 s.

function _vgt07_case(id, hm)
    L, T, P, φ, M = 1.0, 300.0, 1e5, 0.6, 0.032
    D, Y∞ = 2e-5, 0.22
    # Approximately 3× the n=160 errors for weak/strong surface exchange.
    tol_profile, tol_flux, tol_inventory =
        hm == 1e-4 ? (4e-4, 1.1e-3, 3e-3) : (4e-3, 7.5e-3, 8e-3)
    ρg = P*M/(Pyrolysis.Materials.R_GAS*T)
    b = hm/(φ*ρg*D)
    ts = [50.0, 150.0, 400.0]
    fraction(x, t) = Y∞*semiinf_convective_theta(x, t; h = b, k = 1.0, α = D)
    influx(t) = hm*Y∞*erfcx(b*sqrt(D*t))
    inventory(t) = φ*ρg*Y∞*robin_diffusion_uptake(t; D, b)
    function build(n)
        mat = Material(
            name = Symbol(id),
            components = (
                SolidComponent(:matrix; ρ = 400.0, φ, c = 1200.0, k = 0.2, λ = φ*D),
                GaseousComponent(:tracer; M, c = 1000.0, k = 0.03),
                GaseousComponent(:carrier; M, c = 1000.0, k = 0.03),
            ),
        )
        PyrolysisProblem(
            mesh = uniform_mesh(L, n),
            material = mat,
            bc_set = BoundaryConditionSet(
                thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
                mass = (
                    top = AmbientGasBC(
                        GasReservoir(
                            P = P,
                            T = T,
                            mole_fractions = (tracer = Y∞, carrier = 1-Y∞),
                        );
                        film = MassTransferFilm(hm),
                    ),
                    bottom = ImpermeableBC(),
                ),
            ),
            T_initial = z -> T,
            ξ_initial = [z -> 400.0, z -> 0.0, z -> φ*ρg],
            tspan = (0.0, ts[end]),
            output_times = vcat(0.0, ts),
        )
    end
    function profile_error(sol; last_only = false, cell_average = false)
        z = sol.cell_centers[:, end]
        depths = L .- first.(sol.problem.mesh.positions)
        js = last_only ? (length(sol.t):length(sol.t)) : (2:length(sol.t))
        maximum(
            abs(
                sol.ξ[2, i, j]/(φ*ρg) - (
                    cell_average ?
                    Y∞*robin_diffusion_cell_average(
                        depths[i+1],
                        depths[i],
                        sol.t[j];
                        D,
                        b,
                    ) : fraction(L-z[i], sol.t[j])
                ),
            )/Y∞ for j in js for i in eachindex(z)
        )
    end
    function qois(sol)
        top = get_boundary_exchange(sol, :top)
        bottom = get_boundary_exchange(sol, :bottom)
        flux_error =
            maximum(abs(-top.fluxes[1, j]/influx(sol.t[j])-1) for j = 2:length(sol.t))
        mass_error =
            maximum(abs(sum(sol.m[2, :, j])/inventory(sol.t[j])-1) for j = 2:length(sol.t))
        total_error =
            maximum(abs(sum(sol.m[2:3, :, j])/(φ*ρg*L)-1) for j in eachindex(sol.t))
        [
            QoIResult(
                name = "transient tracer profile / reservoir fraction",
                value = profile_error(sol),
                error = profile_error(sol),
                tol = tol_profile,
                comment = "n=160: observed 1.28e-4 / 1.28e-3 for h_m=1e-4 / 1e-2; ~3× margin",
            ),
            QoIResult(
                name = "transient tracer cell average / reservoir fraction",
                value = profile_error(sol; cell_average = true),
                error = profile_error(sol; cell_average = true),
                tol = tol_profile,
                comment = "analytical integral over each cell; same tolerance as the center-value metric",
            ),
            QoIResult(
                name = "surface uptake vs analytical flux",
                value = flux_error,
                error = flux_error,
                tol = tol_flux,
                comment = "h_m Y∞ erfcx(b√(Dt)); observed 3.34e-4 / 2.45e-3; ~3× margin",
            ),
            QoIResult(
                name = "stored tracer vs integrated analytical influx",
                value = mass_error,
                error = mass_error,
                tol = tol_inventory,
                comment = "integrated exact influx; observed 9.25e-4 / 2.55e-3; ~3× margin",
            ),
            QoIResult(
                name = "total gas inventory conservation",
                value = total_error,
                error = total_error,
                tol = 1e-10,
                comment = "equal-mass binary counterdiffusion has zero net mass influx",
            ),
            QoIResult(
                name = "zero net boundary mass flux",
                value = maximum(abs, vec(sum(top.fluxes; dims = 1))),
                error = maximum(abs, vec(sum(top.fluxes; dims = 1)))/(hm*Y∞),
                tol = 3e-8,
                comment = "two species counterdiffuse; flux scale h_m Y∞; enthalpy form observed 1.1e-8 at these tolerances, falling with them to 8.7e-12 at 1e-3 of them",
            ),
            QoIResult(
                name = "sealed back face",
                value = maximum(abs, bottom.fluxes),
                error = maximum(abs, bottom.fluxes),
                tol = 1e-14,
            ),
            QoIResult(
                name = "isobaric invariant",
                value = maximum(abs.(sol.gas.pressure .- P)),
                error = maximum(abs.(sol.gas.pressure .- P))/P,
                tol = 4e-9,
                comment = "observed relative drift ≤1.05e-9 (0.000105 Pa); pressure tolerance 0.001 Pa",
            ),
            QoIResult(
                name = "isothermal invariant",
                value = maximum(abs.(sol.T .- T)),
                error = maximum(abs.(sol.T .- T)),
                tol = 1e-8,
            ),
            QoIResult(
                name = "semi-infinite reference tail at back",
                value = fraction(L, ts[end])/Y∞,
                error = fraction(L, ts[end])/Y∞,
                tol = 1e-12,
                comment = "bounds the neglected far-boundary penetration",
            ),
        ]
    end
    kwargs = (reltol = 1e-9, abstol = 1e-11, pressure_tolerance = 0.001)
    VerificationCase(
        id = id,
        title = "transient binary diffusion through a mass film (h_m = $(hm))",
        tags = [:gas_transport, :diffusion, :mass_film, :transient, :analytical],
        references = [
            "[Gpyro Verification and Validation Manual v1.0 (2026), §4.4, pp. 30–31, Eqs. (4.7)–(4.11)]($(GPYRO_VERIFICATION_GUIDE)#page=32)",
        ],
        build = build,
        solve_kwargs = kwargs,
        qois = qois,
        ci = CIConfig(
            n_cells = 160,
            budget_s = 60.0,
            convergence = ConvergenceSpec(
                n_cells = [40, 80, 160],
                error_fn = sol -> profile_error(sol; last_only = true),
                expected_order = 2.0,
                min_order = 1.8,
                min_mean_order = 1.9,
                solve_kwargs = kwargs,
            ),
        ),
        report = (
            n_cells_ladder = [40, 80, 160],
            cost_n_cells = [320, 640],
            error_norms = sol -> (
                Linf = profile_error(sol; last_only = true),
                cell_average_Linf = profile_error(
                    sol;
                    last_only = true,
                    cell_average = true,
                ),
            ),
            resolution_curves = (ns, qs) -> _vgt07_resolution_curves(id, ns, qs; L),
            comparison_note = "Gpyro V2 (`bff0126c68f2`) uses 161 boundary-centered nodes " *
                              "(160 intervals) and Δt = 0.05 s. Molecular collision parameters " *
                              "give D = 2e-5 m²/s at 300 K and 1 bar. The far-side reservoir " *
                              "contains only carrier gas so Gpyro's positive floor on the zero " *
                              "film coefficient contributes negligible exchange. Its influx is " *
                              "the film flux h_m(Y∞ − Y_s) and its inventory integrates the " *
                              "nodes over half end cells at its native gas density. All codes " *
                              "are compared with the analytical Robin solution, with their " *
                              "native gas constants retained; Gpyro's are point values.",
            curves = sol -> begin
                half_width = L/(2length(sol.cell_centers[:, end]))
                [
                    ReportCurve(
                        name = :profiles,
                        title = "$(id): transient species diffusion",
                        xlabel = "depth (m)",
                        ylabel = "tracer mass fraction",
                        series = [
                            CurveSeries(
                                label = "t = $(sol.t[j]) s",
                                x = L .- sol.cell_centers[:, end],
                                y = vec(sol.ξ[2, :, j]) ./ (φ*ρg),
                                exact = x -> fraction(x, sol.t[j]),
                            ) for j = 2:length(sol.t)
                        ],
                    ),
                    ReportCurve(
                        name = :cell_average_profiles,
                        title = "$(id): tracer cell averages",
                        xlabel = "depth (m)",
                        ylabel = "cell-average tracer mass fraction",
                        series = [
                            CurveSeries(
                                label = "t = $(sol.t[j]) s",
                                x = L .- sol.cell_centers[:, end],
                                y = vec(sol.ξ[2, :, j]) ./ (φ*ρg),
                                exact = x ->
                                    Y∞*robin_diffusion_cell_average(
                                        max(0.0, x-half_width),
                                        min(L, x+half_width),
                                        sol.t[j];
                                        D,
                                        b,
                                    ),
                            ) for j = 2:length(sol.t)
                        ],
                    ),
                    ReportCurve(
                        name = :uptake,
                        title = "$(id): surface influx",
                        xlabel = "t (s)",
                        ylabel = "tracer influx (kg/m²/s)",
                        series = [
                            CurveSeries(
                                label = "surface influx",
                                x = sol.t[2:end],
                                y = -vec(get_boundary_exchange(sol, :top).fluxes[1, 2:end]),
                                exact = influx,
                            ),
                        ],
                    ),
                    ReportCurve(
                        name = :inventory,
                        title = "$(id): tracer inventory",
                        xlabel = "t (s)",
                        ylabel = "tracer mass / area (kg/m²)",
                        series = [
                            CurveSeries(
                                label = "stored tracer",
                                x = sol.t,
                                y = vec(sum(sol.m[2, :, :]; dims = 1)),
                                exact = inventory,
                            ),
                        ],
                    ),
                ]
            end,
        ),
    )
end

# Maximum error over the three output times against mesh spacing, for the
# center-value and cell-average profile references and the boundary ledgers.
function _vgt07_resolution_curves(id, ns, qs; L)
    errors(name) = [only(q.error for q in rung if q.name == name) for rung in qs]
    profiles = ReportCurve(
        name = :resolution_profiles,
        title = "$(id): profile error vs mesh spacing (all output times)",
        xlabel = "mesh spacing (m)",
        ylabel = "maximum error / reservoir fraction",
        xscale = :log10,
        yscale = :log10,
        series = [
            CurveSeries(
                label = "Pyrolysis.jl center values",
                x = L ./ ns,
                y = errors("transient tracer profile / reservoir fraction"),
            ),
            CurveSeries(
                label = "Pyrolysis.jl cell averages",
                x = L ./ ns,
                y = errors("transient tracer cell average / reservoir fraction"),
            ),
        ],
    )
    ledgers = [
        ReportCurve(
            name = slug,
            title = "$(id): $(title) error vs mesh spacing (all output times)",
            xlabel = "mesh spacing (m)",
            ylabel = "maximum relative error",
            xscale = :log10,
            yscale = :log10,
            series = [CurveSeries(label = "Pyrolysis.jl", x = L ./ ns, y = errors(qoi))],
        ) for (slug, title, qoi) in (
            (:resolution_uptake, "surface influx", "surface uptake vs analytical flux"),
            (
                :resolution_inventory,
                "tracer inventory",
                "stored tracer vs integrated analytical influx",
            ),
        )
    ]
    return [profiles; ledgers]
end
push!(VERIFICATION_CASE_BUILDERS, () -> _vgt07_case("V-GT-07a", 1e-4))
push!(VERIFICATION_CASE_BUILDERS, () -> _vgt07_case("V-GT-07b", 1e-2))
