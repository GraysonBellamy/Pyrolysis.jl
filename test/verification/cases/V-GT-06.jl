# Coupled steady gas/energy reference, derived directly from the continuum:
# N'=0, (N cp T-k T')'=0, p'=-N R T/(M Λ).
# At the top, k T'=h β(N cp/h)(Tin-Ts). The resulting exponential T
# profile, integrated pore resistance and film equation determine N and Ts.
function _verification_case_vgt06()
    reference=setprecision(256) do
        TL, Tin=BigFloat(900), BigFloat(300)
        PL, PR=BigFloat(2e5), BigFloat(1e5)
        M, Λ, k, h, cp, L=BigFloat(0.03),
        BigFloat(0.001),
        BigFloat(1),
        BigFloat(200),
        BigFloat(1000),
        BigFloat(0.01)
        R=BigFloat(Pyrolysis.Materials.R_GAS)
        function thermal(N)
            C=N*cp
            b=C/k
            K=C/(-expm1(-b*L))
            hf=C/expm1(C/h)
            Ts=(K*TL+hf*Tin)/(K+hf)
            integral=TL*L+(Ts-TL)*(inv(b)-L/expm1(b*L))
            (; b, Ts, integral, hf)
        end
        lo, hi=BigFloat("1e-12"), BigFloat(1)
        for _ = 1:300
            mid=(lo+hi)/2
            mid*R*thermal(mid).integral/(M*Λ)>PL-PR ? (hi=mid) : (lo=mid)
        end
        N=(lo+hi)/2
        th=thermal(N)
        temperature=z->Float64(TL+(th.Ts-TL)*expm1(th.b*BigFloat(z))/expm1(th.b*L))
        integral=z->TL*BigFloat(z)+(th.Ts-TL)*(expm1(th.b*BigFloat(z))/th.b-BigFloat(z))/expm1(
            th.b*L,
        )
        pressure=z->Float64(PL-N*R*integral(z)/(M*Λ))
        (;
            flux = Float64(N),
            Ts = Float64(th.Ts),
            temperature,
            pressure,
            heat = Float64(th.hf*(Tin-th.Ts)),
        )
    end
    function build(n)
        mat=Material(
            name = :CoupledNonisothermalFilm,
            components = (
                SolidComponent(:matrix; ρ = 500.0, φ = 0.4, c = 1000.0, k = 1.0, λ = 0.001),
                GaseousComponent(:gas; M = 0.03, c = 1000.0, k = 0.03),
            ),
        )
        PyrolysisProblem(
            mesh = create_uniform_mesh(0.01, n),
            material = mat,
            T_initial = reference.temperature,
            ξ_initial = [
                z->500.0,
                z->0.4*reference.pressure(z)*0.03/(
                    Pyrolysis.Materials.R_GAS*reference.temperature(z)
                ),
            ],
            bc_set = BoundaryConditionSet(
                thermal = (
                    bottom = DirichletThermalBC(900.0),
                    top = ConvectiveBC(h = 200.0, T_inf = 300.0, blowing = CouetteFilm()),
                ),
                mass = (
                    bottom = AmbientGasBC(
                        GasReservoir(P = 2e5, T = 900.0, mole_fractions = (gas = 1.0,)),
                    ),
                    top = AmbientGasBC(
                        GasReservoir(P = 1e5, T = 300.0, mole_fractions = (gas = 1.0,));
                        film = MassTransferFilm(0.01),
                    ),
                ),
            ),
            tspan = (0.0, 300.0),
            output_times = [0.0, 150.0, 300.0],
        )
    end
    flux_error(sol) = abs(get_boundary_exchange(sol, :top).fluxes[1, end]/reference.flux-1)
    function qois(sol)
        b=get_boundary_exchange(sol, :top)
        [
            QoIResult(
                name = "coupled continuum flux",
                value = flux_error(sol),
                error = flux_error(sol),
                tol = 0.001,
            ),
            QoIResult(
                name = "coupled surface temperature",
                value = sol.surface_T[end],
                exact = reference.Ts,
                error = abs(sol.surface_T[end]-reference.Ts),
                tol = 0.2,
            ),
            QoIResult(
                name = "exponential thermal field",
                value = maximum(
                    abs.(sol.T[:, end] .- reference.temperature.(get_cell_positions(sol))),
                ),
                error = maximum(
                    abs.(sol.T[:, end] .- reference.temperature.(get_cell_positions(sol))),
                ),
                tol = 0.2,
            ),
            QoIResult(
                name = "continuum reservoir heat",
                value = b.heat[end],
                exact = reference.heat,
                error = abs(b.heat[end]/reference.heat-1),
                tol = 0.001,
            ),
            QoIResult(
                name = "coupled local root residual",
                value = maximum(b.residual),
                error = maximum(b.residual),
                tol = 1e-9,
            ),
        ]
    end
    VerificationCase(
        id = "V-GT-06",
        title = "coupled nonisothermal gas and thermal film",
        tags = [:gas_transport, :analytical, :nonisothermal, :film],
        references = ["integrated steady mass and energy equations with a Couette film"],
        build = build,
        solve_kwargs = (reltol = 1e-9, abstol = 1e-10, pressure_tolerance = 0.001),
        qois = qois,
        report = (
            n_cells_ladder = [10, 20, 40, 80],
            curves = sol->[
                ReportCurve(
                    name = :temperature,
                    title = "Coupled mass and thermal film",
                    xlabel = "z (m)",
                    ylabel = "T (K)",
                    series = [
                        CurveSeries(
                            label = "finite volume",
                            x = get_cell_positions(sol),
                            y = sol.T[:, end],
                            exact = reference.temperature,
                        ),
                    ],
                ),
                ReportCurve(
                    name = :pressure,
                    title = "Coupled nonisothermal pore pressure",
                    xlabel = "z (m)",
                    ylabel = "P (Pa)",
                    series = [
                        CurveSeries(
                            label = "finite volume",
                            x = get_cell_positions(sol),
                            y = get_pressure(sol)[:, end],
                            exact = reference.pressure,
                        ),
                    ],
                ),
            ],
            error_norms = sol->(flux = flux_error(sol),),
        ),
        ci = CIConfig(
            n_cells = 40,
            budget_s = 60.0,
            convergence = ConvergenceSpec(
                n_cells = [10, 20, 40],
                error_fn = flux_error,
                expected_order = 2.0,
                min_order = 1.7,
                min_mean_order = 1.8,
                solve_kwargs = (
                    reltol = 1e-10,
                    abstol = 1e-11,
                    pressure_tolerance = 0.0001,
                ),
            ),
        ),
    )
end
push!(VERIFICATION_CASE_BUILDERS, _verification_case_vgt06)
