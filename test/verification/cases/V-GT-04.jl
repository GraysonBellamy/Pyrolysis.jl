# At constant total pressure, summing the two friction equations gives
# n₁/K₁+n₂/K₂=0. Substitution into species 1 gives
# n₁=-P x₁'/(RT [A-B x₁]), A=1/K₁+tP, B=tP(1-K₂/K₁).
# Thus A-Bx(z) is exponential and the flux follows by a logarithmic integral.
function _verification_case_vgt04()
    reference=setprecision(256) do
        P=BigFloat(1e5)
        T=BigFloat(400)
        R=BigFloat(Pyrolysis.Materials.R_GAS)
        φ=BigFloat(0.5)
        τ=BigFloat(2)
        dp=BigFloat(5e-7)
        L=BigFloat(0.001)
        K1=φ*dp/(3τ)*sqrt(8R*T/(big(π)*BigFloat(0.004)))
        K2=φ*dp/(3τ)*sqrt(8R*T/(big(π)*BigFloat(0.04)))
        resistance=τ/φ*(BigFloat(300)/T)^(BigFloat(7)/4)/(BigFloat(2e-5)*BigFloat(1e5))
        A=inv(K1)+resistance*P
        B=resistance*P*(1-K2/K1)
        xl, xr=BigFloat(0.8), BigFloat(0.2)
        logratio=log((A-B*xr)/(A-B*xl))
        n=P*logratio/(R*T*L*B)
        x=z->Float64((A-(A-B*xl)*exp(logratio*BigFloat(z)/L))/B)
        (; x, n = Float64(n), ratio = Float64(-K2/K1))
    end
    function build(n)
        mat=Material(
            name = :GrahamExchange,
            gas_transport = DustyGas(
                tortuosity = 2.0,
                D_ref = 2e-5,
                T_ref = 300.0,
                P_ref = 1e5,
                pore_diameter = 5e-7,
            ),
            components = (
                SolidComponent(:matrix; ρ = 500.0, φ = 0.5, c = 1000.0, k = 0.2, κ = 0.0),
                GaseousComponent(:light; M = 0.004, c = 0.0, k = 0.03),
                GaseousComponent(:heavy; M = 0.04, c = 0.0, k = 0.03),
            ),
        )
        PyrolysisProblem(
            mesh = create_uniform_mesh(0.001, n),
            material = mat,
            T_initial = z->400.0,
            ξ_initial = [
                z->500.0,
                z->0.5*1e5*0.004*0.5/(Pyrolysis.Materials.R_GAS*400.0),
                z->0.5*1e5*0.04*0.5/(Pyrolysis.Materials.R_GAS*400.0),
            ],
            bc_set = BoundaryConditionSet(
                thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
                mass = (
                    bottom = AmbientGasBC(
                        GasReservoir(
                            P = 1e5,
                            T = 400.0,
                            mole_fractions = (light = 0.8, heavy = 0.2),
                        ),
                    ),
                    top = AmbientGasBC(
                        GasReservoir(
                            P = 1e5,
                            T = 400.0,
                            mole_fractions = (light = 0.2, heavy = 0.8),
                        ),
                    ),
                ),
            ),
            tspan = (0.0, 0.5),
            output_times = [0.0, 0.25, 0.5],
        )
    end
    profile_error(sol) = maximum(
        abs.(get_gas_composition(sol)[1, :, end] .- reference.x.(get_cell_positions(sol))),
    )
    flux_error(sol) =
        abs(get_boundary_exchange(sol, :top).fluxes[1, end]/(0.004*reference.n)-1)
    function qois(sol)
        flux=get_boundary_exchange(sol, :top).fluxes[:, end] ./ [0.004, 0.04]
        [
            QoIResult(
                name = "Graham wall-weighted flux ratio",
                value = flux[2]/flux[1],
                exact = reference.ratio,
                error = abs(flux[2]/flux[1]-reference.ratio),
                tol = 1e-7,
            ),
            QoIResult(
                name = "logarithmic binary flux",
                value = flux[1],
                exact = reference.n,
                error = abs(flux[1]/reference.n-1),
                tol = 0.015,
            ),
            QoIResult(
                name = "exponential composition profile",
                value = profile_error(sol),
                error = profile_error(sol),
                tol = 0.005,
            ),
            QoIResult(
                name = "isobaric solution",
                value = maximum(abs.(get_pressure(sol) .- 1e5)),
                error = maximum(abs.(get_pressure(sol) .- 1e5)),
                tol = 0.01,
            ),
            QoIResult(
                name = "isothermal solution",
                value = maximum(abs.(sol.T .- 400.0)),
                error = maximum(abs.(sol.T .- 400.0)),
                tol = 1e-8,
            ),
        ]
    end
    VerificationCase(
        id = "V-GT-04",
        title = "isobaric dusty-gas binary exchange (Graham relation)",
        tags = [:gas_transport, :dusty_gas, :analytical],
        references = ["integrated binary friction equation at constant pressure"],
        build = build,
        solve_kwargs = (reltol = 1e-9, abstol = 1e-11, pressure_tolerance = 0.001),
        qois = qois,
        report = (
            n_cells_ladder = [20, 40, 80, 160],
            curves = sol->[
                ReportCurve(
                    name = :composition,
                    title = "Isobaric Graham exchange",
                    xlabel = "z (m)",
                    ylabel = "light-gas mole fraction",
                    series = [
                        CurveSeries(
                            label = "finite volume",
                            x = get_cell_positions(sol),
                            y = vec(get_gas_composition(sol)[1, :, end]),
                            exact = reference.x,
                        ),
                    ],
                ),
            ],
            error_norms = sol->(flux = flux_error(sol),),
        ),
        ci = CIConfig(
            n_cells = 80,
            budget_s = 60.0,
            convergence = ConvergenceSpec(
                n_cells = [20, 40, 80],
                error_fn = flux_error,
                expected_order = 1.0,
                min_order = 0.75,
                min_mean_order = 0.85,
                solve_kwargs = (
                    reltol = 1e-10,
                    abstol = 1e-12,
                    pressure_tolerance = 0.0001,
                ),
            ),
        ),
    )
end
push!(VERIFICATION_CASE_BUILDERS, _verification_case_vgt04)
