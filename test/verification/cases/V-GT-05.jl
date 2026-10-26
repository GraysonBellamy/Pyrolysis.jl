# A continuously varying temperature and wall coefficient give an independent
# integrated resistance: n=(PL-PR)/(R ∫T/Λ(T) dz). With Λ=Λ0(T/T0)^2
# and linear T, the integral is T0² L log(TR/TL)/(Λ0 (TR-TL)).
function _verification_case_vgt05()
    reference=setprecision(256) do
        TL, TR=BigFloat(500), BigFloat(900)
        PL, PR=BigFloat(2e5), BigFloat(1e5)
        Λ0=BigFloat(1e-5)
        L=BigFloat(0.002)
        R=BigFloat(Pyrolysis.Materials.R_GAS)
        resistance=TL*TL*L*log(TR/TL)/(Λ0*(TR-TL))
        n=(PL-PR)/(R*resistance)
        T=z->Float64(TL+(TR-TL)*BigFloat(z)/L)
        p=z->Float64(PL-(PL-PR)*log((TL+(TR-TL)*BigFloat(z)/L)/TL)/log(TR/TL))
        (; T, p, flux = Float64(BigFloat(0.03)*n))
    end
    function build(n)
        mat=Material(
            name = :NonisothermalResistance,
            components = (
                SolidComponent(
                    :matrix;
                    ρ = 500.0,
                    φ = 0.4,
                    c = 1000.0,
                    k = 0.2,
                    λ = T->1e-5*(T/500.0)^2,
                ),
                GaseousComponent(:gas; M = 0.03, c = 0.0, k = 0.03),
            ),
        )
        PyrolysisProblem(
            mesh = create_uniform_mesh(0.002, n),
            material = mat,
            T_initial = reference.T,
            ξ_initial = [
                z->500.0,
                z->0.4*(2e5-1e5*z/0.002)*0.03/(Pyrolysis.Materials.R_GAS*reference.T(z)),
            ],
            bc_set = BoundaryConditionSet(
                thermal = (
                    bottom = DirichletThermalBC(500.0),
                    top = DirichletThermalBC(900.0),
                ),
                mass = (
                    bottom = AmbientGasBC(
                        GasReservoir(P = 2e5, T = 500.0, mole_fractions = (gas = 1.0,)),
                    ),
                    top = AmbientGasBC(
                        GasReservoir(P = 1e5, T = 900.0, mole_fractions = (gas = 1.0,)),
                    ),
                ),
            ),
            tspan = (0.0, 2.0),
            output_times = [0.0, 1.0, 2.0],
        )
    end
    flux_error(sol) = abs(get_boundary_exchange(sol, :top).fluxes[1, end]/reference.flux-1)
    profile_error(sol) =
        maximum(abs.(get_pressure(sol)[:, end] .- reference.p.(get_cell_positions(sol))))/1e5
    qois=sol->[
        QoIResult(
            name = "integrated nonisothermal wall resistance",
            value = flux_error(sol),
            error = flux_error(sol),
            tol = 1e-4,
        ),
        QoIResult(
            name = "logarithmic pressure profile",
            value = profile_error(sol),
            error = profile_error(sol),
            tol = 1e-4,
        ),
        QoIResult(
            name = "linear thermal field",
            value = maximum(abs.(sol.T[:, end] .- reference.T.(get_cell_positions(sol)))),
            error = maximum(abs.(sol.T[:, end] .- reference.T.(get_cell_positions(sol)))),
            tol = 1e-8,
        ),
    ]
    VerificationCase(
        id = "V-GT-05",
        title = "nonisothermal integrated wall resistance",
        tags = [:gas_transport, :analytical, :nonisothermal],
        references = ["direct integral of T/Λ(T) along a linear temperature field"],
        build = build,
        solve_kwargs = (reltol = 1e-9, abstol = 1e-11, pressure_tolerance = 0.001),
        qois = qois,
        report = (
            n_cells_ladder = [20, 40, 80, 160],
            curves = sol->[
                ReportCurve(
                    name = :pressure,
                    title = "Nonisothermal integrated resistance",
                    xlabel = "z (m)",
                    ylabel = "P (Pa)",
                    series = [
                        CurveSeries(
                            label = "finite volume",
                            x = get_cell_positions(sol),
                            y = get_pressure(sol)[:, end],
                            exact = reference.p,
                        ),
                    ],
                ),
            ],
            error_norms = sol->(flux = flux_error(sol),),
        ),
        ci = CIConfig(
            n_cells = 80,
            budget_s = 30.0,
            convergence = ConvergenceSpec(
                n_cells = [20, 40, 80],
                error_fn = flux_error,
                expected_order = 2.0,
                min_order = 1.8,
                min_mean_order = 1.9,
                solve_kwargs = (
                    reltol = 1e-10,
                    abstol = 1e-12,
                    pressure_tolerance = 0.0001,
                ),
            ),
        ),
    )
end
push!(VERIFICATION_CASE_BUILDERS, _verification_case_vgt05)
