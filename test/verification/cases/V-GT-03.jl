# Pure-gas slab: integrate (Λ + κ P/μ) dP = -n RT dz in each layer.
# A quadratic potential gives the profile and a scalar monotone interface
# equation fixes the transmitted flux. All reference coefficients are built
# from the original material inputs inside 256-bit arithmetic.
function _vgt03_reference(; dusty = false, layered = false, reverse = false)
    setprecision(256) do
        T=BigFloat(500)
        R=BigFloat(Pyrolysis.Materials.R_GAS)
        M=BigFloat(0.03)
        L=BigFloat(0.002)
        μ=BigFloat(2e-5)
        PL, PR=reverse ? (BigFloat(1e5), BigFloat(2e5)) : (BigFloat(2e5), BigFloat(1e5))
        φl, φr=BigFloat(0.4), BigFloat(0.7)
        λl = dusty ? φl*BigFloat(1e-6)/6*sqrt(8R*T/(big(π)*M)) : BigFloat(1e-5)
        λr = dusty ? φr*BigFloat(1e-6)/6*sqrt(8R*T/(big(π)*M)) : layered ? zero(T) : λl
        ml=layered ? zero(T) : BigFloat(5e-15)/μ
        mr=layered ? BigFloat(5e-15)/μ : ml
        S(p, λ, m) = λ*p+m*p*p/2
        balance=p->(S(PL, λl, ml)-S(p, λl, ml))-(S(p, λr, mr)-S(PR, λr, mr))
        lo, hi=minmax(PL, PR)
        for _ = 1:260
            p=(lo+hi)/2
            balance(p)>0 ? (lo=p) : (hi=p)
        end
        Pi=(lo+hi)/2
        Q=2*(S(PL, λl, ml)-S(Pi, λl, ml))/L
        inverse(S0, λ, m) = 2S0/(λ+sqrt(λ*λ+2m*S0))
        profile=z->Float64(
            z<=L/2 ? inverse(S(PL, λl, ml)-Q*BigFloat(z), λl, ml) :
            inverse(S(Pi, λr, mr)-Q*(BigFloat(z)-L/2), λr, mr),
        )
        return (;
            profile,
            flux = Float64(M*Q/(R*T)),
            PL = Float64(PL),
            PR = Float64(PR),
            T = Float64(T),
            L = Float64(L),
        )
    end
end

function _vgt03_case(id; dusty = false, layered = false, reverse = false)
    ref=_vgt03_reference(; dusty, layered, reverse)
    function build(n)
        model=dusty ?
              DustyGas(
            tortuosity = 2.0,
            D_ref = 2e-5,
            T_ref = 300.0,
            P_ref = 1e5,
            pore_diameter = 1e-6,
        ) : LumpedTransfer()
        material=Material(
            name = Symbol(id),
            gas_transport = model,
            gas_viscosity = 2e-5,
            components = (
                SolidComponent(
                    :left;
                    ρ = 500.0,
                    φ = 0.4,
                    c = 1000.0,
                    k = 0.2,
                    λ = 1e-5,
                    κ = layered ? 0.0 : 5e-15,
                ),
                SolidComponent(
                    :right;
                    ρ = 500.0,
                    φ = 0.7,
                    c = 1000.0,
                    k = 0.2,
                    λ = layered ? 0.0 : 1e-5,
                    κ = 5e-15,
                ),
                GaseousComponent(:gas; M = 0.03, c = 0.0, k = 0.03),
            ),
        )
        left=GasReservoir(P = ref.PL, T = ref.T, mole_fractions = (gas = 1.0,))
        right=GasReservoir(P = ref.PR, T = ref.T, mole_fractions = (gas = 1.0,))
        PyrolysisProblem(
            mesh = create_uniform_mesh(ref.L, n),
            material = material,
            bc_set = BoundaryConditionSet(
                thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
                mass = (bottom = AmbientGasBC(left), top = AmbientGasBC(right)),
            ),
            T_initial = z->ref.T,
            ξ_initial = [
                z->z<ref.L/2 ? 500.0 : 0.0,
                z->z<ref.L/2 ? 0.0 : 500.0,
                z->(z<ref.L/2 ? 0.4 : 0.7)*ref.profile(z)*0.03/(
                    Pyrolysis.Materials.R_GAS*ref.T
                ),
            ],
            tspan = (0.0, 2.0),
            output_times = [0.0, 1.0, 2.0],
        )
    end
    pressure_error(sol) =
        maximum(abs.(get_pressure(sol)[:, end] .- ref.profile.(get_cell_positions(sol))))/abs(
            ref.PL-ref.PR,
        )
    function qois(sol)
        top=get_boundary_exchange(sol, :top)
        bottom=get_boundary_exchange(sol, :bottom)
        flux=top.fluxes[1, end]
        [
            QoIResult(
                name = "continuum pressure profile",
                value = pressure_error(sol),
                error = pressure_error(sol),
                tol = 0.01,
            ),
            QoIResult(
                name = "integrated continuum species flux",
                value = flux,
                exact = ref.flux,
                error = abs(flux/ref.flux-1),
                tol = 0.015,
            ),
            QoIResult(
                name = "steady signed boundary balance",
                value = flux+bottom.fluxes[1, end],
                error = abs((flux+bottom.fluxes[1, end])/ref.flux),
                tol = 1e-7,
            ),
            QoIResult(
                name = "isothermal invariant",
                value = maximum(abs.(sol.T .- ref.T)),
                error = maximum(abs.(sol.T .- ref.T)),
                tol = 1e-8,
            ),
        ]
    end
    VerificationCase(
        id = id,
        title = "pure-gas integrated transport ($(dusty ? "Knudsen + viscous" : layered ? "wall-only / Darcy-only layers" : "lumped + Darcy"), $(reverse ? "reverse" : "forward"))",
        tags = [:gas_transport, :analytical, :interface],
        references = ["integrated pressure potential ΛP + κP²/(2μ)"],
        build = build,
        solve_kwargs = (reltol = 1e-9, abstol = 1e-11, pressure_tolerance = 0.001),
        qois = qois,
        report = (
            n_cells_ladder = [20, 40, 80, 160],
            curves = sol->[
                ReportCurve(
                    name = :pressure,
                    title = "Integrated wall and Darcy pressure",
                    xlabel = "z (m)",
                    ylabel = "P (Pa)",
                    series = [
                        CurveSeries(
                            label = "finite volume",
                            x = get_cell_positions(sol),
                            y = get_pressure(sol)[:, end],
                            exact = ref.profile,
                        ),
                    ],
                ),
            ],
            error_norms = sol->(pressure = pressure_error(sol),),
        ),
        ci = CIConfig(
            n_cells = 80,
            budget_s = 30.0,
            convergence = ConvergenceSpec(
                n_cells = [20, 40, 80],
                error_fn = pressure_error,
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
_verification_case_vgt03() = _vgt03_case("V-GT-03")
_verification_case_vgt03b() = _vgt03_case("V-GT-03b"; layered = true)
_verification_case_vgt03c() = _vgt03_case("V-GT-03c"; dusty = true)
_verification_case_vgt03d() = _vgt03_case("V-GT-03d"; layered = true, reverse = true)
push!(
    VERIFICATION_CASE_BUILDERS,
    _verification_case_vgt03,
    _verification_case_vgt03b,
    _verification_case_vgt03c,
    _verification_case_vgt03d,
)
