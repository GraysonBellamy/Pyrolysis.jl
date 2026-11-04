# The isobaric binary of V-GT-04 carrying three trace species. At constant
# pressure and steady total molar flux N, each species' friction equation is
# linear in its own pressure:
#   -p_j'/(RT) = n_j (1/K_j + tP) - t N p_j,
# so p_j(z) = p*_j + (p_j(0) - p*_j) G^(z/L), with G = exp(RTtNL) and
# p*_j = (p_j(L) - G p_j(0))/(1 - G), n_j = tN p*_j/(1/K_j + tP). Summing the
# species fluxes to N fixes G = Σ p_j(L)/(1+tPK_j) / Σ p_j(0)/(1+tPK_j) over
# every species, traces included. The second trace is the first scaled by
# 1e-12 and must stay proportional to it; the third is absent. Every present
# trace lies far above the pressure budget pressure_tolerance/Ng, below which
# the pore pressures guard the inventories. The friction equations also sum to
# Σ n_j/K_j = 0, the multicomponent Graham relation.
function _verification_case_vgt09()
    T, P, L=400.0, 1e5, 0.001
    masses=(light = 0.004, heavy = 0.04, trace = 0.02, scaled = 0.02, absent = 0.03)
    xa=(1e-6, 4e-6)
    scale=1e-12
    fractions(light, trace) = (
        light = light-trace-scale*trace,
        heavy = 1-light,
        trace = trace,
        scaled = scale*trace,
        absent = 0.0,
    )
    bottom, top=fractions(0.8, xa[1]), fractions(0.2, xa[2])
    reference=setprecision(256) do
        R=BigFloat(Pyrolysis.Materials.R_GAS)
        Tb, Pb, Lb=BigFloat(T), BigFloat(P), BigFloat(L)
        φ, τ, dp=BigFloat(0.5), BigFloat(2), BigFloat(5e-7)
        K=M->φ*dp/(3τ)*sqrt(8R*Tb/(big(π)*BigFloat(M)))
        t=τ/φ*(BigFloat(300)/Tb)^(BigFloat(7)/4)/(BigFloat(2e-5)*BigFloat(1e5))
        weight=x->sum(BigFloat(x[k])/(1+t*Pb*K(masses[k])) for k in keys(masses))
        G=weight(top)/weight(bottom)
        N=log(G)/(R*Tb*t*Lb)
        p0, pL=Pb*BigFloat(bottom.trace), Pb*BigFloat(top.trace)
        pstar=(pL-G*p0)/(1-G)
        x=z->Float64((pstar+(p0-pstar)*G^(BigFloat(z)/Lb))/Pb)
        (; x, n = Float64(t*N*pstar/(inv(K(masses.trace))+t*Pb)))
    end
    function build(n)
        gas(name) = GaseousComponent(name; M = getfield(masses, name), c = 0.0, k = 0.03)
        mat=Material(
            name = :TraceExchange,
            gas_transport = DustyGas(
                tortuosity = 2.0,
                D_ref = 2e-5,
                T_ref = 300.0,
                P_ref = 1e5,
                pore_diameter = 5e-7,
            ),
            components = (
                SolidComponent(:matrix; ρ = 500.0, φ = 0.5, c = 1000.0, k = 0.2, κ = 0.0),
                map(gas, keys(masses))...,
            ),
        )
        # Pore gas at the mean of the two reservoirs.
        initial=fractions(0.5, (xa[1]+xa[2])/2)
        ξ(name) =
            z->getfield(initial, name)*P*getfield(masses, name)*0.5/(
                Pyrolysis.Materials.R_GAS*T
            )
        PyrolysisProblem(
            mesh = uniform_mesh(L, n),
            material = mat,
            T_initial = z->T,
            ξ_initial = [z->500.0, map(ξ, keys(masses))...],
            bc_set = BoundaryConditionSet(
                thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
                mass = (
                    bottom = AmbientGasBC(
                        GasReservoir(P = P, T = T, mole_fractions = bottom),
                    ),
                    top = AmbientGasBC(GasReservoir(P = P, T = T, mole_fractions = top)),
                ),
            ),
            tspan = (0.0, 0.5),
            output_times = [0.0, 0.25, 0.5],
        )
    end
    # Species in the material's gas order; the matrix is component 1.
    moles(sol, j) = sol.ξ[1+j, :, end] ./ masses[j]
    fraction(sol, j) = moles(sol, j) ./ sum(moles(sol, k) for k = 1:length(masses))
    molar_flux(sol) = get_boundary_exchange(sol, :top).fluxes[:, end] ./ collect(masses)
    profile_error(sol) =
        maximum(abs.(fraction(sol, 3) .- reference.x.(sol.cell_centers[:, end]))/xa[2])
    flux_error(sol) = abs(molar_flux(sol)[3]/reference.n-1)
    function qois(sol)
        flux=molar_flux(sol)
        graham=flux .* sqrt.(collect(masses))
        proportional=maximum(abs.(sol.ξ[5, :, end] ./ (scale .* sol.ξ[4, :, end]) .- 1))
        [
            QoIResult(
                name = "trace flux",
                value = flux[3],
                exact = reference.n,
                error = flux_error(sol),
                tol = 0.005,
            ),
            QoIResult(
                name = "trace exponential profile",
                value = profile_error(sol),
                error = profile_error(sol),
                tol = 1e-13,
            ),
            QoIResult(
                name = "trace scaled by 1e-12 stays proportional",
                value = proportional,
                error = proportional,
                tol = 1e-13,
            ),
            QoIResult(
                name = "absent species stays absent",
                value = maximum(abs, sol.ξ[6, :, end]),
                error = maximum(abs, sol.ξ[6, :, end])+abs(flux[5]),
                tol = floatmin(Float64),
            ),
            QoIResult(
                name = "multicomponent Graham relation",
                value = sum(graham)/sum(abs, graham),
                error = abs(sum(graham))/sum(abs, graham),
                tol = 1e-12,
            ),
            QoIResult(
                name = "isobaric solution",
                value = maximum(abs.(sol.gas.pressure .- P)),
                error = maximum(abs.(sol.gas.pressure .- P)),
                tol = 1e-9,
            ),
        ]
    end
    VerificationCase(
        id = "V-GT-09",
        title = "isobaric dusty-gas exchange carrying trace species",
        tags = [:gas_transport, :dusty_gas, :analytical],
        references = ["linear friction equations of an isobaric multicomponent mixture"],
        build = build,
        solve_kwargs = (reltol = 1e-9, abstol = 1e-11, pressure_tolerance = 1e-15),
        qois = qois,
        report = (
            n_cells_ladder = [20, 40, 80, 160],
            curves = sol->[
                ReportCurve(
                    name = :composition,
                    title = "Trace species in an isobaric Graham exchange",
                    xlabel = "z (m)",
                    ylabel = "trace mole fraction",
                    series = [
                        CurveSeries(
                            label = "finite volume",
                            x = sol.cell_centers[:, end],
                            y = fraction(sol, 3),
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
                solve_kwargs = (reltol = 1e-10, abstol = 1e-12, pressure_tolerance = 1e-15),
            ),
        ),
    )
end
push!(VERIFICATION_CASE_BUILDERS, _verification_case_vgt09)
