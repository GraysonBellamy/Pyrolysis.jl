# =============================================================================
# V-HC-04 — Two-layer slab: steady piecewise-linear profile + interface flux
# =============================================================================
#
# FDS `insulated_steel_plate` / `ht_pile` analog. Two inert solids layered via
# composition step ICs; Dirichlet T at both faces; steady exact solution is
# the series-resistance network:
#
#   q = ΔT / (d₁/k₁ + d₂/k₂),   T_interface = T_bot + q·d₁/k₁,  piecewise linear
#
# What is verified: the distance-weighted harmonic-mean face conductivity is
# NODALLY EXACT for piecewise-constant k with the jump on a face — the steady
# cell-center temperatures and the boundary fluxes must match the resistance
# network to integrator precision, for k-ratios 10:1 and 100:1.
#
# Off-face variant (c): composition step ICs are SAMPLED AT CELL CENTERS, so
# a prescribed interface that is not on a face quantizes to the nearest face —
# there are no mixed cells. Measured: the numerics match the
# QUANTIZED-interface network to 1e-11 (i.e. the scheme is still nodally
# exact), and the entire user-facing error is the true-vs-quantized network
# gap (1.64% in flux, 1.47 K in profile at Δz = 0.5 mm) — pinned here against
# its a-priori Δz/2 sensitivity bound.
#
# No convergence spec: (a)/(b) are nodally exact (error is integrator-, not
# mesh-limited); (c) is quantization-dominated and non-monotone in n.

function _vhc04_case(label::String; k1::Float64, k2::Float64, z_i::Float64, n_ci::Int)
    ρ, c = 1000.0, 1000.0
    L = 0.02
    T_bot, T_top = 300.0, 400.0
    # Slowest layer sets the steady-approach horizon: τ ≈ L²/α_min.
    α_min = min(k1, k2) / (ρ * c)
    t_end = 10 * L^2 / α_min

    # Exact network with the TRUE interface position.
    d1, d2 = z_i, L - z_i
    R = d1 / k1 + d2 / k2
    q_exact = (T_top - T_bot) / R
    T_iface = T_bot + q_exact * d1 / k1
    T_exact(z) = z <= z_i ? T_bot + q_exact * z / k1 : T_iface + q_exact * (z - z_i) / k2

    function build(n_cells)
        s1 = SolidComponent(:layer1, ρ = ρ, c = c, k = k1)
        s2 = SolidComponent(:layer2, ρ = ρ, c = c, k = k2)
        mat = Material(name = :VHC04, components = (s1, s2), reactions = ())
        mesh = create_uniform_mesh(L, n_cells)
        bc_set = BoundaryConditionSet(
            thermal = (top = DirichletThermalBC(T_top), bottom = DirichletThermalBC(T_bot)),
        )
        return PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T_bot,
            ξ_initial = [
                z -> step_profile(z; z_interface = z_i, value_left = ρ, value_right = 0.0),
                z -> step_profile(z; z_interface = z_i, value_left = 0.0, value_right = ρ),
            ],
            tspan = (0.0, t_end),
            output_times = [0.0, t_end / 2, t_end],
        )
    end

    function qois(sol)
        out = QoIResult[]
        z = get_cell_positions(sol)
        n = length(z)
        Δz = L / n
        Tn = @view sol.T[:, end]

        # Composition sampling quantizes the interface to the nearest face;
        # the numerics must be nodally exact against THAT network.
        z_q = round(z_i / Δz) * Δz
        Rq = z_q / k1 + (L - z_q) / k2
        q_q = (T_top - T_bot) / Rq
        T_iface_q = T_bot + q_q * z_q / k1
        T_quant(zz) = zz <= z_q ? T_bot + q_q * zz / k1 : T_iface_q + q_q * (zz - z_q) / k2

        err = maximum(abs(Tn[i] - T_quant(z[i])) for i in eachindex(z))
        push!(
            out,
            QoIResult(
                name = "L∞ steady profile vs (quantized) network",
                value = err,
                error = err,
                tol = 1e-8,  # observed 1.5e-11 (a), 3.0e-11 (b), same for c
                comment = "harmonic-mean face k is NODALLY EXACT for face-aligned jumps; " *
                          "k-ratio $(k1 / k2):1",
            ),
        )

        # Boundary fluxes from the discrete half-cell Dirichlet conductance;
        # both must equal the quantized-network flux at steady state.
        q_bot = k1 * (Tn[1] - T_bot) / (Δz / 2)
        q_top = k2 * (T_top - Tn[n]) / (Δz / 2)
        for (name, q_num) in (("bottom", q_bot), ("top", q_top))
            push!(
                out,
                QoIResult(
                    name = "steady $(name)-face flux vs network",
                    value = q_num,
                    exact = q_q,
                    error = abs(q_num - q_q) / q_q,
                    tol = 1e-8,  # observed ≤ 4.1e-10
                    comment = "both faces must carry the exact series flux",
                ),
            )
        end

        # For an off-face interface, the user-facing smearing error is the
        # true-vs-quantized network gap, bounded a priori by the flux
        # sensitivity |dq/dz_i|·Δz/2: |1/k₁ − 1/k₂|·(Δz/2)/R_true.
        if z_q != z_i
            smear = abs(q_q - q_exact) / q_exact
            bound = abs(1 / k1 - 1 / k2) * (Δz / 2) / R
            push!(
                out,
                QoIResult(
                    name = "interface-quantization smearing (documented)",
                    value = smear,
                    exact = 0.0,
                    error = smear,
                    tol = bound,  # a-priori Δz/2 sensitivity bound (2.1% here); observed 1.64%
                    comment = "THE laminate-off-mesh error users see: flux off by " *
                              "$(round(smear * 100, digits = 2))% at Δz = $(Δz) m",
                ),
            )
        end

        steadiness = maximum(abs.(sol.T[:, end] .- sol.T[:, end-1]))
        push!(
            out,
            QoIResult(
                name = "steadiness (last half-interval drift)",
                value = steadiness,
                error = steadiness,
                tol = 1e-6,  # observed ≤ 3.1e-9 K
                comment = "t_end = 10·L²/α_min",
            ),
        )
        return out
    end

    return VerificationCase(
        id = "V-HC-04$(label)",
        title = "two-layer slab, k-ratio $(k1 / k2):1" *
                (label == "c" ? " (off-face interface)" : ""),
        tags = [:heat_conduction, :laminate, :harmonic_mean, :steady],
        references = [
            "series thermal-resistance network (exact)",
            "FDS Verification Guide `insulated_steel_plate` / `ht_pile`",
        ],
        build = build,
        solve_kwargs = (abstol = 1e-10, reltol = 1e-8),
        qois = qois,
        ci = CIConfig(n_cells = n_ci, budget_s = 20.0),
        report = (
            profile_times = [t_end],
            qoi_curves = (:T_profile_steady,),
            curves = sol -> begin
                z = get_cell_positions(sol)
                [
                    ReportCurve(
                        name = :T_profile_steady,
                        title = "V-HC-04$(label): steady two-layer profile " *
                                "(k-ratio $(k1 / k2):1)",
                        xlabel = "z (m)",
                        ylabel = "T (K)",
                        series = [
                            CurveSeries(
                                label = "steady solver profile",
                                x = z,
                                y = sol.T[:, end],
                                # true-interface network; for (c) the visible
                                # offset IS the documented quantization error
                                exact = T_exact,
                            ),
                        ],
                    ),
                ]
            end,
        ),
    )
end

# a/b: interface exactly on a face (z_i = 0.01, n even ⇒ face at L/2).
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _vhc04_case("a"; k1 = 1.0, k2 = 0.1, z_i = 0.01, n_ci = 40),
)
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _vhc04_case("b"; k1 = 1.0, k2 = 0.01, z_i = 0.01, n_ci = 40),
)
# c: interface prescribed 0.2 mm off the nearest face (Δz = 0.5 mm, n = 40):
#    quantizes to the z = 0.010 face; measured against the TRUE-z_i network.
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _vhc04_case("c"; k1 = 1.0, k2 = 0.1, z_i = 0.0102, n_ci = 40),
)
