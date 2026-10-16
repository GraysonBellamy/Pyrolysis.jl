# =============================================================================
# V-CPL-03 — Shrink/swell exact geometric endpoints (FDS `shrink_swell` analog)
# =============================================================================
#
# FDS verifies exact final thickness/density for ρ_product ≷ ρ_reactant
# permutations. Under OUR volumetric-strain model the endpoint follows from
# (technical reference Ch. 10/11; src/Physics/volume_change.jl):
#
#   V̇/V = θ = Σⱼ (γⱼ/ρⱼ) ω̇ⱼ
#
# with ω̇ⱼ the reaction mass rates [kg/m³/s]. For a single reaction
# A → Σⱼ yⱼ Bⱼ (Σyⱼ = 1) at rate r: ω̇_A = −r, ω̇_Bⱼ = +yⱼ r, so per cell
#
#   dV = (γ_A/ρ_A − Σⱼ yⱼ γⱼ/ρⱼ) dm_A        (r V dt = −dm_A)
#
# — LINEAR in the converted mass, hence path-independent: the endpoint is
# exact regardless of rate history, temperature, or spatial distribution.
# Integrating from m_A = ξ_A0·V₀ to 0, with the standard initialization
# ξ_A0 = ρ_A (bulk density) and γ_A = 1:
#
#   V_f/V₀ = 1 − γ_A + Σⱼ yⱼ γⱼ ρ_A/ρⱼ  =  Σⱼ yⱼ γⱼ ρ_A/ρⱼ
#
# and a pure solid product ends exactly at bulk density ρ_B/γ_B. All ρ here
# are the components' BULK densities (the strain kernel uses bulk ρ, not
# ρ_intrinsic — see the SolidComponent docstring note).
#
# NOTE a γ > 1 intumescent branch is NOT representable: γ is
# validated to [0, 1] (universal volume-fraction semantics shared with
# mixture density and conductivity mixing). Swelling is exercised through
# ρ_B < ρ_A (P1) and the partial-contribution branch through γ_B ∈ (0, 1)
# (P6). Documented as a model limitation.
#
# Permutations (ρ_A = 1000, γ_A = 1, isothermal via E = 0, h = 0, everything
# impermeable). Mass QoIs: EVERY permutation asserts total-mass conservation
# to integrator tolerance — every component rides the mesh and carries the
# dilation −ξθ (tech-ref §11.2.3), so the gas held in a collapsing
# impermeable matrix compresses cell by cell (ξ_g ↑ as V ↓) and leaves only
# through the (closed) faces. P4 is the sealed-column anchor for that
# bookkeeping. Gas-producing permutations additionally assert the SOLID
# ledger (char endpoint = y_char·m_A0 exactly).
#
#   P1 swell:            A → B(ρ=500)               V_f/V₀ = 2
#   P2 shrink:           A → B(ρ=2000)              V_f/V₀ = 1/2
#   P3 neutral:          A → B(ρ=1000)              V_f/V₀ = 1
#   P4 char + gas:       A → 0.3 C(ρ=1000) + 0.7 G  V_f/V₀ = 0.3
#   P5 mass loss, ΔV=0:  A → 0.5 C(ρ=500)  + 0.5 G  V_f/V₀ = 1
#   P6 partial γ:        A → B(ρ=1000, γ=0.5)       V_f/V₀ = 1/2
#   P7 WithChi:          P5 chemistry + lateral law A/A₀ = 1 − 0.5·χ̄
#        χ̄_f = 0.5 (released dry-mass fraction) ⇒ A_f/A₀ = 0.75 exactly,
#        L_f/L₀ = (V_f/V₀)/(A_f/A₀) = 4/3.
#
# P1 (the swelling branch, exercising mesh GROWTH) runs through the standard
# harness path; the remaining permutations are solved inside `qois` (each is
# a ~1 s, 8-cell isothermal solve) so the whole family stays one case.
#
# First-order kinetics with A_r·t_end = 16 ⇒ residual conversion e⁻¹⁶ ≈ 1e-7
# bounds the endpoint error; the depletion limiter is disabled so the decay
# stays exactly first-order (with it, the tail below ~1 kg/m³ slows and the
# residual at fixed t_end grows ~1000×).

function _vcpl03_problem(perm; n_cells = 8, n_outputs = 3)
    ρ_A = 1000.0
    T_iso = 400.0
    A_r = 2.0
    t_end = 8.0
    L = 2e-3

    solidA = SolidComponent(:A, ρ = ρ_A, c = 1000.0, k = 0.5)
    gasG(name = :G) = GaseousComponent(name, M = 0.05, c = 1000.0, k = 0.02, λ = 1e-8)

    comps, rxn, exact_Vf, chi_law = if perm == :P1
        (
            (solidA, SolidComponent(:B, ρ = 500.0, c = 1000.0, k = 0.5)),
            Reaction(:r, 1 => 2; A = A_r, E = 0.0, h = 0.0),
            2.0,
            nothing,
        )
    elseif perm == :P2
        (
            (solidA, SolidComponent(:B, ρ = 2000.0, c = 1000.0, k = 0.5)),
            Reaction(:r, 1 => 2; A = A_r, E = 0.0, h = 0.0),
            0.5,
            nothing,
        )
    elseif perm == :P3
        (
            (solidA, SolidComponent(:B, ρ = 1000.0, c = 1000.0, k = 0.5)),
            Reaction(:r, 1 => 2; A = A_r, E = 0.0, h = 0.0),
            1.0,
            nothing,
        )
    elseif perm == :P4
        (
            (solidA, SolidComponent(:C, ρ = 1000.0, c = 1000.0, k = 0.5), gasG()),
            Reaction(:r, 1 => (2, 3); A = A_r, E = 0.0, h = 0.0, yields = (0.3, 0.7)),
            0.3,
            nothing,
        )
    elseif perm == :P5
        (
            (solidA, SolidComponent(:C, ρ = 500.0, c = 1000.0, k = 0.5), gasG()),
            Reaction(:r, 1 => (2, 3); A = A_r, E = 0.0, h = 0.0, yields = (0.5, 0.5)),
            1.0,
            nothing,
        )
    elseif perm == :P6
        (
            (solidA, SolidComponent(:B, ρ = 1000.0, c = 1000.0, k = 0.5, γ = 0.5)),
            Reaction(:r, 1 => 2; A = A_r, E = 0.0, h = 0.0),
            0.5,
            nothing,
        )
    elseif perm == :P7
        (
            (solidA, SolidComponent(:C, ρ = 500.0, c = 1000.0, k = 0.5), gasG()),
            Reaction(:r, 1 => (2, 3); A = A_r, E = 0.0, h = 0.0, yields = (0.5, 0.5)),
            1.0,
            χ̄ -> 1.0 - 0.5 * χ̄,
        )
    else
        error("unknown permutation $perm")
    end

    nc = length(comps)
    mat = Material(
        name = Symbol(:VCPL03_, perm),
        components = comps,
        reactions = (rxn,),
        depletion_limiter = DepletionLimiter(enabled = false),
        lateral_shrinkage_law = chi_law,
    )
    ξinit_t = ntuple(i -> i == 1 ? ρ_A : 0.0, nc)
    mesh = create_uniform_mesh(L, n_cells, Val(nc); T_initial = T_iso, ξ_initial = ξinit_t)
    bc_set = BoundaryConditionSet(
        thermal = (top = DirichletThermalBC(T_iso), bottom = DirichletThermalBC(T_iso)),
    )
    problem = PyrolysisProblem(
        mesh = mesh,
        material = mat,
        bc_set = bc_set,
        T_initial = z -> T_iso,
        ξ_initial = [i == 1 ? (z -> ρ_A) : (z -> 0.0) for i = 1:nc],
        tspan = (0.0, t_end),
        output_times = collect(range(0.0, t_end, length = n_outputs)),
    )
    return (; problem, exact_Vf, L, ρ_A, perm, A_r, t_end)
end

function _vcpl03_solve(perm; n_cells = 8, n_outputs = 3)
    r = _vcpl03_problem(perm; n_cells = n_cells, n_outputs = n_outputs)
    sol = solve(
        r.problem;
        use_ale = true,
        min_thickness = 1e-5,
        progress = false,
        abstol = 1e-10,
        reltol = 1e-8,
    )
    return (; sol, r.exact_Vf, r.L, r.ρ_A, r.perm, r.A_r)
end

function _vcpl03_qois_for(r, out::Vector{QoIResult}; tol_thick, tol_mass)
    sol = r.sol
    th_ratio = sol.extras.thickness_history[end] / r.L
    A_ratio =
        sol.extras.cross_section_area_history[end] /
        sol.extras.cross_section_area_history[1]
    exact_L = r.perm == :P7 ? r.exact_Vf / 0.75 : r.exact_Vf
    err_th = abs(th_ratio - exact_L) / exact_L
    push!(
        out,
        QoIResult(
            name = "$(r.perm): final thickness ratio",
            value = th_ratio,
            exact = exact_L,
            error = err_th,
            tol = tol_thick,
            comment = "path-independent endpoint of θ = Σγⱼω̇ⱼ/ρⱼ (case header)",
        ),
    )
    # Total-mass conservation: exact for every permutation — impermeable
    # everywhere; every component rides the mesh (case header).
    mass_err =
        abs(sol.extras.total_mass[end] - sol.extras.total_mass[1]) /
        sol.extras.total_mass[1]
    push!(
        out,
        QoIResult(
            name = "$(r.perm): total mass conservation",
            value = sol.extras.total_mass[end],
            exact = sol.extras.total_mass[1],
            error = mass_err,
            tol = tol_mass,
            comment = r.perm in (:P4, :P5, :P7) ?
                      "held gas compresses with its cell under collapse; " *
                      "the sealed P4 corner anchors the dilation bookkeeping" :
                      "solid-only permutation; impermeable everywhere",
        ),
    )
    if r.perm in (:P4, :P5, :P7)
        # Solid ledger: final char mass per unit initial area = y_char·m_A0.
        y_char = 0.5 * (r.perm == :P4 ? 0.6 : 1.0)  # P4: 0.3, P5/P7: 0.5
        nodes = filter(!isnan, @view sol.z[:, size(sol.z, 2)])
        A_end = sol.extras.cross_section_area_history[end]
        m_char =
            sum(sol.ξ[2, i, end] * (nodes[i+1] - nodes[i]) for i = 1:(length(nodes)-1)) *
            A_end / sol.extras.cross_section_area_history[1]
        m_exact = y_char * r.ρ_A * r.L
        err = abs(m_char - m_exact) / m_exact
        push!(
            out,
            QoIResult(
                name = "$(r.perm): final char mass (solid ledger)",
                value = m_char,
                exact = m_exact,
                error = err,
                tol = tol_mass,
                comment = "char endpoint is path-independent; complements " *
                          "the total-mass anchor",
            ),
        )
    end
    if r.perm == :P7
        err_A = abs(A_ratio - 0.75) / 0.75
        push!(
            out,
            QoIResult(
                name = "P7: final lateral area ratio A_f/A₀ = law(χ̄ = 1/2)",
                value = A_ratio,
                exact = 0.75,
                error = err_A,
                tol = tol_thick,
                comment = "χ̄_f = released dry-mass fraction = y_gas = 0.5",
            ),
        )
    end
    return out
end

function _verification_case_vcpl03()
    # P1 (swelling branch) goes through the standard harness solve path.
    build(n_cells) = _vcpl03_problem(:P1; n_cells = n_cells).problem

    function qois(sol)
        out = QoIResult[]
        # tolerances: residual conversion e⁻¹⁶ ≈ 1.1e-7 plus integrator error;
        # pin at 3× observed via the print run.
        tol_thick = 3e-6   # 3× observed max: 7.1e-7 (P1)
        tol_mass = 2e-6    # 3× observed max: 6.5e-7 (P1 total mass)
        # P1 from the harness solve
        r1 = (; sol, exact_Vf = 2.0, L = 2e-3, ρ_A = 1000.0, perm = :P1)
        _vcpl03_qois_for(r1, out; tol_thick = tol_thick, tol_mass = tol_mass)
        for perm in (:P2, :P3, :P4, :P5, :P6, :P7)
            r = _vcpl03_solve(perm)
            if r.sol.retcode != :Success
                push!(
                    out,
                    QoIResult(
                        name = "$(perm): retcode",
                        value = 0.0,
                        error = 1.0,
                        tol = 0.5,
                        comment = "solve returned $(r.sol.retcode)",
                    ),
                )
                continue
            end
            _vcpl03_qois_for(r, out; tol_thick = tol_thick, tol_mass = tol_mass)
        end
        return out
    end

    return VerificationCase(
        id = "V-CPL-03",
        title = "shrink/swell exact endpoints (6 permutations + WithChi law)",
        tags = [:ale, :volume_change, :swelling, :chi, :mass_conservation],
        references = [
            "FDS Verification Guide `shrink_swell` (endpoint concept)",
            "Technical reference Ch. 10/11 (volumetric strain model)",
        ],
        build = build,
        solve_kwargs = (
            use_ale = true,
            min_thickness = 1e-5,
            abstol = 1e-10,
            reltol = 1e-8,
        ),
        qois = qois,
        ci = CIConfig(n_cells = 8, budget_s = 60.0),  # 7 solves; dominated by
        # per-permutation compile (each (NC, chi) trait combo specializes)
        report = (
            qoi_curves = (:thickness_histories,),
            # The endpoint derivation is linear in the converted mass (case
            # header), so the WHOLE thickness path is exact, not just the
            # endpoint: L(t)/L₀ = 1 + (V_f/V₀ − 1)·X(t) with the first-order
            # isothermal conversion X(t) = 1 − e^(−A_r·t) (P7 additionally
            # divides by the lateral law A/A₀ = 1 − 0.25·X). Every permutation
            # is re-solved with dense outputs; the harness-provided P1 solve
            # only carries the 3 assert-mode snapshots.
            curves = _ -> begin
                series = CurveSeries[]
                for perm in (:P1, :P2, :P3, :P4, :P5, :P6, :P7)
                    r = _vcpl03_solve(perm; n_outputs = 33)
                    r.sol.retcode == :Success || continue
                    ex = let A_r = r.A_r, Vf = r.exact_Vf, isP7 = (r.perm === :P7)
                        t -> begin
                            X = 1.0 - exp(-A_r * t)
                            V = 1.0 + (Vf - 1.0) * X
                            isP7 ? V / (1.0 - 0.25 * X) : V
                        end
                    end
                    push!(
                        series,
                        CurveSeries(
                            label = string(r.perm),
                            x = r.sol.t,
                            y = r.sol.extras.thickness_history ./ r.L,
                            exact = ex,
                        ),
                    )
                end
                # FDS 6.11.0 overlay (reference/fds/custom/): the shipped
                # shrink_swell case heats with an Arrhenius transient (only
                # endpoints comparable), so vcpl03_shrink_swell_paths.fds
                # reruns FDS isothermally with our P1/P2 densities and
                # A_r = 2 s⁻¹, E = 0 — the whole path then follows the same
                # exact curve. devc values average over the 0.1 s output
                # interval (plot at t − 0.05 s).
                fds = reference_columns(
                    "fds/custom/vcpl03_shrink_swell_paths_devc.csv",
                    "Time",
                    "L_P2",
                    "L_P1",
                )
                t_fds = fds[1][2:end] .- 0.05
                for (perm, col, Vf) in ((:P1, 3, 2.0), (:P2, 2, 0.5))
                    push!(
                        series,
                        CurveSeries(
                            label = "FDS 6.11.0 $(perm)",
                            crosscode = true,
                            x = t_fds,
                            y = fds[col][2:end] ./ 0.002,
                            exact = t -> 1.0 + (Vf - 1.0) * (1.0 - exp(-2.0 * t)),
                        ),
                    )
                end
                [
                    ReportCurve(
                        name = :thickness_histories,
                        title = "V-CPL-03: thickness-ratio histories, all permutations " *
                                "(exact linear-in-conversion paths)",
                        xlabel = "t (s)",
                        ylabel = "L(t)/L₀",
                        series = series,
                    ),
                ]
            end,
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vcpl03)
