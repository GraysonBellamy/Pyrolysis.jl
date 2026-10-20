# =============================================================================
# V-CPL-06 — gas sensible-enthalpy (S_conv) ledger exactness
# =============================================================================
#
# A char-layer column: reactive solid band at z/L ∈ [0.2, 0.4] gasifies with
# h = 0 (isolates the sensible ledger) at ~500–600 K; the gas diffuses up
# through the ~880 K hot layer and vents at the top. Gas pickup between
# formation and exit is ≈ 330 K — the cone char-layer configuration, chosen
# so any gas-enthalpy accounting defect is maximally visible.
#
# V-CPL-06a (Eulerian): at quasi-steady state the discrete S_conv scheme must
# satisfy summation-by-parts EXACTLY (holdup-corrected):
#
#   ∫S_conv dV·dt = −c_g·( ∫J_top·T_ftop dt − ∫Σᵢ(genᵢ − ξ̇_g,i)·Tᵢ·Vᵢ dt )
#
# — the energy S_conv removes from the field is identically the energy the
# discrete species ledger says the gas picks up between formation and exit —
# and match the physical reheat reference −c_g·Σ gen·(T_exit − T_form)·V to
# quasi-steady accuracy. Together these pin the S_conv formulation AND its
# discretization (tech-ref §7).
#
# V-CPL-06b (ALE): same column with the fuel band collapsing as it gasifies
# (γ_gas = 0 ⇒ the overburden subsides — cone-charring kinematics). Pins:
# - the S_conv-vs-reheat defect stays at the (converged) quasi-steady level;
# - the ALE energy budget closes: ΔE = ∫(q_bnd + Q_rxn + S_conv + hand-off +
#   gas hand-off + dilation) dt, with E = Σ ρcp·T·V over every component
#   (the gas hand-off c_g·Σ ξ̇_g·T·V is the sensible heat the stored pore gas
#   gains or loses with its concentration) — and with NO mesh-motion
#   T-advection term: the
#   relative-velocity formulation (tech-ref §11.2.2) has none in production,
#   and this case asserts the residual T-advection accumulator is exactly
#   zero. A full-mesh-velocity `+w·∂T/∂z` term would show up here as an
#   energy sink of ≈ 66 kJ per kg of volatiles (subsiding hot overburden,
#   ΔT ≈ 330 K across the layer).
#
# Ledger accumulation runs in a per-accepted-step DiscreteCallback around the
# pure residual! (wrapped-RHS pattern; no solver changes), with TRAPEZOIDAL
# time quadrature aligned with the E₁ snapshot — a first-order rectangle rule
# leaves an O(dt·Ẋ) budget residual (~30% of ∫q_bnd here) over this transient
# window, which is pure quadrature, not physics. h = 0 ⇒ Q_rxn ≡ 0, kept in
# the ledger as a guard.

using OrdinaryDiffEq
using Pyrolysis.Solver: KenCarp4, AutoFiniteDiff, create_step_update_callback
using Pyrolysis.Problem: apply_initial_conditions!

const _VCPL06_CGAS = 1000.0

function _vcpl06_problem(n_cells; t_end = 240.0)
    L = 0.01
    c_gas = _VCPL06_CGAS
    # λ on ALL components: PARALLEL mixing is volume-fraction weighted, so the
    # solid-dominated column needs solid λ > 0 for interior gas transport
    # (same trick as V-GT-01).
    matrixc = SolidComponent(:matrix, ρ = 2000.0, c = 1000.0, k = 1.0, λ = 1e-4)
    fuel = SolidComponent(:fuel, ρ = 2000.0, c = 1000.0, k = 1.0, λ = 1e-4)
    gas = GaseousComponent(:gas, M = 0.05, c = c_gas, k = 0.02, λ = 1e-4)
    rxn = Reaction(:gasify, 2 => 3; A = 9.3e3, E = 8.0e4, h = 0.0)
    mat = Material(name = :VCPL06, components = (matrixc, fuel, gas), reactions = (rxn,))
    mesh = create_uniform_mesh(
        L,
        n_cells,
        Val(3);
        T_initial = 300.0,
        ξ_initial = (400.0, 0.0, 0.0),
    )
    bc_set = BoundaryConditionSet(
        thermal = (
            top = ConvectiveBC(h = 1000.0, T_inf = 900.0),
            bottom = ConvectiveBC(h = 1000.0, T_inf = 300.0),
        ),
        mass = (top = ConvectiveMassBC(h_m = 1.0, ξ_inf = 0.0), bottom = ImpermeableBC()),
    )
    band(z) = (0.2 * L <= z <= 0.4 * L) ? 100.0 : 0.0
    return PyrolysisProblem(
        mesh = mesh,
        material = mat,
        bc_set = bc_set,
        T_initial = z -> 300.0,
        ξ_initial = [z -> 400.0, band, z -> 0.0],
        tspan = (0.0, t_end),
        output_times = [t_end],
    )
end

# Wrapped-RHS solve with a per-accepted-step energy/mass ledger. Returns a
# NamedTuple consumed by the qois closures (retcode + window integrals).
# `gas` is the gas component's index; the gas and every condensed component
# carry the constant heat capacities `c_gas` and `c_solid`.
function _vcpl06_ledger_solve(
    problem;
    use_ale::Bool,
    t_win = 120.0,
    gas::Int = 3,
    c_gas = _VCPL06_CGAS,
    c_solid = 1000.0,
)
    mat = problem.material
    apply_initial_conditions!(problem)
    def = to_problem_def(problem; use_ale = use_ale)
    ws = build_workspace(problem, def)
    u0 = create_state_vector(ws.layout, ws.mesh)
    n = ws.layout.n_cells
    NC = Pyrolysis.Materials.n_components(mat)
    gidx = gas

    rhs! = (du, u, p, t) -> residual!(du, u, def, t, ws)

    # Integrand slots (fixed order): S_conv, Q_rxn, q_bnd, JtopTt, genT,
    # xidotT, m_gen, m_out, gen_reheat, handoffE, ale_T, dil.
    NL = 12
    acc = zeros(NL)
    cur = zeros(NL)
    prev = zeros(NL)
    st = Dict{Symbol,Float64}(:last_t => 0.0, :E1 => NaN, :E2 => NaN, :in_win => 0.0)
    du_scratch = similar(u0)

    function ledger!(integ)
        t = integ.t
        dt = t - st[:last_t]
        st[:last_t] = t
        dt <= 0 && return nothing
        residual!(du_scratch, integ.u, def, t, ws)  # populate caches at (u, t)
        rc = ws.rhs_cache
        sc = ws.state_cache
        f_top = ws.top_face_ids[1]
        f_bot = ws.bottom_face_ids[1]
        i_top = ws.mesh.faces[f_top].cell_ids[1]
        E = 0.0
        @inbounds for i = 1:n
            E += ws.prop_cache.ρcp[i] * sc.T[i] * sc.volume[i]
        end
        if t >= t_win
            A_top = sc.A_face[f_top]
            # Face temperature at the solved surface slope — the same
            # construction as the S_conv boundary blend.
            s_face = -rc.q_total[f_top] / ws.prop_cache.k[i_top]
            dz_half = sc.volume[i_top] / A_top / 2
            T_ftop = sc.T[i_top] + s_face * dz_half
            J_top = rc.J_gas[gidx, f_top]
            Sc = 0.0
            Qr = 0.0
            genT = 0.0
            xdT = 0.0
            mg = 0.0
            reheat = 0.0
            hoE = 0.0
            aleT = 0.0
            dil = 0.0
            @inbounds for i = 1:n
                V = sc.volume[i]
                T = sc.T[i]
                Sc += rc.S_conv[i] * V
                Qr += rc.Q_rxn[i] * V
                g = rc.dξ_rxn[gidx, i]
                genT += g * T * V
                xdT += ws.dξ_work[gidx, i] * T * V
                mg += g * V
                reheat += g * (T_ftop - T) * V
                # Condensed sensible hand-off: ρcp changes with composition.
                for j = 1:NC
                    j == gidx && continue
                    hoE += c_solid * T * rc.dξ_rxn[j, i] * V
                end
                if use_ale
                    aleT += ws.prop_cache.ρcp[i] * ws.ale_advection[i] * V
                    for j = 1:NC
                        j == gidx && continue
                        dil += -c_solid * T * sc.ξ[j, i] * ws.θ_work[i] * V
                    end
                    dil += ws.prop_cache.ρcp[i] * T * ws.θ_work[i] * V
                end
            end
            cur[1] = Sc
            cur[2] = Qr
            cur[3] = -rc.q_total[f_top] * A_top + rc.q_total[f_bot] * sc.A_face[f_bot]
            cur[4] = J_top * T_ftop * A_top
            cur[5] = genT
            cur[6] = xdT
            cur[7] = mg
            cur[8] = J_top * A_top
            cur[9] = c_gas * reheat
            cur[10] = hoE
            cur[11] = aleT
            cur[12] = dil
            if st[:in_win] == 1.0
                # Trapezoid over [t_prev, t] — both endpoints inside the
                # window, so ∫X dt aligns exactly with E₂ − E₁.
                @inbounds for k = 1:NL
                    acc[k] += 0.5 * (prev[k] + cur[k]) * dt
                end
            else
                # First in-window step: take the E₁ snapshot here and start
                # integrating from this endpoint (the interval ending at this
                # step began before the snapshot — skip it).
                st[:E1] = E
                st[:in_win] = 1.0
            end
            st[:E2] = E
            copyto!(prev, cur)
        end
        return nothing
    end

    ledger_cb =
        DiscreteCallback((u, t, integ) -> true, ledger!; save_positions = (false, false))
    cbs = use_ale ? CallbackSet(create_step_update_callback(ws, mat), ledger_cb) : ledger_cb

    prob = ODEProblem(rhs!, u0, problem.tspan)
    ode_sol = OrdinaryDiffEq.solve(
        prob,
        KenCarp4(
            autodiff = AutoFiniteDiff(),
            linsolve = Pyrolysis.Solver.LinearSolve.LUFactorization(),
        );
        abstol = 1e-8,
        reltol = 1e-6,
        save_everystep = false,
        save_start = false,
        callback = cbs,
        maxiters = 10^7,
    )

    retcode =
        OrdinaryDiffEq.SciMLBase.successful_retcode(ode_sol) ? :Success :
        Symbol(ode_sol.retcode)
    ΔE = st[:E2] - st[:E1]
    acc_nt = (;
        S_conv = acc[1],
        Q_rxn = acc[2],
        q_bnd = acc[3],
        JtopTt = acc[4],
        genT = acc[5],
        xidotT = acc[6],
        m_gen = acc[7],
        m_out = acc[8],
        gen_reheat = acc[9],
        handoffE = acc[10],
        ale_T = acc[11],
        dil = acc[12],
    )
    return (; retcode, acc = acc_nt, ΔE, c_gas)
end

function _verification_case_vcpl06a()
    function qois(sol)
        a = sol.acc
        # SBP identity (holdup-corrected): exact property of the discrete
        # S_conv scheme — machine-precision closure expected at any n.
        S_exact = -sol.c_gas * (a.JtopTt - a.genT + a.xidotT)
        err_sbp = abs(a.S_conv - S_exact) / abs(S_exact)
        # Physical reheat reference: quasi-steady, converged at ~1e-4.
        reheat_ref = -a.gen_reheat
        err_reheat = abs(a.S_conv - reheat_ref) / abs(reheat_ref)
        # Quasi-steady mass throughput (holdup drift bounds the deviation).
        err_mass = abs(a.m_out / a.m_gen - 1.0)
        return [
            QoIResult(
                name = "S_conv summation-by-parts identity (discrete, exact)",
                value = a.S_conv,
                exact = S_exact,
                error = err_sbp,
                tol = 1e-6,  # observed ~1e-9 at n = 40..320
                comment = "energy S_conv removes ≡ energy the discrete gas " *
                          "ledger says the gas picks up",
            ),
            QoIResult(
                name = "S_conv vs physical reheat −c_g·Σ gen·(T_exit−T_form)·V",
                value = a.S_conv,
                exact = reheat_ref,
                error = err_reheat,
                tol = 2.4e-4,  # 3× observed: 7.8e-5 (n=80; stable
                # across the n = 40..320 ladder)
                comment = "quasi-steady physical reference",
            ),
            QoIResult(
                name = "mass throughput: exit flux vs generation",
                value = a.m_out,
                exact = a.m_gen,
                error = err_mass,
                tol = 2.5e-4,  # 3× observed: 7.5e-5 (n=80): holdup drift
                comment = "quasi-steady window guard",
            ),
        ]
    end

    return VerificationCase(
        id = "V-CPL-06a",
        title = "Eulerian gas-enthalpy (S_conv) ledger exactness",
        tags = [:eulerian, :gas_transport, :energy_ledger, :coupled],
        references = [
            "Technical reference §7 (S_conv), §18.4 (conservation-based verification)",
        ],
        build = _vcpl06_problem,
        solve_fn = (problem; kwargs...) -> _vcpl06_ledger_solve(problem; use_ale = false),
        qois = qois,
        ci = CIConfig(n_cells = 80, budget_s = 60.0),
        report = (;),
    )
end

function _verification_case_vcpl06b()
    function qois(sol)
        a = sol.acc
        reheat_ref = -a.gen_reheat
        defect_K = (a.S_conv - reheat_ref) / (sol.c_gas * a.m_gen)
        # ΔE = Σ[ρcp V Ṫ + T V d(ρcp)/dt + ρcp T V̇]: the first is the
        # integrated RHS, the second splits into the condensed reaction
        # hand-off, the gas hand-off c_g·ξ̇_g (full rate: reaction, flux and
        # dilation) and the condensed dilation, the third is the volume term.
        pred =
            a.q_bnd +
            a.Q_rxn +
            a.S_conv +
            a.handoffE +
            sol.c_gas * a.xidotT +
            a.ale_T +
            a.dil
        err_ledger = abs(sol.ΔE - pred) / abs(a.q_bnd)
        return [
            QoIResult(
                name = "residual mesh-motion T-advection accumulator ≡ 0",
                value = a.ale_T,
                exact = 0.0,
                error = abs(a.ale_T),
                tol = 1e-12,
                comment = "relative-velocity formulation: the production ALE " *
                          "energy equation has NO w·∂T/∂z term (a full-w term " *
                          "would sink ≈66 kJ/kg of volatiles here)",
            ),
            QoIResult(
                name = "S_conv-vs-reheat defect per kg of volatiles [K·c_g]",
                value = defect_K,
                exact = 0.0,
                error = abs(defect_K),
                tol = 0.1,  # ~4× observed: 0.026 K (n=80) on a
                # 330 K pickup — the T_s boundary-blend secondary (tech-ref
                # §12.3.2), shrinking with n
                comment = "c_g·0.03 K-scale closure on a 330 K pickup; a " *
                          "spurious full-w T-advection would read −66 K/kg here",
            ),
            QoIResult(
                name = "ALE energy-budget closure (ΔE vs integrated terms)",
                value = sol.ΔE,
                exact = pred,
                error = err_ledger,
                tol = 1.5e-4,  # 3× observed: 4.7e-5 of ∫q_bnd (n=80,
                # trapezoidal ledger; a rectangle rule leaves a ~30%
                # pure-quadrature residual on this transient window)
                comment = "q_bnd + Q_rxn + S_conv + hand-off + gas hand-off + " *
                          "dilation; NO T-advection term",
            ),
        ]
    end

    return VerificationCase(
        id = "V-CPL-06b",
        title = "ALE char-layer gas-enthalpy ledger (subsiding overburden)",
        tags = [:ale, :gas_transport, :energy_ledger, :volume_change, :coupled],
        references = [
            "Technical reference §11.2.2 (relative-velocity ALE energy), §7 (S_conv)",
        ],
        build = _vcpl06_problem,
        solve_fn = (problem; kwargs...) -> _vcpl06_ledger_solve(problem; use_ale = true),
        qois = qois,
        ci = CIConfig(n_cells = 80, budget_s = 90.0),
        report = (;),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vcpl06a)
push!(VERIFICATION_CASE_BUILDERS, _verification_case_vcpl06b)
