# =============================================================================
# V-HC-03 — Finite slab, constant flux, adiabatic back (exact series)
# =============================================================================
#
# Complements the semi-infinite cases at LONG times: the thermal wave reaches
# the back face — the regime cone-calorimeter simulations actually live in.
# Constant q″ at the top face, insulated back; exact solution
# T = T₀ + (q″L/k)·φ(x̃, Fo) with φ from exact/slab_series.jl (C&J §3.8),
# x̃ = z/L from the adiabatic back face.
#
# Two parameter-free exact targets beyond the profiles:
# - front–back plateau: ΔT → q″L/2k for Fo ≫ 1 (100 K here),
# - global energy: mean temperature rises exactly at q″/(ρcL) (0.25 K/s here).

function _verification_case_vhc03()
    ρ, c, k = 1000.0, 1000.0, 0.5
    α = k / (ρ * c)
    L = 0.02
    T0 = 300.0
    q = 5000.0
    Fos = [0.1, 0.5, 1.5]
    times = [Fo * L^2 / α for Fo in Fos]  # 80, 400, 1200 s

    T_exact(x̃, Fo) = T0 + (q * L / k) * slab_constflux_phi(x̃, Fo)

    function build(n_cells)
        solid = SolidComponent(:solid, ρ = ρ, c = c, k = k)
        mat = Material(name = :VHC03, components = (solid,), reactions = ())
        mesh = create_uniform_mesh(L, n_cells, Val(1); T_initial = T0, ξ_initial = (ρ,))
        bc_set =
            BoundaryConditionSet(thermal = (top = HeatFluxBC(q), bottom = AdiabaticBC()))
        return PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T0,
            ξ_initial = [z -> ρ],
            tspan = (0.0, times[end]),
            output_times = vcat(0.0, times),
        )
    end

    function qois(sol)
        out = QoIResult[]
        z = get_cell_positions(sol)
        n = length(z)
        for (j, Fo) in enumerate(Fos)
            Tn = @view sol.T[:, j+1]
            err = maximum(abs(Tn[i] - T_exact(z[i] / L, Fo)) for i in eachindex(z))
            push!(
                out,
                QoIResult(
                    name = "L∞ T-profile error, Fo=$(Fo)",
                    value = err,
                    error = err,
                    tol = 0.006,  # 3× observed max: 2.0e-3 K at Fo=0.1 (n=80, 2026-07-06)
                    comment = "observed 2.0e-3/1.5e-3/1.3e-3 K at Fo=0.1/0.5/1.5",
                ),
            )
        end

        # Front–back ΔT at Fo = 1.5, compared cell-center to cell-center
        # (series term at Fo = 1.5 is < 4e-7, i.e. the plateau q″L/2k up to
        # the finite offset of the first/last cell centers, which the exact
        # curve accounts for).
        ΔT_num = sol.T[n, end] - sol.T[1, end]
        ΔT_exact = T_exact(z[n] / L, Fos[end]) - T_exact(z[1] / L, Fos[end])
        push!(
            out,
            QoIResult(
                name = "front–back ΔT plateau, Fo=1.5",
                value = ΔT_num,
                exact = ΔT_exact,
                error = abs(ΔT_num - ΔT_exact),
                tol = 1e-5,  # observed 3.2e-7 K (n=80, 2026-07-06); generous vs integrator tol
                comment = "steady-shape plateau q″L/2k = 100 K between faces",
            ),
        )

        # Global energy: mean rise between Fo = 0.5 and Fo = 1.5 is exactly
        # q″Δt/(ρcL) = 200 K (uniform mesh → plain mean is the volume mean).
        rise_num = mean(sol.T[:, 4]) - mean(sol.T[:, 3])
        rise_exact = q * (times[3] - times[2]) / (ρ * c * L)
        push!(
            out,
            QoIResult(
                name = "mean-temperature rise (global energy)",
                value = rise_num,
                exact = rise_exact,
                error = abs(rise_num - rise_exact),
                tol = 1e-5,  # observed 5.7e-14 K of 200 K (n=80, 2026-07-06); integrator headroom
                comment = "exact constant q″/(ρcL) = 0.25 K/s; integrator-tolerance limited",
            ),
        )
        return out
    end

    # L2 profile error at Fo = 0.1 (output index 2), the most transient state.
    function l2_error_fo01(sol)
        z = get_cell_positions(sol)
        Tn = @view sol.T[:, 2]
        Δz = L / length(z)
        return sqrt(sum((Tn[i] - T_exact(z[i] / L, 0.1))^2 * Δz for i in eachindex(z)) / L)
    end

    return VerificationCase(
        id = "V-HC-03",
        title = "finite slab, constant flux, adiabatic back",
        tags = [:heat_conduction, :heat_flux_bc, :energy_balance],
        references = ["Carslaw & Jaeger (1959) §3.8 (constant-flux slab series)"],
        build = build,
        solve_kwargs = (abstol = 1e-10, reltol = 1e-8),
        qois = qois,
        ci = CIConfig(
            n_cells = 80,
            budget_s = 20.0,
            convergence = ConvergenceSpec(
                # observed rates 1.9999, 1.9999 (2026-07-06)
                n_cells = [40, 80, 160],
                error_fn = l2_error_fo01,
                expected_order = 2.0,
                min_order = 1.5,
                min_mean_order = 1.8,
                solve_kwargs = (abstol = 1e-12, reltol = 1e-10),
            ),
        ),
        report = (
            n_cells_ladder = [20, 40, 80, 160, 320],
            profile_times = times,
            qoi_curves = (:T_profiles, :front_back_dT_history, :mean_T_history),
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vhc03)
