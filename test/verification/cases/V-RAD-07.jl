# =============================================================================
# V-RAD-07 — in-depth heating with emission (sub-surface peak)
# =============================================================================
#
# A gray slab (δ = 1 cm, κ = 1000 m⁻¹ so κδ = 10, ρ = 1000 kg/m³,
# c = 1400 J/kg·K, k = 0.2 W/m·K, ε = 1) at T₀ = 300 K under a collimated
# heater of 25 kW/m², facing surroundings at T₀ through its exposed face and
# losing heat there by convection (h = 20 W/m²·K to T₀); the back is
# adiabatic and reflects all radiation (a mirror). The heater is absorbed in
# depth as e^{−κx}, every cell emits σT⁴ and absorbs what the others emit,
# and the exposed face both cools by convection and lets the slab's emission
# out, so the temperature peaks below the surface: 8 K above the surface
# temperature at 180 s, 0.6 mm deep, and deepening with time.
#
# Reference: the same discretization solved by an independent method of
# lines (reference/radiation/vrad07_reference.jl, which does not use the
# package): the radiation rows of exact/two_stream.jl eliminated at every
# evaluation, so the reference is an ordinary differential system in the
# cell temperatures, integrated by Kvaerno5 at reltol 1e-12 (its distance to
# FBDF at the same tolerance is at most 2.5e-7 K). On the same mesh the two
# differ only by the time integration. The order is measured against the
# reference's 1600-cell solution, interpolated to the cell centres.
#
# The exposed face is black (ε = 1), so its transmissivity is one both ways:
# the heater enters whole and the slab's emission leaves whole. Under
# radiative transfer the surface balance holds convection only; the
# radiation is in the cells' sources.

function _vrad07_reference(n)
    names, cols = read_reference_csv("radiation/vrad07_n$(n).csv")
    return (t = cols[1], T_s = cols[2], T = permutedims(reduce(hcat, cols[3:end])))
end

# Linear interpolation of the samples (z, y) at z₀.
function _vrad07_interpolate(z, y, z0)
    j = clamp(searchsortedlast(z, z0), 1, length(z) - 1)
    w = (z0 - z[j]) / (z[j+1] - z[j])
    return (1 - w) * y[j] + w * y[j+1]
end

function _verification_case_vrad07()
    ρ, c, k = 1000.0, 1400.0, 0.2
    δ, κ = 0.01, 1000.0
    q0, h, T0 = 25.0e3, 20.0, 300.0
    saves = collect(0.0:5.0:180.0)
    times = [30.0, 60.0, 90.0, 180.0]
    column(t) = findfirst(==(t), saves)

    references = Dict(n => _vrad07_reference(n) for n in (50, 100, 200))
    fine = _vrad07_reference(1600)
    z_fine = [(j - 0.5) * δ / 1600 for j = 1:1600]
    # The continuum proxy at time t (one of `times`) and height z.
    continuum(z, t) = _vrad07_interpolate(z_fine, fine.T[:, findfirst(==(t), fine.t)], z)

    function build(n_cells)
        solid = SolidComponent(:solid, ρ = ρ, c = c, k = k, ε = 1.0, α = κ / ρ)
        material = Material(name = :VRAD07, components = (solid,), reactions = ())
        bc_set = BoundaryConditionSet(
            thermal = (
                top = RadiativeFluxBC(q0; distribution = Collimated()) +
                      RadiativeBC(T_inf = T0) +
                      ConvectiveBC(h = h, T_inf = T0),
                bottom = AdiabaticBC(),
            ),
        )
        return PyrolysisProblem(
            mesh = uniform_mesh(δ, n_cells),
            material = material,
            bc_set = bc_set,
            T_initial = z -> T0,
            ξ_initial = [z -> ρ],
            tspan = (0.0, saves[end]),
            output_times = saves,
        )
    end

    function qois(sol)
        n = size(sol.T, 1)
        ref = references[n]
        @assert ref.t == sol.t
        out = QoIResult[]
        for t in times
            j = column(t)
            err = maximum(abs, sol.T[:, j] .- ref.T[:, j])
            push!(
                out,
                QoIResult(
                    name = "L∞ T against the reference, t=$(t) s",
                    value = err,
                    error = err,
                    tol = 1.2e-6,  # 3× the largest observed, 3.9e-7 K (60 s, temperature form)
                    comment = "same mesh: the two differ by the time integration only",
                ),
            )
        end
        for t in times
            j = column(t)
            push!(
                out,
                QoIResult(
                    name = "front-face temperature, t=$(t) s",
                    value = sol.surface_T[j],
                    exact = ref.T_s[j],
                    error = abs(sol.surface_T[j] - ref.T_s[j]),
                    tol = 1.2e-6,  # 3× the largest observed, 3.9e-7 K (60 s, temperature form)
                    comment = "the surface balance holds convection only",
                ),
            )
        end

        # The sub-surface peak at 90 s, and its cell at every save.
        j90 = column(90.0)
        T_peak, T_peak_ref = maximum(sol.T[:, j90]), maximum(ref.T[:, j90])
        push!(
            out,
            QoIResult(
                name = "sub-surface peak value, t=90 s",
                value = T_peak,
                exact = T_peak_ref,
                error = abs(T_peak - T_peak_ref),
                tol = 1.2e-6,  # 3× the largest observed T difference, 3.9e-7 K
                comment = "the coupled marker: in-depth absorption against surface loss and emission",
            ),
        )
        moved = count(j -> argmax(sol.T[:, j]) != argmax(ref.T[:, j]), 2:length(saves))
        push!(
            out,
            QoIResult(
                name = "peak cell equals the reference's at every save",
                value = Float64(moved),
                error = Float64(moved),
                tol = 0.5,
                comment = "number of saves (5 to 180 s) where the hottest cell differs",
            ),
        )
        exposed = count(j -> argmax(sol.T[:, j]) == n, findall(≥(10.0), saves))
        push!(
            out,
            QoIResult(
                name = "peak below the exposed cell from 10 s",
                value = Float64(exposed),
                error = Float64(exposed),
                tol = 0.5,
                comment = "number of saves from 10 s where the exposed cell is the hottest",
            ),
        )
        return out
    end

    function l2_error_t90(sol)
        z = sol.cell_centers[:, end]
        j = column(90.0)
        e = [sol.T[i, j] - continuum(z[i], 90.0) for i in eachindex(z)]
        return sqrt(sum(abs2, e) / length(e))
    end

    return VerificationCase(
        id = "V-RAD-07",
        title = "in-depth heating with emission (sub-surface peak)",
        tags = [:radiation, :radiative_transfer, :convective_bc, :coupled],
        references = [
            "reference/radiation/vrad07_reference.jl (independent method of lines)",
            "exact/two_stream.jl (the radiation rows, self-tested)",
            "Lautenberger dissertation (2007) §3.4.2 (the slab, without emission)",
        ],
        build = build,
        solve_kwargs = (
            abstol = 1e-10,
            reltol = 1e-8,
            radiation_model = RADIATIVE_TRANSFER,
        ),
        qois = qois,
        ci = CIConfig(
            n_cells = 100,
            budget_s = 2.0,
            convergence = ConvergenceSpec(
                # observed rates 2.004, 2.011 in both energy forms
                n_cells = [50, 100, 200],
                error_fn = l2_error_t90,
                expected_order = 2.0,
                min_order = 1.8,
                min_mean_order = 1.8,
                solve_kwargs = (
                    abstol = 1e-12,
                    reltol = 1e-10,
                    radiation_model = RADIATIVE_TRANSFER,
                ),
            ),
        ),
        report = (
            n_cells_ladder = [50, 100, 200],
            profile_times = times,
            qoi_curves = (:T_profiles, :surface_T_history, :peak_trajectory),
            solve_kwargs = (
                abstol = 1e-12,
                reltol = 1e-10,
                radiation_model = RADIATIVE_TRANSFER,
            ),
            curves = sol -> begin
                n = size(sol.T, 1)
                ref = references[n]
                z = sol.cell_centers[:, end]
                depth(T) = δ - z[argmax(T)]
                [
                    ReportCurve(
                        name = :T_profiles,
                        title = "V-RAD-07: in-depth heating with emission (sub-surface peak)",
                        xlabel = "depth x (m)",
                        ylabel = "T (K)",
                        series = [
                            CurveSeries(
                                label = "t = $(t) s",
                                x = δ .- z,
                                y = sol.T[:, column(t)],
                                exact = x -> continuum(δ - x, t),
                            ) for t in times
                        ],
                    ),
                    ReportCurve(
                        name = :surface_T_history,
                        title = "V-RAD-07: front-face temperature",
                        xlabel = "t (s)",
                        ylabel = "T_s (K)",
                        series = [
                            CurveSeries(label = "front face", x = sol.t, y = sol.surface_T),
                            CurveSeries(
                                label = "method-of-lines reference, same mesh",
                                x = ref.t,
                                y = ref.T_s,
                            ),
                        ],
                    ),
                    ReportCurve(
                        name = :peak_trajectory,
                        title = "V-RAD-07: depth of the hottest cell",
                        xlabel = "t (s)",
                        ylabel = "depth below the surface (m)",
                        series = [
                            CurveSeries(
                                label = "solver",
                                x = sol.t[2:end],
                                y = [depth(sol.T[:, j]) for j = 2:length(sol.t)],
                            ),
                            CurveSeries(
                                label = "method-of-lines reference, same mesh",
                                x = ref.t[2:end],
                                y = [depth(ref.T[:, j]) for j = 2:length(ref.t)],
                            ),
                        ],
                    ),
                ]
            end,
            error_norms = sol -> begin
                z = sol.cell_centers[:, end]
                e = [sol.T[i, column(90.0)] - continuum(z[i], 90.0) for i in eachindex(z)]
                (L2 = sqrt(sum(abs2, e) / length(e)), Linf = maximum(abs, e))
            end,
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vrad07)
