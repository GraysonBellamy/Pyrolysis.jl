# =============================================================================
# V-RAD-04 — a linear emissive power (radiation rows, closed form)
# =============================================================================
#
# A black slab (ε = 1 at both faces, 40 equal cells) whose emissive power
# σT⁴ falls by g = 50 W/m² per cell into the depth, 2×10⁴ − g W/m² in the
# exposed cell (T between 751 and 770 K). Diffuse fluxes arrive on both faces
# at the values of the rows' particular solution, q± = B ∓ const, so every
# face carries the net flux
#
#   q⁻ − q⁺ = g·coth(Dτ/2),   D = 3/2, τ the cell's optical thickness
#
# (`two_stream_linear_flux`), at any τ: the rows are exact across a cell. The
# answer is the solver's consistent start; nothing is integrated in time.
#
# Reported beside, not asserted: the approach to the Rosseland flux. The
# continuous two-stream carries −(2/(Dκ))dB/dx, which at D = 3/2 is the
# diffusion limit; per cell of Dτ = x the rows carry (x/2)coth(x/2) of it,
# 1 + x²/12 for small x. An optically thick cell therefore adds a spurious
# radiative conductivity (×2.53 at Dτ = 5 per cell).

function _verification_case_vrad04()
    D = 1.5
    n, L = 40, 1e-3
    B0, g = 2e4, 50.0
    Δ = L / n
    per_cell = (0.01, 0.1, 1.0, 5.0)   # Dτ of each cell

    # B at the centre of the cell i cells from the exposed face is B0 − g·i.
    B_at(z) = B0 - g * ((L - z) / Δ + 0.5)

    function deck(n_cells, x)
        h = -expm1(-x)
        return _vrad_slab(
            _vrad_nodes(L, n_cells);
            κ = x / (D * Δ),
            ε = 1.0,
            T = z -> (B_at(z) / STEFAN_BOLTZMANN)^0.25,
            top = RadiativeFluxBC(B0 - g + g / h),
            bottom = RadiativeFluxBC(B0 - g * n_cells - g / h),
        )
    end

    build(n_cells) = deck(n_cells, 1.0)

    function qois(_)
        out = QoIResult[]
        for x in per_cell
            f = radiation_fields(deck(n, x))
            err = maximum(abs, (f.qm .- f.qp) .- two_stream_linear_flux(g, x)) / B0
            push!(
                out,
                QoIResult(
                    name = "Dτ = $(x) per cell: the net flux at every face",
                    value = err,
                    error = err,
                    tol = 1e-12,
                    comment = "g·coth(Dτ/2), as a fraction of the largest emissive power",
                ),
            )
        end
        for x in per_cell
            ratio = two_stream_linear_flux(g, x) / (2g / x) - 1
            f = radiation_fields(deck(n, x))
            measured = (f.qm[n÷2] - f.qp[n÷2]) / (2g / x) - 1
            push!(
                out,
                QoIResult(
                    name = "Dτ = $(x) per cell: the flux over Rosseland's, less one",
                    value = measured,
                    exact = x^2 / 12,
                    error = abs(measured - ratio),
                    tol = Inf,
                    comment = "reported: (x/2)coth(x/2) − 1, against x²/12 for small x",
                ),
            )
        end
        return out
    end

    return VerificationCase(
        id = "V-RAD-04",
        title = "a linear emissive power (radiation rows)",
        tags = [:radiation, :radiative_transfer, :closed_form],
        references = [
            "exact/two_stream.jl (`two_stream_linear_flux`, self-tested)",
            "Modest & Mazumder, Radiative Heat Transfer, 4th ed., §15.6 (Rosseland)",
        ],
        build = build,
        solve_fn = radiation_fields,
        qois = qois,
        ci = CIConfig(n_cells = n, budget_s = 5.0),
        report = (
            qoi_curves = (:rosseland,),
            curves = _ -> begin
                xs = exp10.(range(-2, log10(5), length = 15))
                excess(x) = begin
                    f = radiation_fields(deck(n, x))
                    (f.qm[n÷2] - f.qp[n÷2]) / (2g / x) - 1
                end
                [
                    ReportCurve(
                        name = :rosseland,
                        title = "V-RAD-04: the rows' flux over the Rosseland flux, less one",
                        xlabel = "Dτ per cell",
                        ylabel = "q / q_Rosseland − 1",
                        series = [
                            CurveSeries(
                                label = "rows",
                                x = xs,
                                y = excess.(xs),
                                exact = x -> (x / 2) * coth(x / 2) - 1,
                            ),
                        ],
                        xscale = :log10,
                        yscale = :log10,
                    ),
                ]
            end,
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vrad04)
