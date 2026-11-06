# =============================================================================
# V-RAD-03 — an isothermal slab emitting to 0 K (radiation rows, closed forms)
# =============================================================================
#
# A gray slab 1 mm thick at a uniform 900 K with the exposed face's
# transmissivity ε = 0.8 faces surroundings at 0 K (a diffuse incoming flux of
# zero) and loses its emission through that face. Behind it:
#
# - a mirror (an adiabatic back): εB(1 − t²)/(1 − (1 − ε)t²);
# - an open back of the same material, facing 0 K:
#   εB(1 − t)/(1 − (1 − ε)t);
# - an open back whose last cell is another material, of the same absorption
#   coefficient and emissivity ε_l = 0.3:
#   εB(1 − t)(1 + (1 − ε_l)t)/(1 − (1 − ε)(1 − ε_l)t²);
#
# with t = e^{−Dτ_L} and D = 3/2 (`two_stream_slab_emission`). The rows are
# exact across each cell at any optical thickness, so every face flux
# (`two_stream_isothermal_slab` gives the four face values, and each face is
# B + (q − B)e^{−Dτ} of them) is independent of the mesh: uniform, graded and
# irregular meshes of 16 cells, Dτ_L from 10⁻⁴ to 10. The answer is the
# solver's consistent start; nothing is integrated in time.
#
# Reported beside, not asserted: the two-stream's thin-slab deficit. With
# black faces and an open back the slab emits B(1 − e^{−Dτ_L}) against the
# exact B(1 − 2E₃(τ_L)), 25 % low as τ_L → 0 (Dτ against 2τ) and exact as
# the slab becomes opaque.

function _verification_case_vrad03()
    D = 1.5
    ε, ε_l, L, T = 0.8, 0.3, 1e-3, 900.0
    n = 16
    B = STEFAN_BOLTZMANN * T^4
    cold = RadiativeFluxBC(0.0)
    backs = (
        mirror = (bottom = AdiabaticBC(), layer = nothing, ε_b = 0.0),
        open = (bottom = cold, layer = nothing, ε_b = ε),
        open_layer = (bottom = cold, layer = ε_l, ε_b = ε_l),
    )
    meshes = (:uniform, :graded, :irregular)
    thicknesses = (1e-4, 1e-2, 0.3, 1.0, 10.0)   # Dτ_L

    function deck(nodes, x, back; ε = ε, ε_layer = backs[back].layer)
        layer = ε_layer === nothing ? nothing : (; z = nodes[2], ε = ε_layer)
        return _vrad_slab(
            nodes;
            κ = x / (D * L),
            ε,
            T = z -> T,
            top = cold,
            bottom = backs[back].bottom,
            back = layer,
        )
    end

    # The largest difference of every face flux and of the fluxes leaving
    # from the closed forms, over B.
    function face_error(nodes, x, back)
        f = radiation_fields(deck(nodes, x, back))
        ε_b = backs[back].ε_b
        c = two_stream_isothermal_slab(; B, x, ε, ε_b)
        _, depth = _vrad_depths(nodes, x / (D * L))
        qm = B .+ (c.qm_e - B) .* exp.(-D .* depth)
        qp = B .+ (c.qp_b - B) .* exp.(-D .* (depth[end] .- depth))
        return max(
            maximum(abs, f.qm .- qm),
            maximum(abs, f.qp .- qp),
            abs(f.emitted - two_stream_slab_emission(B, x, ε, ε_b)),
            abs(f.back - c.q_back),
            abs(f.absorbed),
        ) / B
    end

    build(n_cells) = deck(_vrad_nodes(L, n_cells), 1.0, :mirror)

    function qois(_)
        out = QoIResult[]
        for back in keys(backs), kind in meshes
            nodes = _vrad_nodes(L, n; kind)
            err = maximum(face_error(nodes, x, back) for x in thicknesses)
            push!(
                out,
                QoIResult(
                    name = "$(back) back, $(kind) mesh: every face against the closed form",
                    value = err,
                    error = err,
                    tol = 1e-12,
                    comment = "largest over Dτ_L = 1e-4 … 10, as a fraction of B",
                ),
            )
        end
        # Reported: the thin-slab deficit, black faces, open back.
        nodes = _vrad_nodes(L, n)
        for x in thicknesses
            f = radiation_fields(deck(nodes, x, :open; ε = 1.0))
            τL = x / D
            exact = -gray_slab_black_faces(fill(τL / n, n), fill(B, n))[1]
            push!(
                out,
                QoIResult(
                    name = "black faces, Dτ_L = $(x): two-stream emission over exact",
                    value = f.emitted / exact,
                    exact = 1.0,
                    error = abs(f.emitted / exact - 1),
                    tol = Inf,
                    comment = "reported: B(1 − e^{−Dτ_L}) against B(1 − 2E₃(τ_L))",
                ),
            )
        end
        return out
    end

    return VerificationCase(
        id = "V-RAD-03",
        title = "an isothermal slab emitting to 0 K (radiation rows)",
        tags = [:radiation, :radiative_transfer, :closed_form],
        references = [
            "exact/two_stream.jl (closed forms of the rows, self-tested)",
            "exact/gray_slab.jl (exact transfer with black faces)",
            "Modest & Mazumder, Radiative Heat Transfer, 4th ed., §14.3",
        ],
        build = build,
        solve_fn = radiation_fields,
        qois = qois,
        ci = CIConfig(n_cells = n, budget_s = 5.0),
        report = (
            qoi_curves = (:slab_emission,),
            curves = _ -> begin
                xs = exp10.(range(-4, 1, length = 21))
                nodes = _vrad_nodes(L, n; kind = :irregular)
                emitted(x, back; kw...) =
                    radiation_fields(deck(nodes, x, back; kw...)).emitted / B
                series = [
                    CurveSeries(
                        label = "$(back) back",
                        x = xs,
                        y = [emitted(x, back) for x in xs],
                        exact = x ->
                            two_stream_slab_emission(B, x, ε, backs[back].ε_b) / B,
                    ) for back in keys(backs)
                ]
                push!(
                    series,
                    CurveSeries(
                        label = "black faces, open back, against exact transfer",
                        x = xs,
                        y = [emitted(x, :open; ε = 1.0) for x in xs],
                        exact = x -> -gray_slab_black_faces([x / D], [1.0])[1],
                    ),
                )
                [
                    ReportCurve(
                        name = :slab_emission,
                        title = "V-RAD-03: emission of an isothermal slab into 0 K",
                        xlabel = "Dτ_L",
                        ylabel = "emitted / σT⁴",
                        series = series,
                        xscale = :log10,
                        yscale = :log10,
                    ),
                ]
            end,
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vrad03)
