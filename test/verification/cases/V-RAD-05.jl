# =============================================================================
# V-RAD-05 — the two-stream against exact transfer on burning profiles
# =============================================================================
#
# The gray two-stream (D = 3/2) is a closure: it carries the hemispherical
# fluxes as if the intensity were isotropic in each hemisphere. This case
# measures what that costs on the temperature and absorption profiles of a
# burning PMMA sample, against exact gray radiative transfer
# (`gray_slab_exact`: the transfer equation of a non-scattering gray medium
# solved per direction, with the same staircase emission and specular faces).
#
# The profiles are frames of the solver's own runs (the MaCFP-2 PMMA
# gasification target with the GIDAZE+ property set: 6 mm, 20 cells/mm, a
# collimated heater of 10 or 25 kW/m², surroundings at 293 K, no convection,
# a mirror back, moving mesh, deactivation) when 50 % and 10 % of the PMMA
# mass is left; reference/radiation/vrad05_profiles.jl writes them. The case
# rebuilds each frame as a deck, takes the radiation of its consistent start
# (nothing is integrated in time), and compares, as fractions of the
# heater's εq:
#
# - the thermal field alone (the solve without the heater: emission, the
#   surroundings, reflection at the faces): the flux leaving the exposed
#   face, the L1 difference of the cell sources over the top 0.5 mm, and the
#   exposed cell's source;
# - the whole field with the heater: the flux leaving and the L1 difference
#   of the sources over the whole sample. The beams themselves are exact; the
#   difference adds the conversion at the mirror (the beam reaching the back
#   travels on as diffuse);
# - the whole field with the same heater diffuse, the default distribution:
#   the same two measures. The two-stream attenuates diffuse incidence as
#   e^{−Dτ} where exact transfer gives 2E₃(τ): it transmits 3 % more through
#   an optical depth of 0.1 and 38 % less through one of 3.
#
# The values are the model's error, not a discretization error: they are
# pinned at 3× the observed as a regression. The exposed cell's two-stream
# source is about three quarters of the exact one (the thin-layer limit, Dτ
# against 2τ). The rebuilt frames are checked against their files: the
# solver's absorption coefficient is ξα of the stored PMMA concentration and
# the exposed cell's emissivity PMMA's 0.96.
#
# GIDAZE+ (MaCFP matl-db, PMMA/Material_Properties/2021/MaCFP_PMMA_GIDAZE+.json):
# ρ = 1210 kg/m³, ε = 0.96, κ = 2870 m⁻¹, k = 0.2(T/300)^−1.2 W/m·K,
# c = 1200(T/300)^1.2 J/kg·K, PMMA → MMA first order with A = 2×10⁸ s⁻¹,
# E = 130 kJ/mol, ΔH = 540 kJ/kg; the gas carries no heat capacity.

const _VRAD_GIDAZE_T0 = 293.0
const _VRAD_GIDAZE_ρ = 1210.0
const _VRAD_GIDAZE_κ = 2870.0

function _vrad_gidaze_material()
    pmma = SolidComponent(
        :PMMA;
        ρ = _VRAD_GIDAZE_ρ,
        φ = 0.01,
        c = CallableProperty(T -> 1200.0 * (T / 300.0)^1.2),
        k = CallableProperty(T -> 0.2 * (T / 300.0)^-1.2),
        ε = 0.96,
        α = _VRAD_GIDAZE_κ / _VRAD_GIDAZE_ρ,
        λ = 2e-5,
    )
    gas = GaseousComponent(:gas; M = 0.10012, c = 0.0, k = 0.02)
    reaction = Reaction(:r1, 1 => 2; A = 2.0e8, E = 1.3e5, h = 5.4e5, n = 1.0)
    return Material(
        name = :GIDAZE,
        components = (pmma, gas),
        reactions = (reaction,),
        mixing = MaterialMixing(k = PARALLEL, λ = PARALLEL),
    )
end

"""
    _vrad_gidaze_problem(nodes, T, ξ; q, distribution = Collimated(),
                         bottom = AdiabaticBC(), tspan, output_times)

The GIDAZE+ PMMA on the line `nodes` (back to exposed face) at the cell
temperatures `T(z)` and PMMA concentrations `ξ(z)`, under a heater `q` [W/m²]
arriving as `distribution` (none when zero) beside surroundings at 293 K, the
gas leaving through the exposed face to the ambient and not through the back.
"""
function _vrad_gidaze_problem(
    nodes,
    T,
    ξ;
    q,
    distribution = Collimated(),
    bottom = AdiabaticBC(),
    tspan = (0.0, 1.0),
    output_times = Float64[],
)
    surroundings = RadiativeBC(T_inf = _VRAD_GIDAZE_T0)
    top = q > 0 ? RadiativeFluxBC(q; distribution) + surroundings : surroundings
    return PyrolysisProblem(
        mesh = line_mesh(nodes),
        material = _vrad_gidaze_material(),
        bc_set = BoundaryConditionSet(
            thermal = (top = top, bottom = bottom),
            mass = (
                top = AmbientGasBC(GasReservoir(P = 0.0, T = _VRAD_GIDAZE_T0)),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = T,
        ξ_initial = [ξ, z -> 0.0],
        tspan = tspan,
        dt_initial = 1e-3,
        output_times = output_times,
    )
end

# A stored frame: the node positions and, per cell from the back, T and the
# PMMA concentration.
function _vrad05_frame(q_kW, left)
    names, cols = read_reference_csv("radiation/vrad05_q$(q_kW)_left$(left).csv")
    z_lower, z_upper, T, ξ = cols
    return (nodes = [z_lower[1]; z_upper], T, ξ)
end

function _verification_case_vrad05()
    D = 1.5
    ε = 0.96
    B_amb = STEFAN_BOLTZMANN * _VRAD_GIDAZE_T0^4
    frames = [(q_kW, left) for q_kW in (10, 25) for left in (50, 10)]
    # The observed model errors over εq, identical in both energy forms; each
    # quantity is pinned at three times its value.
    observed = Dict(
        (10, 50) => (
            leaving = 3.98e-4,
            top = 0.0626,
            exposed = 0.0105,
            heated_leaving = 3.98e-4,
            heated = 0.0732,
            diffuse_leaving = 3.98e-4,
            diffuse = 0.0237,
        ),
        (10, 10) => (
            leaving = 4.43e-3,
            top = 0.0676,
            exposed = 4.89e-3,
            heated_leaving = 0.0135,
            heated = 0.0811,
            diffuse_leaving = 2.75e-4,
            diffuse = 0.0124,
        ),
        (25, 50) => (
            leaving = 3.94e-4,
            top = 0.0387,
            exposed = 2.24e-3,
            heated_leaving = 3.94e-4,
            heated = 0.0457,
            diffuse_leaving = 3.94e-4,
            diffuse = 0.0512,
        ),
        (25, 10) => (
            leaving = 2.30e-3,
            top = 0.0393,
            exposed = 6.73e-4,
            heated_leaving = 0.0158,
            heated = 0.0582,
            diffuse_leaving = 2.44e-3,
            diffuse = 0.0422,
        ),
    )

    function deck(frame, q; distribution = Collimated())
        nodes = frame.nodes
        cell(z) = clamp(searchsortedlast(nodes, z), 1, length(frame.T))
        return _vrad_gidaze_problem(
            nodes,
            z -> frame.T[cell(z)],
            z -> frame.ξ[cell(z)];
            q,
            distribution,
        )
    end

    build(n_cells) = deck(_vrad05_frame(25, 10), 25e3)

    function qois(_)
        out = QoIResult[]
        for (q_kW, left) in frames
            frame = _vrad05_frame(q_kW, left)
            q = 1e3 * q_kW
            scale = ε * q
            label = "$(q_kW) kW/m², $(left) % left"
            Δ = reverse(diff(frame.nodes))               # from the exposed face
            κ = _VRAD_GIDAZE_κ / _VRAD_GIDAZE_ρ .* reverse(frame.ξ)
            τ = κ .* Δ
            B = STEFAN_BOLTZMANN .* reverse(frame.T) .^ 4
            top = cumsum(Δ) .- Δ .< 0.5e-3               # cells starting in the top 0.5 mm

            thermal = radiation_fields(deck(frame, 0.0))
            heated = radiation_fields(deck(frame, q))
            consistency =
                max(maximum(abs, thermal.κ .- κ) / _VRAD_GIDAZE_κ, abs(thermal.ε[1] - ε))
            push!(
                out,
                QoIResult(
                    name = "$(label): the rebuilt frame's κ and exposed ε",
                    value = consistency,
                    error = consistency,
                    tol = 1e-12,
                    comment = "κ = ξα of the stored PMMA concentration; ε of PMMA",
                ),
            )

            g = gray_slab_exact(τ, B; ε, B_amb, ε_b = 0.0)
            S = two_stream_sources(
                thermal.qm,
                thermal.qp,
                thermal.I,
                zeros(size(thermal.I)),
            )
            pinned = observed[(q_kW, left)]
            for (name, value, seen) in (
                (
                    "thermal field: flux leaving",
                    abs(thermal.emitted - g.q_out),
                    pinned.leaving,
                ),
                (
                    "thermal field: L1 of the sources, top 0.5 mm",
                    sum(abs, (S .- g.S)[top]),
                    pinned.top,
                ),
                (
                    "thermal field: the exposed cell's source",
                    abs(S[1] - g.S[1]),
                    pinned.exposed,
                ),
            )
                push!(
                    out,
                    QoIResult(
                        name = "$(label), $(name)",
                        value = value / scale,
                        error = value / scale,
                        tol = 3 * seen,
                        comment = "two-stream against exact transfer, over εq; exposed-cell source ratio $(round(S[1] / g.S[1], digits = 3))",
                    ),
                )
            end

            g = gray_slab_exact(τ, B; ε, B_amb, ε_b = 0.0, q, heater = :collimated)
            S = two_stream_sources(heated.qm, heated.qp, heated.I, zeros(size(heated.I)))
            for (name, value, seen) in (
                (
                    "with the heater: flux leaving",
                    abs(heated.emitted - g.q_out),
                    pinned.heated_leaving,
                ),
                ("with the heater: L1 of the sources", sum(abs, S .- g.S), pinned.heated),
            )
                push!(
                    out,
                    QoIResult(
                        name = "$(label), $(name)",
                        value = value / scale,
                        error = value / scale,
                        tol = 3 * seen,
                        comment = "two-stream with conversion at the mirror against exact transfer, over εq",
                    ),
                )
            end

            diffuse = radiation_fields(deck(frame, q; distribution = Diffuse()))
            g = gray_slab_exact(τ, B; ε, B_amb, ε_b = 0.0, q, heater = :diffuse)
            S = two_stream_sources(
                diffuse.qm,
                diffuse.qp,
                diffuse.I,
                zeros(size(diffuse.I)),
            )
            for (name, value, seen) in (
                (
                    "with a diffuse heater: flux leaving",
                    abs(diffuse.emitted - g.q_out),
                    pinned.diffuse_leaving,
                ),
                (
                    "with a diffuse heater: L1 of the sources",
                    sum(abs, S .- g.S),
                    pinned.diffuse,
                ),
            )
                push!(
                    out,
                    QoIResult(
                        name = "$(label), $(name)",
                        value = value / scale,
                        error = value / scale,
                        tol = 3 * seen,
                        comment = "two-stream, the heater attenuated as e^{−Dτ}, against exact transfer, over εq",
                    ),
                )
            end
        end
        return out
    end

    return VerificationCase(
        id = "V-RAD-05",
        title = "the two-stream against exact transfer on burning profiles",
        tags = [:radiation, :radiative_transfer, :model_error],
        references = [
            "exact/gray_slab.jl (exact gray transfer, self-tested)",
            "reference/radiation/vrad05_profiles.jl (the frames)",
            "MaCFP matl-db, PMMA/Material_Properties/2021/MaCFP_PMMA_GIDAZE+.json",
        ],
        build = build,
        solve_fn = radiation_fields,
        qois = qois,
        ci = CIConfig(n_cells = 1, budget_s = 10.0),
        report = (
            qoi_curves = (:sources_q10_left50, :sources_q25_left10),
            curves = _ -> map(((10, 50), (25, 10))) do (q_kW, left)
                frame = _vrad05_frame(q_kW, left)
                Δ = reverse(diff(frame.nodes))
                depth = cumsum(Δ) .- Δ ./ 2
                τ = _VRAD_GIDAZE_κ / _VRAD_GIDAZE_ρ .* reverse(frame.ξ) .* Δ
                B = STEFAN_BOLTZMANN .* reverse(frame.T) .^ 4
                f = radiation_fields(deck(frame, 0.0))
                S = two_stream_sources(f.qm, f.qp, f.I, zeros(size(f.I))) ./ Δ
                g = gray_slab_exact(τ, B; ε, B_amb, ε_b = 0.0).S ./ Δ
                ReportCurve(
                    name = Symbol("sources_q$(q_kW)_left$(left)"),
                    title = "V-RAD-05: thermal-field sources, $(q_kW) kW/m², $(left) % of the mass left",
                    xlabel = "depth x (m)",
                    ylabel = "radiative source (W/m³)",
                    series = [
                        CurveSeries(
                            label = "two-stream",
                            x = depth,
                            y = S,
                            exact = x ->
                                g[clamp(searchsortedfirst(cumsum(Δ), x), 1, length(g))],
                        ),
                    ],
                )
            end,
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vrad05)
