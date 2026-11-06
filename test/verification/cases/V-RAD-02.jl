# =============================================================================
# V-RAD-02 — a heater through a gray slab (radiation rows, closed forms)
# =============================================================================
#
# The radiation a heater sends through a gray slab of optical thickness
# τ_L = 2 on a graded mesh, with the exposed face's transmissivity ε = 0.8,
# behind it a mirror (an adiabatic back) or an open gray back facing its
# surroundings (`RadiativeBC` on `:bottom`). The answer is the solver's
# consistent start: `setup` solves the radiation rows for the initial state,
# and nothing is integrated in time.
#
# The heater's own field is found by superposition: the rows are linear in
# the fluxes for given temperatures, so the solve with the heater minus the
# solve without it, at the same uniform 300 K, is the field the heater alone
# sends through a cold slab (the solver cannot hold a cell at 0 K). With no
# component absorptivity, every heater enters at the interface's ε.
#
# - Collimated: the beam e^{−τ} at every face, and the thermal field the back
#   reflects into q⁺ (1 − ε_b of the arriving beam) and the exposed face
#   reflects down, in closed form (`two_stream_beam_slab`).
# - Diffuse: the closed form of the rows for a slab of no emission
#   (`two_stream_isothermal_slab` with B = 0).
# - A cone between 20° and 73°, carried by three beams at the Gauss–Legendre
#   nodes of its direction cosine: each beam is exactly e^{−τ/μ_h} of its
#   share, their sum is within 2×10⁻⁴ of the exact downward cone
#   (`gray_slab_cone_transmission`, black faces), and the whole field with
#   the beams converted at the back agrees with `two_stream_solve`.
#
# Every flux is compared at every face, as a fraction of the entering εq.
# Reported beside, not asserted: how far the collimated field behind a mirror
# is from exact radiative transfer with a specular mirror (`gray_slab_exact`),
# the distance the conversion at the back and the two-stream closure make.

const _VRAD_ρ = 1000.0

"""
    _vrad_nodes(L, n; kind = :uniform)

The node positions of a line of thickness `L` and `n` cells, from the back
(z = 0) to the exposed face: `:uniform`; `:graded`, widths growing by 10 %
per cell from the exposed face into the depth; `:irregular`, widths between
0.4 and 1.6 of the mean in a fixed pseudo-random order.
"""
function _vrad_nodes(L, n; kind::Symbol = :uniform)
    widths = if kind === :uniform
        ones(n)
    elseif kind === :graded
        [1.1^(n - i) for i = 1:n]
    elseif kind === :irregular
        [1.0 + 0.6 * sin(12.9898 * i + 78.233 * i^2) for i = 1:n]
    else
        throw(ArgumentError("unknown mesh kind $kind"))
    end
    nodes = [0.0; cumsum(widths)] .* (L / sum(widths))
    nodes[end] = L
    return nodes
end

"""
    _vrad_slab(nodes; κ, ε, T, top, bottom, back = nothing)

A gray slab on the line `nodes`: one solid of absorption coefficient `κ`
[1/m] and emissivity `ε` at the temperatures `T(z)` of its cell centres,
with the thermal conditions `top` and `bottom`. `back = (; z, ε)` puts a
second solid of the same κ and the emissivity `back.ε` in the cells below
`z`.
"""
function _vrad_slab(nodes; κ, ε, T, top, bottom, back = nothing)
    slab = SolidComponent(:slab, ρ = _VRAD_ρ, c = 1000.0, k = 1.0, ε = ε, α = κ / _VRAD_ρ)
    components, fields = (slab,), Any[z->_VRAD_ρ]
    if back !== nothing
        layer = SolidComponent(
            :layer,
            ρ = _VRAD_ρ,
            c = 1000.0,
            k = 1.0,
            ε = back.ε,
            α = κ / _VRAD_ρ,
        )
        components = (slab, layer)
        fields = Any[z->z>back.z ? _VRAD_ρ : 0.0, z->z>back.z ? 0.0 : _VRAD_ρ]
    end
    return PyrolysisProblem(
        mesh = line_mesh(nodes),
        material = Material(name = :VRAD, components = components, reactions = ()),
        bc_set = BoundaryConditionSet(thermal = (top = top, bottom = bottom)),
        T_initial = T,
        ξ_initial = fields,
        tspan = (0.0, 1.0),
    )
end

# The optical thicknesses of the cells of `nodes`, from the exposed face down,
# and the optical depth of each face.
function _vrad_depths(nodes, κ)
    τ = κ .* reverse(diff(nodes))
    return τ, [0.0; cumsum(τ)]
end

function _verification_case_vrad02()
    D = 1.5
    ε, κ, L, T = 0.8, 1000.0, 2e-3, 300.0
    q = 2.5e4
    θ1, θ2 = deg2rad(20.0), deg2rad(73.0)
    heaters = (
        collimated = RadiativeFluxBC(q; distribution = Collimated()),
        diffuse = RadiativeFluxBC(q),
        cone = RadiativeFluxBC(q; distribution = Cone(θ1, θ2; directions = 3)),
    )
    backs = (mirror = AdiabaticBC(), open = RadiativeBC(T_inf = T))
    surroundings = RadiativeBC(T_inf = T)

    deck(nodes, heater, back) = _vrad_slab(
        nodes;
        κ,
        ε,
        T = z -> T,
        top = heater === nothing ? surroundings : heater + surroundings,
        bottom = backs[back],
    )

    # The heater's field: with it minus without it, at the same temperatures.
    function heater_field(nodes, heater, back)
        with = radiation_fields(deck(nodes, heaters[heater], back))
        without = radiation_fields(deck(nodes, nothing, back))
        return (
            qm = with.qm .- without.qm,
            qp = with.qp .- without.qp,
            I = with.I,
            absorbed = with.absorbed - without.absorbed,
            emitted = with.emitted - without.emitted,
            back = with.back - without.back,
        )
    end

    # The closed forms on the slab of `nodes`, `back` open (at ε) or a mirror.
    function closed_form(nodes, heater, back)
        τ, depth = _vrad_depths(nodes, κ)
        τL, t = depth[end], exp(-D * depth[end])
        ε_b = back === :open ? ε : 0.0
        if heater === :diffuse
            c = two_stream_isothermal_slab(; B = 0.0, x = D * τL, ε, Q = q, ε_b)
            return (
                qm = c.qm_e .* exp.(-D .* depth),
                qp = c.qp_b .* exp.(-D .* (τL .- depth)),
                I = zeros(length(depth), 0),
                emitted = c.q_out,
                back = c.q_back,
            )
        end
        beams =
            heater === :collimated ? [TwoStreamBeam(1.0, ε * q)] :
            two_stream_cone_beams(q, θ1, θ2; a = ε)
        I = [b.q0 * exp(-d / b.μ) for d in depth, b in beams]
        arriving = sum(I[end, :])
        # The beams reaching the back: 1 − ε_b reflects into q⁺, which the
        # exposed face reflects down at 1 − ε (the rows' closed form).
        qp_e = t * (1 - ε_b) * arriving / (1 - (1 - ε) * (1 - ε_b) * t^2)
        qm = (1 - ε) * qp_e .* exp.(-D .* depth)
        qp = (qp_e / t) .* exp.(-D .* (τL .- depth))
        return (; qm, qp, I, emitted = ε * qp_e, back = qm[end] + arriving - qp[end])
    end

    build(n_cells) =
        deck(_vrad_nodes(L, n_cells; kind = :graded), heaters.collimated, :mirror)
    n_ci = 20

    function qois(_)
        nodes = _vrad_nodes(L, n_ci; kind = :graded)
        scale = ε * q
        out = QoIResult[]
        for heater in keys(heaters), back in keys(backs)
            f = heater_field(nodes, heater, back)
            c = closed_form(nodes, heater, back)
            err = max(
                maximum(abs, f.qm .- c.qm),
                maximum(abs, f.qp .- c.qp),
                maximum(abs, f.I .- c.I; init = 0.0),
                abs(f.absorbed - scale),
                abs(f.emitted - c.emitted),
                abs(f.back - c.back),
            )
            push!(
                out,
                QoIResult(
                    name = "$(heater) heater, $(back) back: every face against the closed form",
                    value = err / scale,
                    error = err / scale,
                    tol = 1e-12,
                    comment = "largest difference of q⁻, q⁺, the beams and the crossing fluxes, over εq",
                ),
            )
        end

        # The cone's three beams against the exact downward cone (black faces,
        # cold slab), and the whole field against the oracle's sweep.
        τ, depth = _vrad_depths(nodes, κ)
        f = heater_field(nodes, :cone, :mirror)
        exact = [gray_slab_cone_transmission(ε * q, d, θ1, θ2) for d in depth]
        err = maximum(abs, vec(sum(f.I; dims = 2)) .- exact)
        push!(
            out,
            QoIResult(
                name = "cone heater: three beams against the exact downward cone",
                value = err / scale,
                error = err / scale,
                tol = 2e-4,
                comment = "quadrature of the cone in three directions",
            ),
        )
        r = two_stream_solve(
            τ,
            zeros(length(τ));
            D,
            ε,
            ε_b = 0.0,
            beams = two_stream_cone_beams(q, θ1, θ2; a = ε),
        )
        err = max(
            maximum(abs, f.qm .- r.qm),
            maximum(abs, f.qp .- r.qp),
            maximum(abs, f.I .- r.Id),
        )
        push!(
            out,
            QoIResult(
                name = "cone heater, mirror back: the field against two_stream_solve",
                value = err / scale,
                error = err / scale,
                tol = 1e-12,
                comment = "the beams reaching the mirror join q⁺",
            ),
        )

        # Reported: the collimated field behind a mirror against exact
        # transfer with a specular mirror.
        f = heater_field(nodes, :collimated, :mirror)
        g = gray_slab_exact(τ, zeros(length(τ)); ε, ε_b = 0.0, q, heater = :collimated)
        S = two_stream_sources(f.qm, f.qp, f.I, zeros(size(f.I)))
        for (name, value) in (
            ("flux leaving the exposed face", abs(f.emitted - g.q_out)),
            ("L1 of the cell sources", sum(abs, S .- g.S)),
        )
            push!(
                out,
                QoIResult(
                    name = "collimated heater, mirror back: $(name) against exact transfer",
                    value = value / scale,
                    error = value / scale,
                    tol = Inf,
                    comment = "reported: the two-stream closure and the conversion at the back, against a specular mirror",
                ),
            )
        end
        return out
    end

    return VerificationCase(
        id = "V-RAD-02",
        title = "a heater through a gray slab (radiation rows)",
        tags = [:radiation, :radiative_transfer, :closed_form],
        references = [
            "exact/two_stream.jl (closed forms of the rows, self-tested)",
            "exact/gray_slab.jl (exact transfer; the cone in closed form)",
            "Modest & Mazumder, Radiative Heat Transfer, 4th ed., §13.3 and §14.3",
        ],
        build = build,
        solve_fn = radiation_fields,
        qois = qois,
        ci = CIConfig(n_cells = n_ci, budget_s = 5.0),
        report = (
            qoi_curves = (:collimated_mirror, :cone_beams),
            curves = _ -> begin
                nodes = _vrad_nodes(L, n_ci; kind = :graded)
                x = reverse(nodes)[1] .- reverse(nodes)   # depth of each face
                f = heater_field(nodes, :collimated, :mirror)
                c = closed_form(nodes, :collimated, :mirror)
                τL, scale = κ * L, ε * q
                qp_b = c.qp[end]
                cone = heater_field(nodes, :cone, :mirror)
                [
                    ReportCurve(
                        name = :collimated_mirror,
                        title = "V-RAD-02: a collimated heater through a slab, mirror back",
                        xlabel = "depth x (m)",
                        ylabel = "flux / εq",
                        series = [
                            CurveSeries(
                                label = "beam",
                                x = x,
                                y = vec(f.I) ./ scale,
                                exact = d -> exp(-κ * d),
                            ),
                            CurveSeries(
                                label = "q⁺ (from the mirror)",
                                x = x,
                                y = f.qp ./ scale,
                                exact = d -> qp_b * exp(-D * (τL - κ * d)) / scale,
                            ),
                            CurveSeries(
                                label = "q⁻ (reflected at the exposed face)",
                                x = x,
                                y = f.qm ./ scale,
                                exact = d -> c.qm[1] * exp(-D * κ * d) / scale,
                            ),
                        ],
                    ),
                    ReportCurve(
                        name = :cone_beams,
                        title = "V-RAD-02: a 20°–73° cone carried by three beams",
                        xlabel = "depth x (m)",
                        ylabel = "downward flux / εq",
                        series = [
                            CurveSeries(
                                label = "three beams",
                                x = x,
                                y = vec(sum(cone.I; dims = 2)) ./ scale,
                                exact = d ->
                                    gray_slab_cone_transmission(1.0, κ * d, θ1, θ2),
                            ),
                        ],
                    ),
                ]
            end,
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vrad02)
