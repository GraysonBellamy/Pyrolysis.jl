# =============================================================================
# V-RAD-01 — a gray slab in a radiative enclosure at its own temperature
# =============================================================================
#
# Under radiative transfer every cell emits and absorbs with one absorption
# coefficient, and a face crosses radiation at its cell's emissivity in both
# directions, so a slab at T_h that faces surroundings at T_h exchanges
# nothing: at any optical thickness, on any mesh, behind a mirror or an open
# back.
#
# The rows: a slab 1 mm thick at a uniform 700 K whose exposed face looks at
# surroundings at 700 K (`RadiativeBC`), behind a mirror (an adiabatic back) or
# an open back facing 700 K, with κL from 10⁻³ to 10³, the emissivity 0.3,
# 0.96 or 1, on uniform, graded and irregular meshes of 16 cells. Every face
# flux is σT_h⁴ and every cell's radiative source zero, to 10⁻¹² σT_h⁴. The
# answer is the solver's consistent start; nothing is integrated in time.
#
# The solver: a slab of 40 cells with PMMA's density, absorption coefficient
# and emissivity, κL from 0.03 to 10, integrated from T_h to 30 of its
# radiative time constants ρcL/(4εσT_h³) in three enclosures:
#
# - E1: diffuse radiation σT_h⁴ onto the exposed face, an open back facing T_h;
# - E2: the same exposed face, a mirror back;
# - E3: as E2, behind an opaque 10 µm layer of another material.
#
# The exposed face's radiation is a diffuse `RadiativeFluxBC(σT_h⁴)` beside a
# `RadiativeBC` facing 1 K, so the flux and the surroundings both enter the
# rows; the 1 K surroundings add σ·1⁴, which settles the slab 1/(4T_h³) =
# 7.3×10⁻¹⁰ K above T_h behind a mirror and half that with E1's open back.
# Every cell must stay within 10⁻⁸ K of T_h at every saved time.

# The slabs, meshes and closed forms below are shared by V-RAD-02 to V-RAD-07
# (their files load after this one).

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

# The enclosure decks of the solver part: `case` E1, E2 or E3 at optical
# thickness `κL` on `n` cells, saved at the start and at a decade of times
# up to 30 radiative time constants.
function _vrad01_enclosure(κL, case; n = 40, T_h = 700.0)
    ρ, c, κ, ε = 1210.0, 1500.0, 2870.0, 0.96
    L = κL / κ
    slab = SolidComponent(:slab, ρ = ρ, c = c, k = 0.2, ε = ε, α = κ / ρ)
    components = (slab,)
    nodes = collect(range(0.0, L; length = n + 1))
    fields = Any[z->ρ]
    if case === :E3
        back_layer = 1e-5
        back = SolidComponent(:back, ρ = 1000.0, c = 1000.0, k = 1.0, ε = ε, α = 1e6)
        components = (slab, back)
        nodes = [0.0; back_layer .+ nodes]
        inside(z) = z > back_layer - 1e-12
        fields = Any[z->inside(z) ? ρ : 0.0, z->inside(z) ? 0.0 : 1000.0]
    end
    top = RadiativeFluxBC(STEFAN_BOLTZMANN * T_h^4) + RadiativeBC(T_inf = 1.0)
    bottom = case === :E1 ? RadiativeBC(T_inf = T_h) : AdiabaticBC()
    t_end = 30 * ρ * c * max(L, 1e-4) / (4ε * STEFAN_BOLTZMANN * T_h^3)
    return PyrolysisProblem(
        mesh = line_mesh(nodes),
        material = Material(name = :enclosure, components = components, reactions = ()),
        bc_set = BoundaryConditionSet(thermal = (top = top, bottom = bottom)),
        T_initial = z -> T_h,
        ξ_initial = fields,
        tspan = (0.0, t_end),
        output_times = [0.0; t_end .* exp10.(-6:0)],
    )
end

function _verification_case_vrad01()
    T_h, L, n_rows = 700.0, 1e-3, 16
    B = STEFAN_BOLTZMANN * T_h^4
    surroundings = RadiativeBC(T_inf = T_h)
    backs = (mirror = AdiabaticBC(), open = surroundings)
    meshes = (:uniform, :graded, :irregular)
    emissivities = (0.3, 0.96, 1.0)
    thicknesses = (1e-3, 1e-1, 1.0, 10.0, 1e3)   # κL
    enclosures = (:E1, :E2, :E3)
    depths = (0.03, 0.1, 0.3, 1.0, 3.0, 10.0)   # κL
    kwargs = (radiation_model = RADIATIVE_TRANSFER,)

    # The largest departure of a face flux from σT_h⁴, of a cell's radiative
    # source from zero, and of the radiation crossing the faces from zero,
    # over σT_h⁴.
    function row_error(nodes, κL, ε, back)
        f = radiation_fields(
            _vrad_slab(
                nodes;
                κ = κL / L,
                ε,
                T = z -> T_h,
                top = surroundings,
                bottom = backs[back],
            ),
        )
        net = f.qm .- f.qp
        source = net[1:(end-1)] .- net[2:end]
        return max(
            maximum(abs, f.qm .- B),
            maximum(abs, f.qp .- B),
            maximum(abs, source),
            abs(f.absorbed - f.emitted - f.back),
            abs(f.back),
        ) / B
    end

    deviation(sol) = maximum(abs, sol.T .- T_h)
    function deviations(case, main)
        return [
            (case === :E2 && κL == 1.0) ? deviation(main) :
            deviation(suite_solve(_vrad01_enclosure(κL, case); kwargs...)) for κL in depths
        ]
    end

    function qois(sol)
        out = QoIResult[]
        for back in keys(backs), kind in meshes
            nodes = _vrad_nodes(L, n_rows; kind)
            err = maximum(
                row_error(nodes, κL, ε, back) for κL in thicknesses, ε in emissivities
            )
            push!(
                out,
                QoIResult(
                    name = "rows, $(back) back, $(kind) mesh: every face flux σT_h⁴, every source zero",
                    value = err,
                    error = err,
                    tol = 1e-12,
                    comment = "largest over κL = 1e-3 … 1e3 and ε = 0.3, 0.96, 1, as a fraction of σT_h⁴",
                ),
            )
        end
        for case in enclosures
            d = maximum(deviations(case, sol))
            push!(
                out,
                QoIResult(
                    name = "solver, $(case): every cell stays at T_h",
                    value = d,
                    error = d,
                    tol = 1e-8,
                    comment = "K, largest over κL = 0.03 … 10, every cell and saved time; the 1 K surroundings settle the slab 7.3e-10 K above T_h",
                ),
            )
        end
        return out
    end

    return VerificationCase(
        id = "V-RAD-01",
        title = "a gray slab in a radiative enclosure at its own temperature",
        tags = [:radiation, :radiative_transfer, :equilibrium],
        references = [
            "Kirchhoff's law: the cells absorb and emit with one coefficient, the faces transmit with one emissivity both ways",
        ],
        build = n -> _vrad01_enclosure(1.0, :E2; n),
        solve_kwargs = kwargs,
        qois = qois,
        ci = CIConfig(n_cells = 40, budget_s = 10.0),
        report = (
            qoi_curves = (:enclosure_deviation,),
            curves = sol -> [
                ReportCurve(
                    name = :enclosure_deviation,
                    title = "V-RAD-01: largest departure from T_h",
                    xlabel = "κL",
                    ylabel = "max |T − T_h| (K)",
                    series = [
                        CurveSeries(
                            label = "$(case)",
                            x = collect(depths),
                            y = deviations(case, sol),
                        ) for case in enclosures
                    ],
                    xscale = :log10,
                ),
            ],
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vrad01)
