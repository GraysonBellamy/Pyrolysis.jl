# Leftover 1 of the Phase 1 plan (§7.1): the lifted assembly's measured cost.
#
# On the audit decks (`eulerian`, `ale`, `ale_darcy` at 25, 100 and 400
# cells: `problem_for` of the Phase 0 harness), one colored Jacobian of the
# finite-volume core against one `Structured` assembly of the old core
# (`DirectSolve` and `ApproxSolve`), warm medians, in one session. Each core
# is evaluated at the same physical state, the harness's mixed reacting
# profile (`hot_state!`): temperatures from 450 to 650 K, the solid reduced
# by 5 to 35 %, a gas concentration from 0.02 to 0.1 kg/m³.
#
# The two factors of Phase 0's projection are reported separately:
#
#   colors          the new core's color count against the lifted model's
#                   count at c = 0 (`lifted_pattern`, Phase 0 item 4);
#   per color       one dual partial in residual-equivalents,
#                   (t_dual(C) − t_res) / (C t_res), on each core's residual
#                   at its batch width (the old core at C = 16, as Phase 0),
#                   and the whole colored Jacobian per color, t_jac / (K t_res).
#
# Phase 0's projection, ⌈K_lifted / 16⌉ t_dual,old(16), is recomputed in the
# same session beside the measurement. The new core runs in both energy forms.
#
# Run from the repository root when the machine is least loaded:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/leftover1_jacobian_cost.jl
# Writes test/newcore/gates/measurements/leftover1_jacobian_cost.{csv,md}.

const _HARNESS =
    normpath(joinpath(@__DIR__, "..", "..", "..", "scripts", "mesh_phase0", "baselines"))
include(joinpath(_HARNESS, "jacobian_common.jl"))

using ForwardDiff
using ForwardDiff: Dual, Partials
using Pyrolysis.States: n_cells, n_components, mass_index, energy_index, energy_state
using Pyrolysis.States: energy_state_value
using Pyrolysis.ColoredJacobian: jac_prototype, jacobian!
using Pyrolysis.Integration: setup
const new_residual! = Pyrolysis.Assembly.residual!

const OUT = joinpath(@__DIR__, "measurements")
mkpath(OUT)
const DECKS = (:eulerian, :ale, :ale_darcy)
const SIZES = (25, 100, 400)
const ENERGY_STATES = (:temperature, :enthalpy)
const T_EVAL = 1.0
const OLD_CHUNK = 16
const REPEATS = 11

struct OldTag end

# One-hot partials in slot `mod1(i, C)`, as the Phase 0 dual timing seeds them.
function _seeded(::Type{D}, u) where {D}
    C = ForwardDiff.npartials(D)
    return D[
        D(u[i], Partials{C,Float64}(ntuple(m -> m == mod1(i, C) ? 1.0 : 0.0, Val(C)))) for
        i in eachindex(u)
    ]
end

# The new core's state at the harness's mixed reacting profile: the masses
# from the concentrations and the reference volumes, the energy state from the
# temperature; positions, velocities and intensities as set up.
function _hot_state(cp)
    layout = cp.def.layout
    n = n_cells(layout)
    V = cp.def.setup.reference.volume
    m0 = cp.m0
    u = copy(cp.u0)
    for i = 1:n
        ξ_solid = m0[1, i] / V[i] * (0.65 + 0.3 * i / n)
        m = (ξ_solid * V[i], (0.02 + 0.08 * i / n) * V[i])
        u[mass_index(layout, i, 1)] = m[1]
        u[mass_index(layout, i, 2)] = m[2]
        u[energy_index(layout, i)] = energy_state_value(
            energy_state(layout),
            cp.def.setup,
            450.0 + 200.0 * i / n,
            m,
            i,
        )
    end
    return u
end

function _old_core(deck, n)
    def, ws, u, NC, N_cells = build_deck(deck, n)
    hot_state!(u, N_cells)
    N = length(u)
    du = similar(u)
    t_res = time_op(() -> residual!(du, u, def, T_EVAL, ws); repeats = REPEATS)
    structured = map((DirectSolve(), ApproxSolve())) do strategy
        spec = JacobianSpec(Structured(; strategy), def, ws)
        J = Pyrolysis.Jacobian._alloc_target(spec, N)
        compute_jacobian!(J, spec, u, def, T_EVAL, ws)
        all(isfinite, nonzeros(J)) ||
            error("non-finite old-core Jacobian on $deck at $n cells")
        time_op(() -> compute_jacobian!(J, spec, u, def, T_EVAL, ws); repeats = REPEATS)
    end
    Tu = Dual{OldTag,Float64,OLD_CHUNK}
    ws_ad = build_ad_workspace(def, ws, Tu)
    ws_ad.I_radiation = ws.I_radiation
    u_d = _seeded(Tu, u)
    du_d = similar(u_d)
    t_dual = time_op(() -> residual!(du_d, u_d, def, T_EVAL, ws_ad); repeats = REPEATS)
    mode = def.mode
    chi = chi_mode(mode) isa WithChi
    bl = radiation_mode(mode) isa BeerLambert
    K_lifted =
        color_pattern(lifted_pattern(N_cells, NC; chi, bl, species_stencil = 3)[1]).colors
    return (;
        NC,
        N,
        t_res,
        direct = structured[1],
        approx = structured[2],
        t_dual,
        K_lifted,
    )
end

function _new_core(deck, n, energy)
    problem, kw = problem_for(deck, n)
    cp = setup(problem; energy_state = energy, use_ale = kw.use_ale, kw.radiation_model)
    u = _hot_state(cp)
    du = similar(u)
    cj = cp.jacobian
    J = jac_prototype(cj)
    jacobian!(J, cj, u, T_EVAL)
    all(isfinite, nonzeros(J)) || error("non-finite new-core Jacobian on $deck at $n cells")
    t_res = time_op(() -> new_residual!(du, u, cp.def, T_EVAL, cp.ws); repeats = REPEATS)
    t_jac = time_op(() -> jacobian!(J, cj, u, T_EVAL); repeats = REPEATS)
    D = only(keys(cp.ws.dual))
    C = ForwardDiff.npartials(D)
    u_d = _seeded(D, u)
    du_d = similar(u_d)
    t_dual =
        time_op(() -> new_residual!(du_d, u_d, cp.def, T_EVAL, cp.ws); repeats = REPEATS)
    bytes = @allocated jacobian!(J, cj, u, T_EVAL)
    return (;
        N = length(u),
        K = cj.n_colors,
        C,
        passes = cld(cj.n_colors, C),
        t_res,
        t_jac,
        t_dual,
        bytes,
    )
end

_us(t) = round(t * 1e6; digits = 1)
_r(x; d = 2) = round(x; digits = d)

include(joinpath(@__DIR__, "stamp.jl"))
load_start = Sys.loadavg()
started = time()
rows = []
for deck in DECKS, n in SIZES
    old = _old_core(deck, n)
    for energy in ENERGY_STATES
        new = _new_core(deck, n, energy)
        per_color_new = (new.t_dual.median - new.t_res.median) / (new.C * new.t_res.median)
        per_color_old =
            (old.t_dual.median - old.t_res.median) / (OLD_CHUNK * old.t_res.median)
        projected = cld(old.K_lifted, OLD_CHUNK) * old.t_dual.median
        row = (;
            deck,
            n,
            NC = old.NC,
            energy,
            state_old = old.N,
            state_new = new.N,
            colors = new.K,
            colors_lifted = old.K_lifted,
            batch = new.C,
            passes = new.passes,
            res_old_us = _us(old.t_res.median),
            res_new_us = _us(new.t_res.median),
            direct_us = _us(old.direct.median),
            approx_us = _us(old.approx.median),
            colored_us = _us(new.t_jac.median),
            colored_over_direct = _r(new.t_jac.median / old.direct.median),
            colored_over_approx = _r(new.t_jac.median / old.approx.median),
            projected_us = _us(projected),
            projected_over_direct = _r(projected / old.direct.median),
            projected_over_approx = _r(projected / old.approx.median),
            per_color_new = _r(per_color_new),
            per_color_old = _r(per_color_old),
            jacobian_per_color = _r(new.t_jac.median / (new.K * new.t_res.median)),
            colored_bytes = new.bytes,
        )
        push!(rows, row)
        println(join(("$k=$v" for (k, v) in pairs(row)), " "))
        flush(stdout)
    end
end
load_end = Sys.loadavg()

open(joinpath(OUT, "leftover1_jacobian_cost.csv"), "w") do io
    csv_row(io, keys(first(rows))...)
    for row in rows
        csv_row(io, values(row)...)
    end
end

open(joinpath(OUT, "leftover1_jacobian_cost.md"), "w") do io
    println(io, "# Leftover 1: one colored Jacobian against one structured assembly\n")
    println(
        io,
        "Produced by `test/newcore/gates/leftover1_jacobian_cost.jl` at $(COMMIT_STAMP). Julia $(VERSION), ",
        "$(Threads.nthreads()) thread, $(Sys.CPU_THREADS) logical CPUs ($(Sys.cpu_info()[1].model)); ",
        "load average $(round.(load_start; digits = 1)) at the start and ",
        "$(round.(load_end; digits = 1)) at the end. Warm medians of $(REPEATS) samples ",
        "(`time_op` of the Phase 0 harness), every timing in one session.\n",
    )
    println(
        io,
        "| deck | cells | energy | colors (lifted model) | batch × passes | residual old / new [µs] | ",
        "DirectSolve [µs] | ApproxSolve [µs] | colored [µs] | colored / Direct | colored / Approx | ",
        "projected / Direct | projected / Approx | per color new (old) | Jacobian per color |",
    )
    println(io, "|" * repeat("---|", 15))
    for r in rows
        println(
            io,
            "| $(r.deck) | $(r.n) | $(r.energy) | $(r.colors) ($(r.colors_lifted)) | ",
            "$(r.batch) × $(r.passes) | $(r.res_old_us) / $(r.res_new_us) | $(r.direct_us) | ",
            "$(r.approx_us) | $(r.colored_us) | $(r.colored_over_direct) | ",
            "$(r.colored_over_approx) | $(r.projected_over_direct) | $(r.projected_over_approx) | ",
            "$(r.per_color_new) ($(r.per_color_old)) | $(r.jacobian_per_color) |",
        )
    end
end
println(
    "wrote $(OUT); load $(load_start) -> $(load_end); $(round(Int, time() - started)) s",
)
