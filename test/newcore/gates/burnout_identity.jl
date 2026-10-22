# Gate item 5 of the Phase 2 plan (§7.3): the decks that must not change under
# the correction, hashed on the finite-volume core. Each deck is solved as the
# suite solves it (its builder at the CI resolution, its solve keywords) and
# the digest of `(T, m, z, t)` of its solution recorded; the Douglas fir deck
# of examples/07_douglas_fir_cone.jl over its first 60 s joins them.
#
# The decks split by the correction's installation rule. (A) A deck the rule
# leaves without a correction (a fixed mesh, no deactivation, or no layer that
# can vanish) must equal its reference; on a moving mesh it is solved again
# with an equidistribution correction given explicitly, which must leave every
# temperature and mass within 1e-10 of the first solve's, relative to the
# field's largest, on the same accepted steps, with the correction exactly
# zero. (B) A deck the rule installs the correction on records its digest at
# the package defaults as the reference, beside its steps without the
# correction and the largest correction velocity; on V-BO-04's uniform column
# that velocity must be exactly zero.
#
# Per moving deck, the thinnest live cell over its initial thickness at the
# saved states of both of its solves (item 5(A)'s smallest `h/h₀`). Per class
# B deck, from the saved positions and the correction the rule binds, the
# largest node activation and the most nodes active at one saved time (a node
# between two live cells, with nonzero gain and nonzero activation; item
# 5(B)).
#
# With `--record` the digests become the reference
# (measurements/burnout_identity_reference.csv; plan §12 records the tree it
# was taken on); without it every digest is compared with the reference and
# the table says which decks are bit-identical.
#
# Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/burnout_identity.jl [--record]
# Writes test/newcore/gates/measurements/burnout_identity.md.

ENV["PYROLYSIS_CORE"] = "new"
ENV["PYROLYSIS_ENERGY_STATE"] = "temperature"
using Pyrolysis
using Printf
include(joinpath(@__DIR__, "stamp.jl"))

const OUT = joinpath(@__DIR__, "measurements")
mkpath(OUT)
const REFERENCE = joinpath(OUT, "burnout_identity_reference.csv")
const RECORD = "--record" in ARGS
const DECKS = (
    "V-CPL-02",
    "V-CPL-02b",
    "V-BO-04",
    "V-CPL-03",
    "V-CPL-06a",
    "V-HC-05",
    "V-RAD-01",
    "V-RX-02d",
)
const DF_SECONDS = 60.0

function load_cases()
    m = Module(:cases_identity)
    Core.eval(m, :(include(p) = Base.include($m, p)))
    Base.include(m, joinpath(@__DIR__, "..", "..", "verification", "casedefs.jl"))
    return m, Base.invokelatest(m.load_verification_cases)
end

# The old core's depletion keywords have no counterpart on the finite-volume
# core; `core_solve` strips them for the routed cases, and this does the same
# for a case the suite does not route.
function _core_keywords(kw)
    drop = (:depletion_threshold, :min_cells, :jacobian, :linear_solver)
    kept = [k => v for (k, v) in pairs(kw) if !(k in drop)]
    nt = (; kept...)
    if haskey(nt, :handle_depletion)
        d = nt.handle_depletion
        nt = (;
            (k => v for (k, v) in pairs(nt) if k != :handle_depletion)...,
            deactivation = d,
        )
    end
    return nt
end

digest(sol) = string(hash((sol.T, sol.m, sol.z, sol.t)); base = 16)

# The Douglas fir deck: the material and boundary conditions of the example,
# without its printing, plotting and solve, on the native mesh of its
# thickness and cell count.
function douglas_fir_problem()
    lines = readlines(
        joinpath(@__DIR__, "..", "..", "..", "examples", "07_douglas_fir_cone.jl"),
    )
    a = findfirst(startswith("const T_HFM_MAX"), lines)
    b = findfirst(startswith("solution = solve("), lines)
    out = String[]
    skip = false
    for line in lines[a:(b-1)]
        if skip
            startswith(line, ")") && (skip = false)
            continue
        end
        if startswith(line, "println(")
            skip = !endswith(rstrip(line), ")")
            continue
        end
        startswith(line, "flush(") && continue
        push!(out, line)
    end
    m = Module(:douglas_fir_deck)
    Core.eval(m, :(using Pyrolysis))
    Core.eval(m, Meta.parseall(join(out, "\n")))
    # Built by a second `eval`, which runs in the latest world, so the deck's
    # fresh bindings are visible.
    return Core.eval(
        m,
        :(PyrolysisProblem(
            mesh = create_uniform_mesh(sample_thickness, n_cells),
            material = douglas_fir,
            bc_set = bc_set,
            T_initial = T_initial_func,
            ξ_initial = ξ_initial_funcs,
            tspan = (0.0, $DF_SECONDS),
            dt_initial = dt_initial,
            output_times = collect(0.0:1.0:($DF_SECONDS)),
        )),
    )
end

# The correction the installation rule binds on the deck solved with the case
# keywords `kw`, with the cells' initial thicknesses; `nothing` on a fixed
# mesh.
function bound(problem, kw)
    use_ale = get(kw, :use_ale, false)
    deactivation = get(kw, :handle_depletion, get(kw, :deactivation, false))
    use_ale || return nothing
    cp = Base.invokelatest(Pyrolysis.Integration.setup, problem; use_ale, deactivation)
    return (; correction = cp.def.closure.correction, h0 = cp.def.setup.reference.thickness)
end

installed(b) = b !== nothing && b.correction isa Pyrolysis.Closures.LayerCorrection

# The largest difference of two solutions' temperatures and masses, each
# relative to that field's largest value.
function difference(a, b)
    size(a.T) == size(b.T) || return Inf
    d = 0.0
    for (x, y) in ((a.T, b.T), (a.m, b.m))
        keep = [k for k in eachindex(x) if isfinite(x[k]) && isfinite(y[k])]
        d = max(d, maximum(abs, x[keep] .- y[keep]; init = 0.0) / maximum(abs, y[keep]))
    end
    return d
end

largest_c(sol) = maximum(abs, filter(isfinite, sol.correction_velocity); init = 0.0)

# The thinnest live cell over its initial thickness at the saved states (the
# positions are `NaN` beyond the outermost live node).
function thinnest(sol, h0)
    z, r = sol.z, Inf
    for j in axes(z, 2), i in eachindex(h0)
        h = z[i+1, j] - z[i, j]
        isfinite(h) && (r = min(r, h / h0[i]))
    end
    return r
end

# The largest node activation of the bound correction `corr` at the saved
# states, and the most nodes active at one saved time with that time: a node
# between two live cells, with nonzero gain and nonzero activation, the smooth
# "or" of its two cells' switches.
function activation_extent(sol, corr, h0)
    z = sol.z
    largest, most, t_most = 0.0, 0, NaN
    for j in axes(z, 2)
        active = 0
        for nd = 2:(size(z, 1)-1)
            h_in, h_out = z[nd, j] - z[nd-1, j], z[nd+1, j] - z[nd, j]
            (iszero(corr.gain[nd]) || !isfinite(h_in) || !isfinite(h_out)) && continue
            a_in = Pyrolysis.Closures.activation(corr, h_in, h0[nd-1])
            a_out = Pyrolysis.Closures.activation(corr, h_out, h0[nd])
            a = 1 - (1 - a_in) * (1 - a_out)
            largest = max(largest, a)
            a > 0 && (active += 1)
        end
        active > most && ((most, t_most) = (active, sol.t[j]))
    end
    return (; largest, most, t_most)
end

function main()
    m, cases = load_cases()
    rows = []
    record(deck, problem, sol, b, solve) = begin
        moving, rule = b !== nothing, installed(b)
        second =
            rule ? solve(Pyrolysis.Closures.NoCorrection()) :
            moving ? solve(Pyrolysis.Closures.EquidistributionCorrection()) : nothing
        r = (;
            deck,
            cells = size(sol.T, 1),
            digest = digest(sol),
            steps = sol.statistics.accepted,
            retcode = sol.retcode,
            moving,
            installed = rule,
            check = if rule
                (; steps = second.statistics.accepted, largest = largest_c(sol))
            elseif moving
                (;
                    same = second.statistics.accepted == sol.statistics.accepted &&
                           second.statistics.rejected == sol.statistics.rejected,
                    difference = difference(second, sol),
                    largest = largest_c(second),
                )
            else
                nothing
            end,
            thinnest = moving ? min(thinnest(sol, b.h0), thinnest(second, b.h0)) : NaN,
            extent = rule ? activation_extent(sol, b.correction, b.h0) : nothing,
        )
        push!(rows, r)
        println(r)
        flush(stdout)
    end
    for id in DECKS
        case = only(filter(c -> c.id == id, cases))
        problem = Base.invokelatest(case.build, case.ci.n_cells)
        solve(extra...) =
            Base.invokelatest(m.runs_on_core, case) ?
            Base.invokelatest(
                m._case_solve,
                case,
                problem;
                case.solve_kwargs...,
                extra...,
            ) :
            Base.invokelatest(
                Pyrolysis.Integration.solve,
                problem;
                _core_keywords(case.solve_kwargs)...,
                extra...,
            )
        kw = _core_keywords(case.solve_kwargs)
        record(id, problem, solve(), bound(problem, kw), c -> solve(:correction => c))
    end
    df = douglas_fir_problem()
    df_solve(extra...) = Base.invokelatest(
        Pyrolysis.Integration.solve,
        df;
        use_ale = true,
        radiation_model = SURFACE_ABSORPTION,
        abstol = 1e-9,
        reltol = 1e-7,
        dtmin = 1e-14,
        extra...,
    )
    record(
        "Douglas fir, 60 s",
        df,
        df_solve(),
        bound(df, (; use_ale = true)),
        c -> df_solve(:correction => c),
    )

    reference = Dict{String,String}()
    if isfile(REFERENCE)
        for line in readlines(REFERENCE)[2:end]
            deck, dig = rsplit(line, ","; limit = 2)
            reference[strip(deck, '"')] = dig
        end
    end
    if RECORD
        open(REFERENCE, "w") do io
            println(io, "deck,digest")
            for r in rows
                println(io, "\"", r.deck, "\",", r.digest)
            end
        end
    end
    _g(x) = string(round(x; sigdigits = 3))
    open(joinpath(OUT, "burnout_identity.md"), "w") do io
        println(io, "# Burnout identity: the decks the correction must not change\n")
        println(
            io,
            "Generated by `test/newcore/gates/burnout_identity.jl` at $(COMMIT_STAMP), finite-volume core, temperature form, each deck at its CI resolution with its solve keywords; the digest is `hash((T, m, z, t))` of the solution. The reference digests are `burnout_identity_reference.csv`$(RECORD ? ", recorded by this run" : ""). Class A: the installation rule leaves the deck without a correction; on a moving mesh the deck solved again with an equidistribution correction given explicitly must take the same steps, stay within 1e-10 (relative) and keep the correction at zero. Class B: the rule installs the correction; its steps without it and the largest correction velocity are recorded. Thinnest h/h₀: per moving deck, the thinnest live cell over its initial thickness at the saved states of both solves. Activation (class B): the largest node activation at the saved states, and the most nodes active at one saved time (between two live cells, with nonzero gain and activation) with that time [s].\n",
        )
        println(
            io,
            "| deck | cells | steps | retcode | class | digest | equal to the reference | explicit correction (A) or without the correction (B) | thinnest h/h₀ | activation: largest; most active nodes (t) |",
        )
        println(io, "|---|---|---|---|---|---|---|---|---|---|")
        for r in rows
            ref = get(reference, r.deck, "")
            eq =
                RECORD ? "recorded" :
                (isempty(ref) ? "no reference" : (ref == r.digest ? "yes" : "**no**"))
            class = r.installed ? "B" : (r.moving ? "A" : "A (fixed mesh)")
            check =
                r.check === nothing ? "" :
                r.installed ?
                "$(r.check.steps) steps without; largest |c| $(_g(r.check.largest)) m/s" :
                "same steps: $(r.check.same); difference $(_g(r.check.difference)); largest |c| $(_g(r.check.largest))"
            thin = r.moving ? _g(r.thinnest) : ""
            extent =
                r.extent === nothing ? "" :
                "$(_g(r.extent.largest)); $(r.extent.most) " *
                (r.extent.most > 0 ? "($(_g(r.extent.t_most)))" : "")
            println(
                io,
                "| $(r.deck) | $(r.cells) | $(r.steps) | $(r.retcode) | $(class) | `$(r.digest)` | $(eq) | $(check) | $(thin) | $(extent) |",
            )
        end
    end
    return rows
end

main()
