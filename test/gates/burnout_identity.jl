# The decks whose solutions the layer correction must leave unchanged, and
# the digests that identify those solutions. Each deck is solved as the
# verification suite solves it (its builder at the CI resolution, its solve
# keywords, the suite's solve) and the digest of `(T, m, z, t)` of its
# solution recorded; the Douglas fir deck of examples/07_douglas_fir_cone.jl
# over its first 60 s joins them, and so does V-CPL-06a, the char-layer
# column of `char_layer_problem` below at 80 cells, which no verification
# case builds and which is solved in the suite's energy state with no other
# keyword.
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
# saved states of both of its solves (class A's smallest `h/h₀`). Per class
# B deck, from the saved positions and the correction the rule binds, the
# largest node activation and the most nodes active at one saved time (a node
# between two live cells, with nonzero gain and nonzero activation).
#
# With `--record` the digests become the reference
# (measurements/burnout_identity_reference.csv); without it every digest is
# compared with the reference and the table says which decks are
# bit-identical. Including this file defines the decks and the functions
# without solving anything (`deletion_identity.jl` and `ttfx.jl` do).
#
# Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 test/gates/burnout_identity.jl [--record]
# Writes test/gates/measurements/burnout_identity.md.

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
const CHAR_LAYER = "V-CPL-06a"
const CHAR_LAYER_CELLS = 80
const DF_SECONDS = 60.0

# The binding `name` of the module `m` in the latest world: the cases are
# loaded into `m` at run time, in a newer world than the caller's.
latest(m, name) = Base.invokelatest(getglobal, m, name)

# The verification cases loaded into a module of their own; `casedefs.jl`
# reads the energy state from `PYROLYSIS_ENERGY_STATE` when it loads.
function load_cases()
    m = Module(:cases_identity)
    Core.eval(m, :(include(p) = Base.include($m, p)))
    Base.include(m, joinpath(@__DIR__, "..", "verification", "casedefs.jl"))
    return m, Base.invokelatest(latest(m, :load_verification_cases))
end

digest(sol) = string(hash((sol.T, sol.m, sol.z, sol.t)); base = 16)

# The char-layer column V-CPL-06a: a 10 mm column of an inert matrix and,
# at z/L ∈ [0.2, 0.4], a fuel band that gasifies with no heat of reaction at
# about 500–600 K; the gas diffuses up through the hotter layer above (the
# top held near 900 K by convection) and vents at the top, picking up about
# 330 K on its way out. Every component carries a gas-transfer coefficient,
# so the solid-dominated column transports gas in its interior.
function char_layer_problem(n_cells; t_end = 240.0)
    L = 0.01
    matrixc = SolidComponent(:matrix, ρ = 2000.0, c = 1000.0, k = 1.0, λ = 1e-4)
    fuel = SolidComponent(:fuel, ρ = 2000.0, c = 1000.0, k = 1.0, λ = 1e-4)
    gas = GaseousComponent(:gas, M = 0.05, c = 1000.0, k = 0.02, λ = 1e-4)
    rxn = Reaction(:gasify, 2 => 3; A = 9.3e3, E = 8.0e4, h = 0.0)
    mat = Material(name = :VCPL06, components = (matrixc, fuel, gas), reactions = (rxn,))
    mesh = create_uniform_mesh(L, n_cells)
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

# The deck `id` with the cases `cases` of the module `m`: its problem, its
# solve with further keywords (`solve(:correction => c)`) and its solve
# keywords. A case deck is built and solved as the suite builds and solves
# it; the char-layer column through the suite's solve, its energy state
# alone.
function deck(m, cases, id)
    if id == CHAR_LAYER
        problem = char_layer_problem(CHAR_LAYER_CELLS)
        solve = (extra...) -> Base.invokelatest(latest(m, :suite_solve), problem; extra...)
        return problem, solve, (;)
    end
    case = only(filter(c -> c.id == id, cases))
    problem = Base.invokelatest(case.build, case.ci.n_cells)
    solve =
        (extra...) -> Base.invokelatest(
            latest(m, :_case_solve),
            case,
            problem;
            case.solve_kwargs...,
            extra...,
        )
    return problem, solve, case.solve_kwargs
end

# The Douglas fir deck: the material and boundary conditions of the example,
# without its printing, plotting and solve, on the native mesh of its
# thickness and cell count.
function douglas_fir_problem()
    lines = readlines(joinpath(@__DIR__, "..", "..", "examples", "07_douglas_fir_cone.jl"))
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

# The Douglas fir deck's solve, with further keywords.
douglas_fir_solve(df, extra...) = Base.invokelatest(
    Pyrolysis.solve,
    df;
    use_ale = true,
    radiation_model = SURFACE_ABSORPTION,
    abstol = 1e-9,
    reltol = 1e-7,
    dtmin = 1e-14,
    extra...,
)

# The correction the installation rule binds on the deck solved with the case
# keywords `kw`, with the cells' initial thicknesses; `nothing` on a fixed
# mesh.
function bound(problem, kw)
    use_ale = get(kw, :use_ale, false)
    deactivation = get(kw, :deactivation, false)
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
        problem, solve, kw = deck(m, cases, id)
        record(id, problem, solve(), bound(problem, kw), c -> solve(:correction => c))
    end
    df = douglas_fir_problem()
    record(
        "Douglas fir, 60 s",
        df,
        douglas_fir_solve(df),
        bound(df, (; use_ale = true)),
        c -> douglas_fir_solve(df, :correction => c),
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
            "Generated by `test/gates/burnout_identity.jl` at $(COMMIT_STAMP), temperature form, each deck at its CI resolution with its solve keywords; the digest is `hash((T, m, z, t))` of the solution. The reference digests are `burnout_identity_reference.csv`$(RECORD ? ", recorded by this run" : ""). Class A: the installation rule leaves the deck without a correction; on a moving mesh the deck solved again with an equidistribution correction given explicitly must take the same steps, stay within 1e-10 (relative) and keep the correction at zero. Class B: the rule installs the correction; its steps without it and the largest correction velocity are recorded. Thinnest h/h₀: per moving deck, the thinnest live cell over its initial thickness at the saved states of both solves. Activation (class B): the largest node activation at the saved states, and the most nodes active at one saved time (between two live cells, with nonzero gain and activation) with that time [s].\n",
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

abspath(PROGRAM_FILE) == @__FILE__() && main()
