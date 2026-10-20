# Gate item 7 of the Phase 1 plan (§7.3) and item 3 of the Phase 2 plan: the
# energy balance of the temperature form against the old core's ledger closure
# on the same decks at the same tolerance, its convergence under refinement,
# and the inventory form's exactness.
#
# The decks are V-CPL-06's char-layer column (80 cells, 240 s, the window
# from 120 s on), Eulerian and on a moving mesh: a reactive band gasifies and
# the gas picks up about 330 K on its way out through the top, so a gas
# enthalpy defect shows; and V-BO-04's burning column (40 cells, moving mesh)
# over the window before the first deactivation event on either core (0 to
# 290 s; the old core parks its first cells at 299.6 s, the new core
# deactivates the column at 465 s), solved without deactivation, and across
# that event (3). V-BO-04 is adiabatic, so its closures are over the change of the total enthalpy over
# the window, each core's own.
#
# 1. Closure at the case's tolerance (abstol 1e-8, reltol 1e-6 on V-CPL-06,
#    abstol 1e-9, reltol 1e-7 on V-BO-04). The old core's
#    is V-CPL-06b's budget, ΔE against the trapezoidal integral of its terms
#    over the accepted steps of the window, over the integrated boundary heat
#    (the Eulerian column has no dilation terms). The new core's is the same
#    construction on its own ledger: the total sensible enthalpy (Σ m h(T) in
#    the temperature form, the enthalpy inventories in the inventory form)
#    against the trapezoidal integral of what enters through the boundary
#    faces, the reaction heat and the formation terms, over the integrated
#    boundary heat (the boundary energy flux less the enthalpy its mass flux
#    carries).
# 2. Refinement. The new core at fixed steps with the ledger integral carried
#    as an extra state, so the integrator's own quadrature accounts for it and
#    Σ H − E is a linear invariant: its largest defect over the solve, at
#    three steps in the temperature form (the observed orders) and at the
#    coarsest in the inventory form (round-off). 20 cells, the default
#    integrator, a dense automatic-differentiation Jacobian of the augmented
#    system. V-BO-04 runs at 0.1, 0.05 and 0.025 s steps (at 0.2 s and above
#    the gas venting at the start of the burn drives a gas mass negative,
#    which a fixed step cannot reject); its column stays at its initial
#    temperature, so its enthalpy is linear in the masses and the defect is
#    round-off in both forms.
# 3. V-BO-04 across its event: the case's run with deactivation, 40 cells, at
#    its tolerance, with event location off and on. The inventory is the live
#    cells' extended energy (the inventories in the inventory form), the
#    boundary terms are taken at the live boundary faces and the sources over
#    the live cells; each invocation's lost enthalpy and regularization leave
#    the inventory at the event, its lost masses the mass inventory. A probe
#    located where the smallest threshold margin is `PRE_EVENT` gives the
#    ledger its last node before each event; the steps after it to the event
#    are what location shortens. The old core's budget does not park cells, so
#    its closure before the event is the comparison.
#
# Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/energy_ledger.jl
# Writes test/newcore/gates/measurements/energy_ledger.md.

ENV["PYROLYSIS_CORE"] = "old"
ENV["PYROLYSIS_ENERGY_STATE"] = "temperature"
include(joinpath(@__DIR__, "..", "..", "verification", "casedefs.jl"))

using ForwardDiff
using LinearAlgebra: I
using OrdinaryDiffEq: OrdinaryDiffEq, ODEFunction, ODEProblem, DiscreteCallback, Rodas5P
using OrdinaryDiffEq: ContinuousCallback, CallbackSet
using Pyrolysis.Deactivation: threshold_margin
# The new core's names are qualified: V-CPL-06's ledger, included into this
# module, calls the old core's `residual!`.
const Meshes = Pyrolysis.Meshes
const States = Pyrolysis.States
const Assembly = Pyrolysis.Assembly
const Integration = Pyrolysis.Integration

const OUT = joinpath(@__DIR__, "measurements")
mkpath(OUT)
const N_CELLS = 80
const T_WINDOW = 120.0
const TOL = (abstol = 1e-8, reltol = 1e-6)
const REFINEMENT_CELLS = 20
const STEPS = (2.0, 1.0, 0.5)
const VBO04_CELLS = 40
const VBO04_WINDOW = 290.0
const VBO04_TOL = (abstol = 1e-9, reltol = 1e-7)
const VBO04_STEPS = (0.1, 0.05, 0.025)
const ROUND_OFF = 1e-12
const PRE_EVENT = 1e-8

# The total-enthalpy rate the inventories imply: what enters through the
# boundary faces (the condition's heat flux and the enthalpy the boundary mass
# fluxes carry), the reaction heat and the formation terms.
function ledger_rate(def, c, ds)
    topo = def.setup.topology
    g = ds.geometry
    L = zero(eltype(c.F_energy))
    for f = 1:Meshes.n_faces(topo)
        topo.faces[f].cells[2] == 0 && (L -= c.F_energy[f] * g.face_area[f])
    end
    for i = 1:Meshes.n_cells(topo), j in axes(c.S_mass, 1)
        L += c.S_mass[j, i] * ds.h[j, i]
    end
    return L + sum(c.S_energy)
end

# The heat entering through the boundary faces, without the enthalpy their
# mass fluxes carry [W].
function boundary_heat(def, c, ds)
    topo = def.setup.topology
    g = ds.geometry
    q = 0.0
    for f = 1:Meshes.n_faces(topo)
        topo.faces[f].cells[2] == 0 || continue
        h = States.component_enthalpies(def.setup.laws, c.T_face[f])
        carried = sum(c.F_mass[j, f] * h[j] for j in axes(c.F_mass, 1))
        q -= (c.F_energy[f] - carried) * g.face_area[f]
    end
    return q
end

# The total sensible enthalpy of the state `u` [J].
total_enthalpy(def, ws, u) = total_enthalpy(States.energy_state(def.layout), def, ws, u)
total_enthalpy(::States.EnthalpyState, def, ws, u) =
    sum(u[States.energy_index(def.layout, i)] for i = 1:Meshes.n_cells(def.layout))
function total_enthalpy(::States.TemperatureState, def, ws, u)
    ds = ws.derived
    States.derive!(ds, u, def.layout, def.setup)
    return sum(ds.m[j, i] * ds.h[j, i] for j in axes(ds.m, 1), i in axes(ds.m, 2))
end

# The new core's terms over the window from `t_win`: ΔH, its defect against
# the trapezoidal ledger, and the integrated boundary heat.
function new_core_closure(problem; use_ale, energy_state, t_win = T_WINDOW, tol = TOL)
    cp = Integration.setup(problem; use_ale, energy_state, tol...)
    def, ws = cp.def, cp.ws
    du = zeros(length(cp.u0))
    acc, prev = zeros(2), zeros(2)
    E1, E2, last_t, started = Ref(NaN), Ref(NaN), Ref(0.0), Ref(false)
    function ledger!(integ)
        t = integ.t
        t >= t_win || return nothing
        Assembly.residual!(du, integ.u, def, t, ws)
        ds, c = Assembly.buffers(ws, def, Float64)
        cur = (ledger_rate(def, c, ds), boundary_heat(def, c, ds))
        E = total_enthalpy(def, ws, integ.u)
        if started[]
            acc .+= 0.5 .* (prev .+ cur) .* (t - last_t[])
        else
            E1[], started[] = E, true
        end
        E2[], last_t[] = E, t
        prev .= cur
        OrdinaryDiffEq.SciMLBase.derivative_discontinuity!(integ, false)
        return nothing
    end
    callback =
        DiscreteCallback((u, t, integ) -> true, ledger!; save_positions = (false, false))
    wall = @elapsed sol = Integration.integrate(cp; callback)
    return (;
        retcode = sol.retcode,
        defect = E2[] - E1[] - acc[1],
        ΔE = E2[] - E1[],
        boundary = acc[2],
        steps = sol.statistics.accepted,
        wall,
    )
end

# The inventory rates of the live region: the energy rate of `ledger_rate`
# and the mass rate, with the boundary terms at the live boundary faces and
# the sources over the live cells.
function live_rates(def, c, ds)
    topo = def.setup.topology
    g = ds.geometry
    L, M = 0.0, 0.0
    for f0 = 1:Meshes.n_faces(topo)
        topo.faces[f0].cells[2] == 0 || continue
        f = Meshes.live_face(topo, ds.active, f0)
        f == 0 && continue
        L -= c.F_energy[f] * g.face_area[f]
        M -= sum(c.F_mass[j, f] for j in axes(c.F_mass, 1)) * g.face_area[f]
    end
    for i = 1:Meshes.n_cells(topo)
        ds.active[i] || continue
        L += c.S_energy[i]
        for j in axes(c.S_mass, 1)
            L += c.S_mass[j, i] * ds.h[j, i]
        end
    end
    return L, M
end

# The live cells' extended energy [J] and mass [kg] in the state `u`, from
# its derived state `ds`.
function live_inventories(def, cells, ds, u)
    layout, laws = def.layout, def.setup.laws
    E, m = 0.0, 0.0
    for i = 1:Meshes.n_cells(layout)
        ds.active[i] || continue
        masses = ntuple(j -> u[States.mass_index(layout, i, j)], size(ds.m, 1))
        E +=
            States.energy_state(layout) isa States.EnthalpyState ?
            u[States.energy_index(layout, i)] :
            States.extended_energy(cells.caloric, laws, ds.T[i], masses, i)
        m += sum(masses)
    end
    return E, m
end

# The new core's closures across the events of a deactivating solve: the
# energy and mass defects against the trapezoidal ledger with the events'
# terms, the inventory changes, and the events' terms.
function new_core_event_closure(problem; energy_state, located, tol = VBO04_TOL)
    cp = Integration.setup(
        problem;
        use_ale = true,
        energy_state,
        deactivation = true,
        tol...,
    )
    def, ws, cells = cp.def, cp.ws, cp.deactivation
    du = zeros(length(cp.u0))
    acc, prev = zeros(2), zeros(2)
    first_node, last_node, last_t, started = zeros(2), zeros(2), Ref(0.0), Ref(false)
    function node!(integ)
        t = integ.t
        Assembly.residual!(du, integ.u, def, t, ws)
        ds, c = Assembly.buffers(ws, def, Float64)
        cur = live_rates(def, c, ds)
        inventories = live_inventories(def, cells, ds, integ.u)
        if started[]
            acc .+= 0.5 .* (prev .+ cur) .* (t - last_t[])
        else
            first_node .= inventories
            started[] = true
        end
        last_node .= inventories
        last_t[] = t
        prev .= cur
        OrdinaryDiffEq.SciMLBase.derivative_discontinuity!(integ, false)
        return nothing
    end
    margin(u) = threshold_margin(cells, def.layout, u)
    pre_event = ContinuousCallback(
        (u, t, integ) -> margin(u) - PRE_EVENT,
        nothing;
        affect_neg! = node!,
        save_positions = (false, false),
    )
    step = DiscreteCallback((u, t, integ) -> true, node!; save_positions = (false, false))
    sol = Integration.integrate(
        cp;
        callback = CallbackSet(pre_event, step),
        located_events = located,
        min_thickness = 1e-6,
    )
    events = sol.events
    lost_energy = sum(e.lost_enthalpy + e.regularization for e in events; init = 0.0)
    lost_mass = sum(sum(e.lost_mass) for e in events; init = 0.0)
    ΔE, Δm = last_node .- first_node
    return (;
        retcode = sol.retcode,
        outcome = sol.outcome,
        invocations = sol.statistics.invocations,
        events = length(events),
        t_event = isempty(events) ? NaN : first(events).t,
        energy_defect = ΔE - acc[1] + lost_energy,
        mass_defect = Δm - acc[2] + lost_mass,
        ΔE,
        Δm,
        lost_energy,
        regularization = sum(e.regularization for e in events; init = 0.0),
        created_volume = sum(e.created_volume for e in events; init = 0.0),
        steps = sol.statistics.accepted,
    )
end

# The old core's terms: V-CPL-06b's budget on the same deck (`ledger` its
# window and the deck's gas index and heat capacities).
function old_core_closure(problem; use_ale, ledger...)
    wall = @elapsed r = _vcpl06_ledger_solve(problem; use_ale, ledger...)
    a = r.acc
    pred = a.q_bnd + a.Q_rxn + a.S_conv + a.handoffE + r.c_gas * a.xidotT + a.ale_T + a.dil
    return (;
        retcode = r.retcode,
        defect = r.ΔE - pred,
        ΔE = r.ΔE,
        boundary = a.q_bnd,
        wall,
    )
end

# The largest defect of Σ H − E over a fixed-step solve with the ledger
# integral E carried as an extra state.
function ledger_defect(cp, dt, t_end)
    def, ws = cp.def, cp.ws
    n = length(cp.u0)
    function rhs!(dv, v, p, t)
        u = view(v, 1:n)
        Assembly.residual!(view(dv, 1:n), u, def, t, ws)
        ds, c = Assembly.buffers(ws, def, eltype(v))
        dv[n+1] = ledger_rate(def, c, ds)
        return nothing
    end
    function jac!(J, v, p, t)
        ForwardDiff.jacobian!(J, (dv, x) -> rhs!(dv, x, p, t), zeros(n + 1), v)
        return nothing
    end
    M =
        all(States.is_differential(def.layout)) ? I :
        cat(Matrix(States.mass_matrix(def.layout)), 1.0; dims = (1, 2))
    f = ODEFunction(rhs!; jac = jac!, mass_matrix = M)
    sol = OrdinaryDiffEq.solve(
        ODEProblem(f, vcat(cp.u0, 0.0), (0.0, t_end)),
        Rodas5P(autodiff = OrdinaryDiffEq.AutoFiniteDiff());
        adaptive = false,
        dt,
    )
    H0 = total_enthalpy(def, ws, cp.u0)
    defect = maximum(abs(total_enthalpy(def, ws, v[1:n]) - H0 - v[n+1]) for v in sol.u)
    return sol.retcode, defect, maximum(abs(v[n+1]) for v in sol.u)
end

include(joinpath(@__DIR__, "stamp.jl"))
load_start = Sys.loadavg()
started = time()
case = only(filter(c -> c.id == "V-CPL-06b", load_verification_cases()))
closures = []
for use_ale in (false, true)
    deck = use_ale ? "moving" : "Eulerian"
    old = Base.invokelatest(old_core_closure, case.build(N_CELLS); use_ale)
    push!(
        closures,
        (;
            deck,
            run = "old core",
            old.retcode,
            closure = abs(old.defect) / abs(old.boundary),
            old.wall,
        ),
    )
    println(closures[end])
    for energy_state in (:temperature, :enthalpy)
        r = new_core_closure(case.build(N_CELLS); use_ale, energy_state)
        push!(
            closures,
            (;
                deck,
                run = "new core, $(energy_state) form",
                r.retcode,
                closure = abs(r.defect) / abs(r.boundary),
                r.wall,
            ),
        )
        println(closures[end], " steps=", r.steps)
        flush(stdout)
    end
end

# V-BO-04 over the window before its first event, each core's closure over
# its own enthalpy change (the old core's E = Σ ρc_p T V, the new core's the
# sensible enthalpy of its laws).
vbo04_problem(n) = Base.invokelatest(_vbo04_setup, n; t_end = VBO04_WINDOW).problem
vbo04 = []
old = Base.invokelatest(
    old_core_closure,
    vbo04_problem(VBO04_CELLS);
    use_ale = true,
    t_win = 0.0,
    gas = 2,
    c_gas = 1000.0,
    c_solid = 1000.0,
)
push!(vbo04, (; run = "old core", old.retcode, closure = abs(old.defect) / abs(old.ΔE)))
println(vbo04[end])
for energy_state in (:temperature, :enthalpy)
    r = new_core_closure(
        vbo04_problem(VBO04_CELLS);
        use_ale = true,
        energy_state,
        t_win = 0.0,
        tol = VBO04_TOL,
    )
    push!(
        vbo04,
        (;
            run = "new core, $(energy_state) form",
            r.retcode,
            closure = abs(r.defect) / abs(r.ΔE),
        ),
    )
    println(vbo04[end], " steps=", r.steps)
    flush(stdout)
end

# V-BO-04 across its event, with event location off and on.
vbo04_event = []
for energy_state in (:temperature, :enthalpy), located in (false, true)
    r = new_core_event_closure(
        Base.invokelatest(_vbo04_setup, VBO04_CELLS).problem;
        energy_state,
        located,
    )
    push!(
        vbo04_event,
        (;
            energy_state,
            location = located ? "on" : "off",
            r...,
            energy_closure = abs(r.energy_defect) / abs(r.ΔE),
            mass_closure = abs(r.mass_defect) / abs(r.Δm),
        ),
    )
    println(vbo04_event[end])
    flush(stdout)
end

refinement = []
for use_ale in (false, true)
    deck = use_ale ? "moving" : "Eulerian"
    problem = case.build(REFINEMENT_CELLS)
    t_end = problem.tspan[2]
    for energy_state in (:temperature, :enthalpy)
        steps = energy_state === :temperature ? STEPS : STEPS[1:1]
        for dt in steps
            cp = Integration.setup(problem; use_ale, energy_state, TOL...)
            retcode, defect, scale = ledger_defect(cp, dt, t_end)
            push!(refinement, (; deck, energy_state, dt, retcode, defect, scale))
            println(refinement[end])
            flush(stdout)
        end
    end
end
vbo04_refinement = []
for energy_state in (:temperature, :enthalpy)
    steps = energy_state === :temperature ? VBO04_STEPS : VBO04_STEPS[1:1]
    for dt in steps
        cp = Integration.setup(
            vbo04_problem(REFINEMENT_CELLS);
            use_ale = true,
            energy_state,
            VBO04_TOL...,
        )
        retcode, defect, scale = ledger_defect(cp, dt, VBO04_WINDOW)
        push!(
            vbo04_refinement,
            (; deck = "V-BO-04", energy_state, dt, retcode, defect, scale),
        )
        println(vbo04_refinement[end])
        flush(stdout)
    end
end
load_end = Sys.loadavg()

_g(x; d = 3) = x isa AbstractFloat ? string(round(x; sigdigits = d)) : string(x)
open(joinpath(OUT, "energy_ledger.md"), "w") do io
    println(io, "# The energy balance against the ledger\n")
    println(
        io,
        "Produced by `test/newcore/gates/energy_ledger.jl` at $(COMMIT_STAMP). Julia $(VERSION), ",
        "$(Threads.nthreads()) thread, $(Sys.CPU_THREADS) logical CPUs; load average ",
        "$(round.(load_start; digits = 1)) at the start and $(round.(load_end; digits = 1)) ",
        "at the end. Phase 1's gate item 7 on V-CPL-06's char-layer column, and Phase 2's item 3 ",
        "on V-BO-04's burning column before and across its event.\n",
    )
    println(
        io,
        "## Closure at the case's tolerance\n\n$(N_CELLS) cells, abstol $(TOL.abstol), ",
        "reltol $(TOL.reltol), the window t ≥ $(T_WINDOW) s; |ΔE − ∫ ledger| over the ",
        "integrated boundary heat, trapezoidal over the accepted steps on both cores.\n",
    )
    println(io, "| deck | run | retcode | closure |")
    println(io, "|---|---|---|---|")
    for r in closures
        println(io, "| $(r.deck) | $(r.run) | $(r.retcode) | $(_g(r.closure)) |")
    end
    println(
        io,
        "\n## Refinement with the ledger integrated by the integrator\n\n",
        "$(REFINEMENT_CELLS) cells, 0 to 240 s at fixed steps; the largest |Σ H − H₀ − E| ",
        "over the solve [J] against the largest |E|.\n",
    )
    println(io, "| deck | form | dt [s] | retcode | defect [J] | scale [J] | order |")
    println(io, "|---|---|---|---|---|---|---|")
    for (k, r) in enumerate(refinement)
        p = k > 1 ? refinement[k-1] : nothing
        order =
            p !== nothing && p.deck == r.deck && p.energy_state == r.energy_state ?
            _g(log2(p.defect / r.defect)) : ""
        println(
            io,
            "| $(r.deck) | $(r.energy_state) | $(r.dt) | $(r.retcode) | $(_g(r.defect)) | ",
            "$(_g(r.scale)) | $(order) |",
        )
    end
    println(
        io,
        "\n## V-BO-04 before its first event\n\n$(VBO04_CELLS) cells on a moving mesh, abstol ",
        "$(VBO04_TOL.abstol), reltol $(VBO04_TOL.reltol) (the case's), 0 to $(VBO04_WINDOW) s ",
        "without deactivation; |ΔE − ∫ ledger| over |ΔE|, each core's own, trapezoidal over the ",
        "accepted steps (the old core's budget is V-CPL-06b's).\n",
    )
    println(io, "| run | retcode | closure |")
    println(io, "|---|---|---|")
    for r in vbo04
        println(io, "| $(r.run) | $(r.retcode) | $(_g(r.closure)) |")
    end
    println(
        io,
        "\n$(REFINEMENT_CELLS) cells, 0 to $(VBO04_WINDOW) s at fixed steps; the largest ",
        "|Σ H − H₀ − E| over the solve [J] against the largest |E|, and their ratio. The column ",
        "stays at its initial temperature, so its enthalpy is linear in the masses in both forms ",
        "and the ledger an invariant the integrator keeps to round-off; an order is given only ",
        "between defects above $(ROUND_OFF) of the scale.\n",
    )
    println(io, "| form | dt [s] | retcode | defect [J] | scale [J] | ratio | order |")
    println(io, "|---|---|---|---|---|---|---|")
    for (k, r) in enumerate(vbo04_refinement)
        p = k > 1 ? vbo04_refinement[k-1] : nothing
        resolved(x) = x.defect > ROUND_OFF * x.scale
        order =
            p !== nothing &&
            p.energy_state == r.energy_state &&
            resolved(p) &&
            resolved(r) ? _g(log2(p.defect / r.defect)) : ""
        println(
            io,
            "| $(r.energy_state) | $(r.dt) | $(r.retcode) | $(_g(r.defect)) | ",
            "$(_g(r.scale)) | $(_g(r.defect / r.scale)) | $(order) |",
        )
    end
    println(
        io,
        "\n## V-BO-04 across its event\n\nThe case's run, $(VBO04_CELLS) cells on a moving ",
        "mesh with deactivation, abstol $(VBO04_TOL.abstol), reltol $(VBO04_TOL.reltol), ",
        "event location off and on. Inventories over the live cells (the extended energy, ",
        "the inventory itself in the enthalpy form, and the mass), the boundary terms at the ",
        "live boundary faces, the sources over the live cells, trapezoidal over the accepted ",
        "steps with a node where the smallest threshold margin is $(PRE_EVENT); each ",
        "invocation's lost enthalpy and regularization and its lost masses leave the ",
        "inventories at the event. Closure: |Δ − ∫ ledger + lost| over |Δ|, energy and mass. ",
        "Event: the first invocation's time [s]; lost energy, regularization: the events' ",
        "sums over |ΔE|; created volume: the events' sum [m³]. The old core's budget does ",
        "not park cells; its closure before the event is the table above.\n",
    )
    println(
        io,
        "| form | location | retcode | outcome | invocations | events | event [s] | energy closure | mass closure | lost energy | regularization | created volume | steps |",
    )
    println(io, "|" * repeat("---|", 13))
    for r in vbo04_event
        println(
            io,
            "| $(r.energy_state) | $(r.location) | $(r.retcode) | $(r.outcome) | $(r.invocations) | ",
            "$(r.events) | $(_g(r.t_event; d = 7)) | $(_g(r.energy_closure)) | $(_g(r.mass_closure)) | ",
            "$(_g(abs(r.lost_energy) / abs(r.ΔE))) | $(_g(abs(r.regularization) / abs(r.ΔE))) | ",
            "$(_g(r.created_volume)) | $(r.steps) |",
        )
    end
end
println(
    "wrote $(OUT); load $(load_start) -> $(load_end); $(round(Int, time() - started)) s",
)
