# Gate item 7 of the Phase 1 plan (§7.3): the energy balance of the
# temperature form against the old core's ledger closure on the same decks at
# the same tolerance, its convergence under refinement, and the inventory
# form's exactness.
#
# The decks are V-CPL-06's char-layer column (80 cells, 240 s, the window
# from 120 s on), Eulerian and on a moving mesh: a reactive band gasifies and
# the gas picks up about 330 K on its way out through the top, so a gas
# enthalpy defect shows.
#
# 1. Closure at the case's tolerance (abstol 1e-8, reltol 1e-6). The old core's
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
#    system.
#
# Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/energy_ledger.jl
# Writes test/newcore/gates/measurements/energy_ledger.md.

ENV["PYROLYSIS_CORE"] = "old"
ENV["PYROLYSIS_ENERGY_STATE"] = "temperature"
include(joinpath(@__DIR__, "..", "..", "verification", "casedefs.jl"))

using Dates
using ForwardDiff
using LinearAlgebra: I
using OrdinaryDiffEq: OrdinaryDiffEq, ODEFunction, ODEProblem, DiscreteCallback, Rodas5P
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

# The new core's closure over the window: ΔH against the trapezoidal ledger,
# over the integrated boundary heat.
function new_core_closure(problem; use_ale, energy_state)
    cp = Integration.setup(problem; use_ale, energy_state, TOL...)
    def, ws = cp.def, cp.ws
    du = zeros(length(cp.u0))
    acc, prev = zeros(2), zeros(2)
    E1, E2, last_t, started = Ref(NaN), Ref(NaN), Ref(0.0), Ref(false)
    function ledger!(integ)
        t = integ.t
        t >= T_WINDOW || return nothing
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
        closure = abs(E2[] - E1[] - acc[1]) / abs(acc[2]),
        boundary = acc[2],
        steps = sol.statistics.accepted,
        wall,
    )
end

# The old core's closure: V-CPL-06b's budget on the same column.
function old_core_closure(problem; use_ale)
    wall = @elapsed r = _vcpl06_ledger_solve(problem; use_ale)
    a = r.acc
    pred = a.q_bnd + a.Q_rxn + a.S_conv + a.handoffE + r.c_gas * a.xidotT + a.ale_T + a.dil
    return (; retcode = r.retcode, closure = abs(r.ΔE - pred) / abs(a.q_bnd), wall)
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

load_start = Sys.loadavg()
started = now()
case = only(filter(c -> c.id == "V-CPL-06b", load_verification_cases()))
closures = []
for use_ale in (false, true)
    deck = use_ale ? "moving" : "Eulerian"
    old = Base.invokelatest(old_core_closure, case.build(N_CELLS); use_ale)
    push!(closures, (; deck, run = "old core", old...))
    println(closures[end])
    for energy_state in (:temperature, :enthalpy)
        r = new_core_closure(case.build(N_CELLS); use_ale, energy_state)
        push!(
            closures,
            (; deck, run = "new core, $(energy_state) form", r.retcode, r.closure, r.wall),
        )
        println(closures[end], " steps=", r.steps)
        flush(stdout)
    end
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
load_end = Sys.loadavg()

_g(x; d = 3) = x isa AbstractFloat ? string(round(x; sigdigits = d)) : string(x)
open(joinpath(OUT, "energy_ledger.md"), "w") do io
    println(io, "# Gate item 7: the energy balance against the ledger\n")
    println(
        io,
        "Produced by `test/newcore/gates/energy_ledger.jl`. Julia $(VERSION), ",
        "$(Threads.nthreads()) thread, $(Sys.CPU_THREADS) logical CPUs; load average ",
        "$(round.(load_start; digits = 1)) at the start and $(round.(load_end; digits = 1)) ",
        "at the end. V-CPL-06's char-layer column.\n",
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
end
println(
    "wrote $(OUT); load $(load_start) -> $(load_end); $(round(now() - started, Second))",
)
