# ============================================================================
# The integrator's domain test: admissibility and constraint accuracy
# ============================================================================

"""
	inadmissible(u, layout, setup, abstol) -> Bool

Whether the state `u` lies outside the physical domain: a live component mass
more negative than its absolute tolerance, a live cell at or below absolute
zero ([`below_absolute_zero`](@ref)), or on a moving mesh a live cell whose
nodes have crossed (a non-positive volume) or a non-positive cross-section
area. Read from the raw state, never through a guard: the trial temperatures
keep a state outside the domain evaluable, and this test keeps it from being
accepted.

Each live cell is held to its own tolerances, a merge's receiver included:
it holds the merged cells' masses against its own `abstol`, which is what the
domain test holds every later proposal to ([`DomainCheck`](@ref)), so a state
after an event that this test accepts is one the integrator can continue
from.
"""
function inadmissible(
    u::AbstractVector,
    layout::StateLayout{NC},
    setup,
    abstol::Vector{Float64},
) where {NC}
    @inbounds for i = 1:n_cells(layout)
        layout.mask[energy_index(layout, i)] || continue
        for j = 1:NC
            k = mass_index(layout, i, j)
            u[k] < -abstol[k] && return true
        end
        if is_moving(layout)
            u[position_index(layout, i+1)] <= u[position_index(layout, i)] && return true
        end
    end
    has_chi(layout) && u[area_index(layout)] <= 0 && return true
    return below_absolute_zero(u, layout, setup) != 0
end

"""
	below_absolute_zero(u, layout, setup) -> Int

The first live cell of the state `u` whose temperature is at or below
absolute zero, or zero when there is none: its energy state at or below that
of zero kelvin at its masses, the relation the caloric inversion inverts
([`energy_state_value`](@ref)), read from the raw state without inverting it.
"""
function below_absolute_zero(u::AbstractVector, layout::StateLayout{NC}, setup) where {NC}
    energy = energy_state(layout)
    @inbounds for i = 1:n_cells(layout)
        k = energy_index(layout, i)
        layout.mask[k] || continue
        m = ntuple(j -> u[mass_index(layout, i, j)], Val(NC))
        u[k] <= energy_state_value(energy, setup, 0.0, m, i) && return i
    end
    return 0
end

"""
	ConstraintMonitor

What the domain test saw: the largest scaled constraint residual of the last
proposed state it tested (`last`) and of every accepted one (`largest`), the
steps it rejected on admissibility or on a constraint (`rejections`), the
current run of consecutive constraint rejections, the live cell at or below
absolute zero in the last proposal if it was rejected for one (`cold`, zero
otherwise) and its time, and, once the solve cannot go on, the diagnostic the
step observer ends it with.
"""
mutable struct ConstraintMonitor
    last::Float64
    largest::Float64
    rejections::Int
    consecutive::Int
    cold::Int
    cold_t::Float64
    failed::Bool
    message::String
end

ConstraintMonitor() = ConstraintMonitor(0.0, 0.0, 0, 0, 0, NaN, false, "")

"""
	MAX_CONSTRAINT_REJECTIONS

Consecutive steps rejected on a constraint residual after which the solve
ends with `ReturnCode.Failure`.
"""
const MAX_CONSTRAINT_REJECTIONS = 50

# Whether the integrator's error test accepts its current proposal.
const _ErrorTestWrapper = FunctionWrapper{Bool,Tuple{}}

"""
	DomainCheck(cp::CoreProblem; max_consecutive = MAX_CONSTRAINT_REJECTIONS)

The integrator's `isoutofdomain` test. A proposed state is rejected, and the
step shrunk by the controller's minimum factor, when

- it is [`inadmissible`](@ref), a live cell at or below absolute zero
  included (recorded in the monitor, so a solve that cannot pass it says
  so);
- under Darcy–Fick transport a live cell's porosity
  `1 − Σⱼ mⱼ/(ρⱼ(T) V)`, read from the raw masses and volume, leaves `[0, 1]`
  by more than the skeletal volume of its mass tolerances;
- an algebraic row's residual exceeds `rtol · S + a`, with `S` its own
  equation's scale and `a` the absolute tolerance of its unknown in the
  residual's units ([`constraint_residual`](@ref)).

The test reads only the proposals the integrator's error test accepts, once
bound to the integrator ([`bind_error_test!`](@ref)): the integrator evaluates
it after the error estimate, and a step it rejects is cut by the controller's
smallest factor, so a proposal the error test rejects is left to that test
and its controller's own cut. The residual is evaluated once per tested
state. After more than `max_consecutive` consecutive constraint rejections
the state is accepted and the monitor marked failed, and the step observer
ends the solve.
"""
struct DomainCheck{D<:CoreDefinition,W<:CoreWorkspace}
    def::D
    ws::W
    du::Vector{Float64}
    abstol::Vector{Float64}
    reltol::Float64
    intensity_scale::Float64
    balance_cell::Vector{Int}
    max_consecutive::Int
    monitor::ConstraintMonitor
    error_test::Base.RefValue{_ErrorTestWrapper}
end

function DomainCheck(cp; max_consecutive::Int = MAX_CONSTRAINT_REJECTIONS)
    def = cp.def
    layout = def.layout
    balance = zeros(Int, n_nodes(layout))
    cl = def.closure
    if cl isa BalanceRowClosure && !(cl.correction isa PrescribedCorrection)
        for k in cl.order
            balance[cl.row_node[k]] = k
        end
    end
    I_scale =
        has_intensities(layout) ? cp.abstol[intensity_index(layout, 1)] / cp.reltol : 1.0
    return DomainCheck(
        def,
        cp.ws,
        zeros(state_length(layout)),
        cp.abstol,
        cp.reltol,
        I_scale,
        balance,
        max_consecutive,
        ConstraintMonitor(),
        Ref(_ErrorTestWrapper(() -> true)),
    )
end

"""
	bind_error_test!(check::DomainCheck, integrator) -> check

Read the error test of `integrator` (`accept_step_controller`) before
testing a proposal; unbound, every proposal is tested.
"""
function bind_error_test!(check::DomainCheck, integrator)
    check.error_test[] = _ErrorTestWrapper(
        () -> OrdinaryDiffEqCore.accept_step_controller(integrator, integrator.alg),
    )
    return check
end

function (check::DomainCheck)(u, p, t)
    check.error_test[]() || return false
    layout, setup = check.def.layout, check.def.setup
    m = check.monitor
    if inadmissible(u, layout, setup, check.abstol)
        m.cold = below_absolute_zero(u, layout, setup)
        m.cold_t = t
        return _reject!(m)
    end
    m.cold = 0
    residual!(check.du, u, check.def, t, check.ws)
    ds, c = check.ws.derived, check.ws.contributions
    _porosity_out(check, ds, setup.gas) && return _reject!(m)
    r = constraint_residual(check, ds, c, t)
    isfinite(r) || return _reject!(m)
    m.last = r
    if r > 1
        m.consecutive += 1
        if m.consecutive > check.max_consecutive
            m.failed = true
            m.message =
                "a constraint residual stayed above its tolerance for $(m.consecutive) " *
                "consecutive steps at t = $t (scaled residual $r)"
            return false
        end
        return _reject!(m)
    end
    m.consecutive = 0
    return false
end

@inline function _reject!(m::ConstraintMonitor)
    m.rejections += 1
    return true
end

# Porosity outside [0, 1] beyond the skeletal volume of the mass tolerances;
# only Darcy–Fick transport reads a porosity.
@inline _porosity_out(check, ds, ::GasTransport) = false

function _porosity_out(
    check::DomainCheck,
    ds::DerivedState{NC},
    ::DarcyFickTransport,
) where {NC}
    layout = check.def.layout
    material = check.def.setup.material
    @inbounds for i = 1:n_cells(layout)
        layout.mask[energy_index(layout, i)] || continue
        V = ds.V[i]
        T = ds.T_trial[i]
        φ = 1 - skeletal_volume_fraction(ntuple(j -> ds.m[j, i] / V, Val(NC)), T, material)
        tol = skeletal_volume_fraction(
            ntuple(j -> check.abstol[mass_index(layout, i, j)] / V, Val(NC)),
            T,
            material,
        )
        (φ < -tol || φ > 1 + tol) && return true
    end
    return false
end

"""
	constraint_residual(check, ds, c, t) -> Float64

The largest scaled residual `|r| / (rtol · S + a)` over the live algebraic
rows, at the derived state `ds` and contributions `c` of the residual just
evaluated into `check.du`:

- a balance row of cell `k`: `S = |θₖ| Vₖ + |θ_A| Vₖ + Σ_f A_f |w_f| + Σ_f A_f |c_f|`,
  the last the sweep of the node corrections, `a` the velocity tolerance of
  its node times its face area;
- the pinned node's row and a prescribed-velocity row: `S = |wₙ|`, `a` the
  node's velocity tolerance;
- a face-intensity row: `S` the incident-intensity scale, `a` its tolerance;
- the area row: `S = |A|`; the lateral-rate row: `S = |θ_A|` plus the
  volume-weighted mean `|θ|`; `a` their tolerances.

A value above one fails the check.
"""
function constraint_residual(check::DomainCheck, ds, c, t)
    layout = check.def.layout
    du, tol, rtol = check.du, check.abstol, check.reltol
    mask = layout.mask
    worst = 0.0
    if is_moving(layout)
        g = ds.geometry
        topo = check.def.setup.topology
        @inbounds for nd = 1:n_nodes(layout)
            row = velocity_index(layout, nd)
            mask[row] || continue
            k = check.balance_cell[nd]
            if k == 0
                S, a = abs(ds.w[nd]), tol[row]
            else
                S = (abs(c.θ[k]) + abs(ds.θ_A[])) * ds.V_eval[k]
                for f in topo.cells[k].faces
                    S += g.face_area[f] * abs(ds.w[topo.faces[f].nodes[1]])
                end
                for f in topo.cells[k].faces
                    S += g.face_area[f] * abs(ds.c[topo.faces[f].nodes[1]])
                end
                # The face at node `nd` (face `f` sits at node `f` in 1D).
                a = tol[row] * g.face_area[nd]
            end
            worst = max(worst, abs(du[row]) / (rtol * S + a))
        end
    end
    if has_intensities(layout)
        @inbounds for f = 1:n_nodes(layout)
            row = intensity_index(layout, f)
            mask[row] || continue
            worst = max(worst, abs(du[row]) / (rtol * check.intensity_scale + tol[row]))
        end
    end
    if has_chi(layout)
        A = area_index(layout)
        worst = max(worst, abs(du[A]) / (rtol * abs(ds.A[]) + tol[A]))
        L = lateral_rate_index(layout)
        θ̄ = sum(i -> abs(c.θ[i]) * ds.V_eval[i], 1:n_cells(layout)) / sum(ds.V_eval)
        worst = max(worst, abs(du[L]) / (rtol * (abs(ds.θ_A[]) + θ̄) + tol[L]))
    end
    return worst
end
