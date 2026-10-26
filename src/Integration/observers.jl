# ============================================================================
# Read-only observers: the progress bar and the conservation recorder
# ============================================================================

# The progress bar, advanced after every accepted step to the percentage of
# the time span covered, under the policy's flag rule.
function _progress_observer(policy, tspan)
    t0, t1 = Float64(tspan[1]), Float64(tspan[2])
    bar = Progress(100; desc = "Solving: ", showspeed = true)
    shown = Ref(0)
    return DiscreteCallback(
        (u, t, integrator) -> true,
        integrator -> begin
            p = floor(Int, 100 * (integrator.t - t0) / (t1 - t0))
            if p > shown[]
                ProgressMeter.update!(bar, p)
                shown[] = p
            end
            apply_observer_flag!(integrator, policy)
            return nothing
        end;
        save_positions = (false, false),
        finalize = (c, u, t, integrator) -> ProgressMeter.finish!(bar),
    )
end

# The three-point Gauss–Legendre nodes and weights on [−1, 1].
const _GAUSS3 = ((-sqrt(3 / 5), 5 / 9), (0.0, 8 / 9), (sqrt(3 / 5), 5 / 9))

const _SpeciesWrapper = FunctionWrapper{Vector{Float64},Tuple{_StateVector,Float64}}
const _PairWrapper = FunctionWrapper{NTuple{2,Float64},Tuple{_StateVector,Float64}}
const _ScalarWrapper = FunctionWrapper{Float64,Tuple{_StateVector}}

# A conservation recorder bound to a set-up problem: the live region's rates
# and inventories at a state, each from a residual evaluation in a workspace
# of its own so the integrator's is untouched, and the conditions of the
# events that can end a step early (a located deactivation, the thickness
# stop), each crossing zero downward.
struct _BoundRecorder
    rec::ConservationRecorder
    rates::_PairWrapper
    inventories::_PairWrapper
    crossings::Vector{_ScalarWrapper}
    species_rates::_SpeciesWrapper
    species_inventories::_SpeciesWrapper
end

function _BoundRecorder(rec::ConservationRecorder, cp::CoreProblem, crossings)
    def = cp.def
    ws = CoreWorkspace(def, cp.T0)
    ext = _extension(def)
    du = zeros(length(cp.u0))
    rates = _PairWrapper((u, t) -> begin
        residual!(du, u, def, t, ws)
        ds, c = buffers(ws, def, Float64)
        return ledger_rates(def, c, ds)
    end)
    inventories = _PairWrapper((u, t) -> begin
        residual!(du, u, def, t, ws)
        ds, _ = buffers(ws, def, Float64)
        return live_inventories(def, ext, ds, u)
    end)
    ids=gas_component_indices(def.setup.material.components)
    species_rates=_SpeciesWrapper(
        (u, t)->begin
            residual!(du, u, def, t, ws)
            ds, c=buffers(ws, def, Float64)
            v=zeros(3length(ids))
            for (g, j) in enumerate(ids)
                v[g]=sum(c.S_mass[j, i] for i in eachindex(ds.active) if ds.active[i])
                for f0 = 1:n_faces(def.setup.topology)
                    def.setup.topology.faces[f0].cells[2]==0 || continue
                    f=live_face(def.setup.topology, ds.active, f0)
                    f==0 && continue
                    N=Outputs._outward_sign(def.setup.topology, ds.active, f)*c.F_mass[
                        j,
                        f,
                    ]*ds.geometry.face_area[f]
                    v[length(ids)+g]+=max(N, 0.0)
                    v[2length(ids)+g]+=max(-N, 0.0)
                end
            end
            v
        end,
    )
    species_inventories=_SpeciesWrapper(
        (u, t)->begin
            residual!(du, u, def, t, ws)
            ds, _=buffers(ws, def, Float64)
            [sum(ds.m[j, i] for i in eachindex(ds.active) if ds.active[i]) for j in ids]
        end,
    )
    rec.species=[def.setup.material.components[j].name for j in ids]
    rec.species_production=zeros(length(ids))
    rec.species_escape=zeros(length(ids))
    rec.species_ingress=zeros(length(ids))
    rec.species_initial=species_inventories(cp.u0, Float64(cp.problem.tspan[1]))
    rec.species_final=copy(rec.species_initial)
    rec.step_species=zeros(3length(ids))
    rec.cut_species=copy(rec.step_species)
    rec.mass_outflow, rec.energy_inflow, rec.steps = 0.0, 0.0, 0
    rec.reading, rec.cut = false, NaN
    t0 = Float64(cp.problem.tspan[1])
    rec.t_done = t0
    rec.initial = rec.final = inventories(cp.u0, t0)
    return _BoundRecorder(
        rec,
        rates,
        inventories,
        collect(_ScalarWrapper, crossings),
        species_rates,
        species_inventories,
    )
end

# The integrals of the outflow and the enthalpy rate over [a, b] on the
# step's interpolant.
function _step_integrals(b::_BoundRecorder, integrator, t0, t1, species)
    fill!(species, 0.0)
    E, out = 0.0, 0.0
    for (x, w) in _GAUSS3
        τ = (t0 + t1) / 2 + (t1 - t0) / 2 * x
        u=integrator(τ)
        e, o = b.rates(u, τ)
        species .+= (w*(t1-t0)/2) .* b.species_rates(u, τ)
        E += w * e
        out += w * o
    end
    return ((t1 - t0) / 2 * E, (t1 - t0) / 2 * out)
end

# Reads the step that just ended, once, before its events: its integrals, and
# where an event condition first crosses zero inside it the integrals up to
# that time, found by bisection on the interpolant.
function _read_step!(b::_BoundRecorder, u, t, integrator)
    rec = b.rec
    (t == integrator.t > integrator.tprev) || return nothing
    rec.reading && integrator.tprev == rec.tprev && return nothing
    rec.reading, rec.t, rec.tprev = true, t, integrator.tprev
    rec.step = _step_integrals(b, integrator, integrator.tprev, t, rec.step_species)
    rec.cut = NaN
    for g in b.crossings
        (g(u) <= 0 && g(integrator.uprev) > 0) || continue
        lo, hi = integrator.tprev, t
        for _ = 1:60
            mid = (lo + hi) / 2
            g(integrator(mid)) > 0 ? (lo = mid) : (hi = mid)
        end
        (isnan(rec.cut) || hi < rec.cut) && (rec.cut = hi)
    end
    isnan(rec.cut) || (
        rec.to_cut =
            _step_integrals(b, integrator, integrator.tprev, rec.cut, rec.cut_species)
    )
    return nothing
end

# Commits the step after its events: its whole integrals when it ran to its
# end, those up to the crossing when an event ended it there, `NaN` when
# neither accounts for it; then the inventories after it.
function _commit_step!(b::_BoundRecorder, integrator)
    rec = b.rec
    t = integrator.t
    on = rec.reading && integrator.tprev == rec.t_done
    if on && t == rec.t
        add = rec.step
        sp=rec.step_species
    elseif on && t < rec.t && !isnan(rec.cut) && isapprox(t, rec.cut; rtol = 1e-10)
        add = rec.to_cut
        sp=rec.cut_species
    else
        add = (NaN, NaN)
        sp=fill(NaN, length(rec.step_species))
    end
    ng=length(rec.species)
    rec.species_production .+= @view sp[1:ng]
    rec.species_escape .+= @view sp[(ng+1):2ng]
    rec.species_ingress .+= @view sp[(2ng+1):3ng]
    rec.species_final=b.species_inventories(integrator.u, t)
    rec.energy_inflow += add[1]
    rec.mass_outflow += add[2]
    rec.steps += 1
    rec.t_done, rec.reading = t, false
    rec.final = b.inventories(integrator.u, t)
    return nothing
end

# The reader: a continuous callback whose condition never crosses zero, so it
# sees every step's end before the step's events and never fires.
_reader(b::_BoundRecorder) = ContinuousCallback(
    (u, t, integrator) -> (_read_step!(b, u, t, integrator); 1.0),
    integrator -> nothing;
    save_positions = (false, false),
    interp_points = 0,
)

# The recorders among a solve's `callback`, and the rest of it.
_split_recorders(::Nothing) = (ConservationRecorder[], nothing)
_split_recorders(rec::ConservationRecorder) = ([rec], nothing)
_split_recorders(callback) = (ConservationRecorder[], callback)
function _split_recorders(callbacks::Tuple)
    recs = ConservationRecorder[c for c in callbacks if c isa ConservationRecorder]
    rest = Tuple(c for c in callbacks if !(c isa ConservationRecorder))
    return recs,
    isempty(rest) ? nothing : length(rest) == 1 ? only(rest) : CallbackSet(rest...)
end
