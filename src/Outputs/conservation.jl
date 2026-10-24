# ============================================================================
# Conservation: the ledger of the live region
# ============================================================================

"""
	ledger_rates(def, c, ds) -> (energy, outflow)

The rates the live region's inventories follow, from the derived state `ds`
and the kernel contributions `c` of one residual evaluation: the enthalpy rate
[W] (the energy entering through the boundary faces of the live cells, by
conduction or the boundary condition and with the enthalpy their mass fluxes
carry, plus the live cells' reaction heat and formation terms `Σⱼ ωⱼ hⱼ`) and
the mass outflow through those faces [kg/s].
"""
function ledger_rates(def::CoreDefinition, c, ds)
    topo = def.setup.topology
    g = ds.geometry
    E = zero(eltype(c.F_energy))
    out = zero(eltype(c.F_mass))
    for f0 = 1:n_faces(topo)
        topo.faces[f0].cells[2] == 0 || continue
        f = live_face(topo, ds.active, f0)
        f == 0 && continue
        E -= c.F_energy[f] * g.face_area[f]
        for j in axes(c.F_mass, 1)
            out += c.F_mass[j, f] * g.face_area[f]
        end
    end
    for i = 1:n_cells(topo)
        ds.active[i] || continue
        E += c.S_energy[i]
        for j in axes(c.S_mass, 1)
            E += c.S_mass[j, i] * ds.h[j, i]
        end
    end
    return E, out
end

"""
	live_inventories(def, ext, ds, u) -> (energy, mass)

The live cells' energy [J] and mass [kg] in the state `u` with its derived
state `ds`: the enthalpy inventories in the inventory form; in the temperature
form the extended energy of the cells' masses at their temperatures under the
caloric extension `ext`, which equals their sensible enthalpy `Σⱼ mⱼ hⱼ(T)`
away from empty cells and is what a deactivation event conserves.
"""
function live_inventories(def::CoreDefinition, ext::CaloricExtension, ds, u)
    layout, laws = def.layout, def.setup.laws
    NC = n_components(layout)
    E, M = 0.0, 0.0
    for i = 1:n_cells(layout)
        ds.active[i] || continue
        m = ntuple(j -> u[mass_index(layout, i, j)], NC)
        E +=
            energy_state(layout) isa EnthalpyState ? u[energy_index(layout, i)] :
            extended_energy(ext, laws, ds.T[i], m, i)
        M += sum(m)
    end
    return E, M
end

# The caloric extension of a solve's cells.
_extension(def::CoreDefinition) = CaloricExtension(
    def.setup.laws,
    def.setup.reference.mass_scale,
    def.setup.reference.capacity_floor,
)

"""
	ConservationRecorder()

A read-only observer to pass as `solve`'s `callback` (alone, or in a tuple
with other callbacks): over each accepted step it integrates the live
region's mass outflow and enthalpy rate ([`ledger_rates`](@ref)) by
three-point Gauss–Legendre quadrature on the step's interpolant, up to the
time where an event ends the step early, and it records the live inventories
([`live_inventories`](@ref)) at the start and after every step.
[`conservation_audit`](@ref) closes them against each other. A step it
cannot account for (one a callback of the caller's ends early) leaves its
integrals `NaN`.

Fields: `mass_outflow` [kg], `energy_inflow` [J], `initial` and `final`
(each `(energy, mass)`, [J] and [kg]), `steps`.
"""
mutable struct ConservationRecorder
    mass_outflow::Float64
    energy_inflow::Float64
    initial::NTuple{2,Float64}
    final::NTuple{2,Float64}
    steps::Int
    # The step read before its events: its ends, its integrals, and the time
    # and integrals up to the first crossing of an event condition inside it.
    reading::Bool
    t::Float64
    tprev::Float64
    step::NTuple{2,Float64}
    cut::Float64
    to_cut::NTuple{2,Float64}
    t_done::Float64
end

ConservationRecorder() = ConservationRecorder(
    0.0,
    0.0,
    (NaN, NaN),
    (NaN, NaN),
    0,
    false,
    NaN,
    NaN,
    (0.0, 0.0),
    NaN,
    (0.0, 0.0),
    NaN,
)

"""
	conservation_audit(sol[, recorder]) -> NamedTuple

The mass and energy closures of the solve `sol`: per quantity `(; initial,
final, change, inflow, lost, defect, relative)`, with `change = final −
initial` of the live inventories ([`live_inventories`](@ref)), `inflow` the
time integral of the rate entering them (the negated mass outflow, the
enthalpy rate), `lost` what deactivation events took out of them (their lost
masses; their lost enthalpy and regularization), `defect = change − inflow +
lost`, and `relative` the defect over the largest of `|change|`, `|inflow|`
and `|lost|`. With the [`ConservationRecorder`](@ref) that observed the solve,
from its start to its end over the accepted steps; without one, from the
first to the last saved state by the trapezoid over the saved states, whose
error is the saving interval's. `method` names which.
"""
function conservation_audit(
    sol::PyrolysisSolution,
    rec::Union{Nothing,ConservationRecorder} = nothing,
)
    x = sol.extras
    nt = length(sol.t)
    nt >= 1 || throw(ArgumentError("the solution saved no state"))
    if rec === nothing
        out, inflow = 0.0, 0.0
        for k = 2:nt
            dt = sol.t[k] - sol.t[k-1]
            out += (x.mass_outflow[k-1] + x.mass_outflow[k]) * dt / 2
            inflow += (x.energy_rate[k-1] + x.energy_rate[k]) * dt / 2
        end
        initial = (x.live_energy[1], x.live_mass[1])
        final = (x.live_energy[end], x.live_mass[end])
        events = [e for e in sol.events if sol.t[1] < e.t <= sol.t[end]]
        method = :trapezoid
    else
        out, inflow = rec.mass_outflow, rec.energy_inflow
        initial, final = rec.initial, rec.final
        events = sol.events
        method = :recorder
    end
    lost_mass = sum((sum(e.lost_mass) for e in events); init = 0.0)
    lost_energy = sum((e.lost_enthalpy + e.regularization for e in events); init = 0.0)
    closure(a, b, inflow, lost) = begin
        change = b - a
        defect = change - inflow + lost
        scale = max(abs(change), abs(inflow), abs(lost))
        (; initial = a, final = b, change, inflow, lost, defect, relative = defect / scale)
    end
    return (;
        method,
        mass = closure(initial[2], final[2], -out, lost_mass),
        energy = closure(initial[1], final[1], inflow, lost_energy),
    )
end
