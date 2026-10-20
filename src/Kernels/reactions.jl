# ============================================================================
# Reactions
# ============================================================================

"""
	SWITCH_WIDTH

The reactivity switch is one above `SWITCH_WIDTH` times the deactivation
threshold of a cell's remaining condensed mass fraction. A wider switch slows
the outermost cell so much that the cell inward of it depletes first.
"""
const SWITCH_WIDTH = 1.5

"""
	ReactivitySwitch{NC}

The smooth switch-off of a cell's reactions as its condensed mass approaches
the deactivation threshold: the rates are multiplied by
`s = x² (3 − 2x)` with `x = mˢ / (SWITCH_WIDTH · ε_m · mˢ_ref)`, clamped to
`[0, 1]`, so `s` is one above `SWITCH_WIDTH · ε_m`, `C¹`, about three quarters at
the threshold (which the cell therefore still reaches) and zero at no
condensed mass. On a moving mesh a cell's thickness falls with its mass and its
concentrations do not, so the concentration limiter alone never slows a
vanishing cell.

`reference` holds the cells' initial condensed masses, `solid` flags the
condensed components.
"""
struct ReactivitySwitch{NC}
    threshold::Float64
    reference::Vector{Float64}
    solid::NTuple{NC,Bool}
end

"""
	reactivity(switch, ds, i)

The factor the rates of cell `i` are multiplied by.
"""
@inline function reactivity(sw::ReactivitySwitch{NC}, ds::DerivedState{NC}, i) where {NC}
    ref = sw.reference[i]
    ref > 0 || return one(eltype(ds.m))
    ms = zero(eltype(ds.m))
    @inbounds for j = 1:NC
        sw.solid[j] && (ms += ds.m[j, i])
    end
    x = ms / (SWITCH_WIDTH * sw.threshold * ref)
    x >= 1 && return one(x)
    x <= 0 && return zero(x)
    return x * x * (3 - 2x)
end

"""
	ReactionKernel{S} <: Kernel

Arrhenius kinetics evaluated at the trial temperature: per cell the reaction
mass sources `ωᵢⱼ Vᵢ` and the reaction heat `Q_rxn(Tᵢ) Vᵢ`, with the heats quoted
at the local temperature, times the cell's [`reactivity`](@ref).

The rates carry a [`OneSidedLimiter`](@ref), so an absent or negative trial
concentration is never consumed whatever the reaction order, and every rate is
finite with a finite derivative at any real state; a cell in which every
reaction lacks a reactant contributes nothing and is not evaluated. `switch`
is a [`ReactivitySwitch`](@ref) when cells can be deactivated, `nothing`
otherwise.
"""
struct ReactionKernel{S} <: Kernel
    limiter::OneSidedLimiter
    switch::S
end

"""
	ReactionKernel(material; switch = nothing)

The reaction kernel with the one-sided limiter of the material's threshold and
enabled flag; above the limiter's join its rates equal those of the material's
own limiter.
"""
ReactionKernel(material::Material; switch = nothing) =
    ReactionKernel(OneSidedLimiter(material.depletion_limiter), switch)

# Whether every reaction has a reactant at or below zero in `ξ`: under an
# enabled one-sided limiter its rate and every derivative of it are then
# exactly zero.
@inline _all_quenched(ξ, reactions::Tuple) =
    _quenched(first(reactions), ξ) && _all_quenched(ξ, Base.tail(reactions))
@inline _all_quenched(ξ, ::Tuple{}) = true
@inline _quenched(rxn, ξ) = any(r -> ξ[first(r)] <= 0, rxn.reactants)

function contribute!(
    c::Contributions{NC},
    kernel::ReactionKernel,
    ds::DerivedState{NC},
    topo::MeshTopology,
    material::Material{NC},
    t,
) where {NC}
    @inbounds for i = 1:n_cells(topo)
        ξ = ntuple(j -> ds.ξ[j, i], Val(NC))
        kernel.limiter.enabled && _all_quenched(ξ, material.reactions) && continue
        T = ds.T_trial[i]
        V = ds.V_eval[i]
        evaluate_reactions!(c.rates, c.heats, ξ, T, material, kernel.limiter)
        species_source_terms!(c.species, c.rates, material)
        sV = kernel.switch === nothing ? V : reactivity(kernel.switch, ds, i) * V
        for j = 1:NC
            c.S_mass[j, i] += c.species[j] * sV
        end
        c.S_energy[i] += reaction_heat_source(c.rates, c.heats) * sV
    end
    return c
end
