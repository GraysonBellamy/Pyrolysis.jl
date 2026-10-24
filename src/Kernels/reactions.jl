# ============================================================================
# Reactions
# ============================================================================

"""
	ReactionKernel <: Kernel

Arrhenius kinetics evaluated at the trial temperature: per cell the reaction
mass sources `ωᵢⱼ Vᵢ` and the reaction heat `Q_rxn(Tᵢ) Vᵢ`, with the heats quoted
at the local temperature.

The rates carry the material's [`DepletionLimiter`](@ref), so an absent or negative trial
concentration is never consumed whatever the reaction order, and every rate is
finite with a finite derivative at any real state; a cell in which every
reaction lacks a reactant contributes nothing and is not evaluated.
"""
struct ReactionKernel <: Kernel
    limiter::DepletionLimiter
end

"""
	ReactionKernel(material)

The reaction kernel with the material's depletion limiter.
"""
ReactionKernel(material::Material) = ReactionKernel(material.depletion_limiter)

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
        for j = 1:NC
            c.S_mass[j, i] += c.species[j] * V
        end
        c.S_energy[i] += reaction_heat_source(c.rates, c.heats) * V
    end
    return c
end
