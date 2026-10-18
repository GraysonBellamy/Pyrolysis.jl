# ============================================================================
# Reactions
# ============================================================================

"""
	ReactionKernel <: Kernel

Arrhenius kinetics evaluated at the trial temperature: per cell the reaction
mass sources `ωᵢⱼ Vᵢ` and the reaction heat `Q_rxn(Tᵢ) Vᵢ`, with the heats quoted
at the local temperature.

The rates carry a [`OneSidedLimiter`](@ref), so an absent or negative trial
concentration is never consumed whatever the reaction order, and every rate is
finite with a finite derivative at any real state.
"""
struct ReactionKernel <: Kernel
    limiter::OneSidedLimiter
end

"""
	ReactionKernel(material)

The reaction kernel with the one-sided limiter of the material's threshold and
enabled flag; above the limiter's join its rates equal those of the material's
own limiter.
"""
ReactionKernel(material::Material) =
    ReactionKernel(OneSidedLimiter(material.depletion_limiter))

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
