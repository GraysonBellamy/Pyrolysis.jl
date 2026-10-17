# ============================================================================
# Reactions
# ============================================================================

"""
	ReactionKernel <: Kernel

Arrhenius kinetics evaluated at the trial temperature: per cell the reaction
mass sources `ωᵢⱼ Vᵢ` and the reaction heat `Q_rxn(Tᵢ) Vᵢ`, with the heats quoted
at the local temperature and the material's depletion limiter on the rates.
"""
struct ReactionKernel <: Kernel end

function contribute!(
    c::Contributions{NC},
    ::ReactionKernel,
    ds::DerivedState{NC},
    topo::MeshTopology,
    material::Material{NC},
    t,
) where {NC}
    @inbounds for i = 1:n_cells(topo)
        ξ = ntuple(j -> ds.ξ[j, i], Val(NC))
        T = ds.T_trial[i]
        V = ds.V_eval[i]
        evaluate_reactions!(c.rates, c.heats, ξ, T, material)
        species_source_terms!(c.species, c.rates, material)
        for j = 1:NC
            c.S_mass[j, i] += c.species[j] * V
        end
        c.S_energy[i] += reaction_heat_source(c.rates, c.heats) * V
    end
    return c
end
