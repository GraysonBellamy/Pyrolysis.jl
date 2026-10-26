# Exact series barriers are discontinuous when their support appears or
# disappears. Check the channels actually used, while retaining invariant
# barriers, identically closed mixtures, and purely parallel weighted rules.
function validate_transport_support(material, ξ0; moving_support = false)
    has_gas(material) || return nothing
    participating=[
        j for j in eachindex(material.components) if
        !is_gas(material.components[j]) && material.components[j].swelling_factor>0
    ]
    for channel in (:λ, :κ)
        channel===:λ && material.gas_transport isa DustyGas && continue
        spec=getproperty(material.mixing, channel)
        (spec.rule==SERIES || (spec.rule==WEIGHTED && spec.β<1)) || continue
        laws=AbstractPropertyFunction[]
        for j in participating
            c=material.components[j]
            law=channel===:λ ? c.gas_transfer_coeff : c.permeability
            law === nothing && throw(
                ArgumentError(
                    "active lumped transport requires an explicit condensed λ coefficient",
                ),
            )
            push!(laws, law)
        end
        isempty(laws) && continue
        zeroids=[participating[k] for k in eachindex(laws) if identically_zero(laws[k])]
        isempty(zeroids) && continue
        length(zeroids)==length(laws) && continue
        for i in axes(ξ0, 2)
            # Reachability, not merely membership in the reaction network:
            # with the depletion limiter an absent precursor cannot react.
            # Moving transport may import any initially present component.
            available=Set(
                j for j in eachindex(material.components) if
                (moving_support ? any(>(0), view(ξ0, j, :)) : ξ0[j, i]>0) ||
                    is_gas(material.components[j])
            )
            changed=Set{Int}()
            again=true
            while again
                again=false
                for rxn in material.reactions
                    rxn.A==0 && continue
                    (
                        material.depletion_limiter.enabled &&
                        !all(pair->first(pair) in available, rxn.reactants)
                    ) && continue
                    for pair in rxn.reactants
                        push!(changed, first(pair))
                    end
                    for (j, yield) in rxn.products
                        yield>0 || continue
                        push!(changed, j)
                        if !(j in available)
                            push!(available, j)
                            again=true
                        end
                    end
                end
            end
            persistent=any(j->ξ0[j, i]>0 && !(j in changed), zeroids) && !moving_support
            persistent && continue
            risky=filter(
                j->(j in changed) || (moving_support && any(>(0), view(ξ0, j, :))),
                zeroids,
            )
            isempty(risky) && continue
            names=map(j->material.components[j].name, risky)
            throw(
                ArgumentError(
                    "active series $channel channel can change support of zero-coefficient component(s) $names in cell $i; specify a continuous physical mixing closure or an explicit support event",
                ),
            )
        end
    end
    return nothing
end
