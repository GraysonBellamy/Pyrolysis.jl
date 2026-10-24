# ============================================================================
# Solid-carried fields
# ============================================================================

"""
	CarriedField{S}

A field carried by the solid. `frame` is `:reference` for a constant per-cell
quantity (an orientation) or `:material` for one convected with the solid;
`packed = true` stores an intensive field `a` as the inventory `Q = mₛ a` so
the solid mass flux transports it and mass-weighted means are conserved;
`source` names the per-cell source of the packed inventory: the carrier's
reaction source `a (dmₛ/dt)_rxn`, and for the deformation gradient also
`mₛ (∇vₛ) F`.
"""
struct CarriedField{S}
    name::Symbol
    frame::Symbol
    packed::Bool
    source::S
    function CarriedField(name::Symbol, frame::Symbol, packed::Bool, source::S) where {S}
        frame in (:reference, :material) ||
            throw(ArgumentError("frame must be :reference or :material, got :$frame"))
        return new{S}(name, frame, packed, source)
    end
end

"""
	CONVERSION_HISTORY

The conversion-history pair of the variable cross section: `G`, the cumulative
dry pyrolysis gas a cell has produced, sourced by the dry-gas reaction rates,
and `M`, the cell's initial dry solid mass, without a source. Both are
extensive and carried by the solid; `χᵢ = Gᵢ / Mᵢ` and the column mean is
`χ̄ = Σ G / Σ M`.
"""
const CONVERSION_HISTORY = (
    CarriedField(:G, :material, false, :dry_gas_production),
    CarriedField(:M, :material, false, :none),
)

"""
	material_axes(F, reference) -> axes

The physical directions of a material's principal axes after the deformation
`F`, for a reference frame whose columns are the axes in priority order
(wood: grain, radial, tangential). The first axis convects exactly as
`F a₁`; the second is `F a₂` projected orthogonal to it; the third completes a
right-handed frame. An orthogonal frame cannot keep two convected directions,
so the priority decides which one is kept.
"""
function material_axes(F::AbstractMatrix, reference::AbstractMatrix)
    d = size(F, 1)
    size(F) == (d, d) == size(reference) ||
        throw(DimensionMismatch("F and the reference frame must be square and alike"))
    if d == 1
        return reshape([one(eltype(F))], 1, 1)
    elseif d == 2
        a1 = F * reference[:, 1]
        a1 /= norm(a1)
        a2 = [-a1[2], a1[1]]
        return hcat(a1, a2)
    elseif d == 3
        a1 = F * reference[:, 1]
        a1 /= norm(a1)
        a2 = F * reference[:, 2]
        a2 -= dot(a2, a1) * a1
        a2 /= norm(a2)
        return hcat(a1, a2, cross(a1, a2))
    end
    throw(ArgumentError("material frames are defined in one to three dimensions"))
end

"""
	physical_tensor(principal, axes) -> SMatrix

The physical tensor `Σᵢ kᵢ êᵢ êᵢᵀ` of a material with principal values
`principal` along the columns of `axes`.
"""
function physical_tensor(principal::AbstractVector, axes::AbstractMatrix)
    d = length(principal)
    size(axes) == (d, d) || throw(DimensionMismatch("one principal value per axis"))
    return SMatrix{d,d}(axes * Diagonal(principal) * axes')
end
