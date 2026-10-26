_property_primal(x) = x
_property_primal(x::ForwardDiff.Dual) = _property_primal(ForwardDiff.value(x))
# Matrix transport data use gas-filtered bulk occupancy. This choice does not
# change the independent thermal-property or skeletal-porosity mixture rules.
@inline _transport_occupancy(c::GaseousComponent, ξ, T) = zero(ξ+T)
@inline _transport_occupancy(c::Union{SolidComponent,LiquidComponent}, ξ, T) =
    c.swelling_factor*ξ/c.density(T)+zero(T)
@inline _matrix_λ(c::GaseousComponent, T) = zero(T)
@inline _matrix_λ(c::Union{SolidComponent,LiquidComponent}, T) =
    c.gas_transfer_coeff === nothing ? zero(T) : c.gas_transfer_coeff(T)

"""Return `(q,w)` for gas-filtered γ-weighted condensed bulk occupancy."""
@inline function transport_matrix_state(ξ::NTuple{NC}, T, material) where {NC}
    qj = map((c, x)->_transport_occupancy(c, x, T), material.components, ξ)
    q = sum(qj)
    w = iszero(_property_primal(q)) ? map(zero, qj) : map(x->x/q, qj)
    return q, w
end

@inline function _transport_mix(w, values, spec)
    parallel = sum(map(*, w, values))
    spec.rule == PARALLEL && return parallel
    spec.rule == WEIGHTED && spec.β == 1 && return parallel
    invseries = zero(parallel)
    closed = false
    for j in eachindex(w)
        if w[j] > 0
            if values[j] == 0
                closed = true
            else
                invseries += w[j]/values[j]
            end
        end
    end
    series = closed || invseries == 0 ? zero(parallel) : inv(invseries)
    spec.rule == SERIES && return series
    spec.rule == WEIGHTED && return spec.β*parallel+(1-spec.β)*series
    throw(ArgumentError("unsupported matrix transport mixing rule $(spec.rule)"))
end

@inline function _matrix_coefficients(w, T, material, φ)
    λ = _transport_mix(w, map(c->_matrix_λ(c, T), material.components), material.mixing.λ)
    κ =
        material.mixing.κ.rule == CARMAN_KOZENY ?
        _carman_kozeny_permeability(φ, material.mixing.κ.ℓ) :
        _transport_mix(
            w,
            map(c->get_permeability(c, T), material.components),
            material.mixing.κ,
        )
    return λ+zero(T+φ), κ+zero(T+φ)
end

@inline function _matrix_blend(q, model)
    lo, hi = model.blend
    q >= hi && return one(q)
    q <= lo && return zero(q)
    s = (q-lo)/(hi-lo)
    return s*s*(3-2s)
end

"""
    transport_domain(ξ,T,material)

Check the physical matrix domain independently of its finite trial-state
extension. Returns `admissible`, `outcome`, matrix occupancy and skeletal
occupancy. Without endpoint data the domain ends at the upper blend boundary.
"""
function transport_domain(ξ, T, material)
    q, _ = transport_matrix_state(ξ, T, material)
    s = skeletal_volume_fraction(ξ, T, material)
    model = material.gas_transport
    !has_gas(material) &&
        return (admissible = true, outcome = :success, q = q, skeletal = s)
    ep = model.endpoint
    outcome =
        q >= model.blend[2] ? :success :
        ep === nothing ? :matrix_exhausted :
        s <= ep.max_skeletal_occupancy ? :success : :unsupported_endpoint_occupancy
    return (admissible = outcome===:success, outcome = outcome, q = q, skeletal = s)
end

# In a trial outside a declared matrix domain the normalized matrix remains
# evaluable, and is zero at no occupancy. Accepted states are checked separately.
@inline _endpoint_coefficients(::Nothing, λ, κ, w) = (λ, κ)
@inline _endpoint_coefficients(e::LumpedEndpoint, λ, κ, w) = (w*λ+(1-w)*e.λ, w*κ+(1-w)*e.κ)
@inline _endpoint_coefficients(e::DustyEndpoint, λ, κ, w) =
    (zero(λ), w*κ+(1-w)*e.permeability)

@inline function transport_coefficients(ξ::NTuple{NC}, T, material, φ) where {NC}
    q, weights = transport_matrix_state(ξ, T, material)
    λm, κm = _matrix_coefficients(weights, T, material, φ)
    model = material.gas_transport
    w = _matrix_blend(q, model)
    λ, κ = _endpoint_coefficients(model.endpoint, λm, κm, w)
    return _transport_coefficients(model, λ, κ, κm, φ, T, w, material.components)
end
@inline function _transport_coefficients(model::LumpedTransfer, λ, κ, κm, φ, T, w, comps)
    K = map(c->is_gas(c) ? λ : zero(λ), comps)
    return (; λ, κ, K, t = zero(λ))
end
@inline function _capillary_diameter(model, κ, φ)
    ε = model.capillary_regularization
    ε > 0 && return sqrt(32model.tortuosity/φ)*κ/sqrt(κ+ε)
    # Exact zero is allowed only for an identically closed channel (setup checks
    # paths crossing zero). Its parameter derivative is intentionally undefined.
    _property_primal(κ) == 0 && return zero(κ+φ)
    return sqrt(32model.tortuosity*κ/φ)
end
@inline function _transport_coefficients(model::DustyGas, λ, κ, κm, φ, T, w, comps)
    τ = model.tortuosity
    dp =
        model.pore_diameter === nothing ? _capillary_diameter(model, κm, φ) :
        model.pore_diameter(T)
    matrix_geometry = φ*dp/(3τ)
    molecular_factor = (model.T_ref/T)^1.75/(model.D_ref*model.P_ref)
    ep = model.endpoint
    geom, t = if ep === nothing
        matrix_geometry, τ/φ*molecular_factor
    else
        w*matrix_geometry+(1-w)*φ*ep.pore_diameter/(3ep.tortuosity),
        (w*τ+(1-w)*ep.tortuosity)/φ*molecular_factor
    end
    K = map(c->is_gas(c) ? geom*sqrt(8R_GAS*T/(π*molar_mass(c))) : zero(geom+T), comps)
    return (; λ = zero(κ), κ, K, t)
end

"""Superficial matrix transfer [m²/s], including the declared endpoint blend."""
@inline effective_gas_transfer(ξ::NTuple, T, material) =
    transport_coefficients(ξ, T, material, 1-skeletal_volume_fraction(ξ, T, material)).λ
"""Effective permeability [m²], including the declared endpoint blend."""
@inline effective_permeability(ξ::NTuple, T, material, φ) =
    transport_coefficients(ξ, T, material, φ).κ

"""Whether both constitutive channels are identically closed for all states."""
function closed_transport_channel(material)
    has_gas(material) || return true
    model=material.gas_transport
    ep=model.endpoint
    ep!==nothing &&
        (
            ep isa LumpedEndpoint ? ep.λ>0 || ep.κ>0 :
            ep.pore_diameter>0 || ep.permeability>0
        ) &&
        return false
    material.mixing.κ.rule==CARMAN_KOZENY && return false
    cs=filter(c->!is_gas(c) && c.swelling_factor>0, material.components)
    all(c->identically_zero(c.permeability), cs) || return false
    if model isa DustyGas
        return model.pore_diameter===nothing || identically_zero(model.pore_diameter)
    end
    return all(
        c->c.gas_transfer_coeff===nothing || identically_zero(c.gas_transfer_coeff),
        cs,
    )
end
