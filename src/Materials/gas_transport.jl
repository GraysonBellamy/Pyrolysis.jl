"""The pore-gas constitutive model, independent of solver state and kernels."""
abstract type GasTransport end

"""Internal specialization for a material with no gaseous components."""
struct NoGasTransport <: GasTransport end

function _nonnegative_finite(x, name)
    isfinite(x) && x >= 0 || throw(ArgumentError("$name must be finite and nonnegative"))
    return Float64(x)
end
function _positive_finite(x, name)
    isfinite(x) && x > 0 || throw(ArgumentError("$name must be finite and positive"))
    return Float64(x)
end

"""
    LumpedEndpoint(; λ, κ, max_skeletal_occupancy=0)

Finite superficial transfer [m²/s] and permeability [m²] when the modeled
matrix disappears. `max_skeletal_occupancy` declares how much remaining
condensed skeletal volume this closure supports; zero describes an empty
channel. The coefficients are physical data, independent of numerical cell size.
"""
struct LumpedEndpoint
    λ::Float64
    κ::Float64
    max_skeletal_occupancy::Float64
    function LumpedEndpoint(; λ, κ, max_skeletal_occupancy = 0.0)
        s = _nonnegative_finite(max_skeletal_occupancy, "max_skeletal_occupancy")
        s < 1 || throw(ArgumentError("max_skeletal_occupancy must be less than one"))
        new(_nonnegative_finite(λ, "λ"), _nonnegative_finite(κ, "κ"), s)
    end
end

"""
    DustyEndpoint(; pore_diameter, permeability, tortuosity=1, max_skeletal_occupancy=0)

Finite channel geometry and permeability for the matrix-absent limit. A
tortuosity of one describes an unobstructed channel only within the declared
skeletal occupancy domain. The material's molecular reference data remain in use.
"""
struct DustyEndpoint
    pore_diameter::Float64
    permeability::Float64
    tortuosity::Float64
    max_skeletal_occupancy::Float64
    function DustyEndpoint(;
        pore_diameter,
        permeability,
        tortuosity = 1.0,
        max_skeletal_occupancy = 0.0,
    )
        s = _nonnegative_finite(max_skeletal_occupancy, "max_skeletal_occupancy")
        s < 1 || throw(ArgumentError("max_skeletal_occupancy must be less than one"))
        τ = _positive_finite(tortuosity, "tortuosity")
        τ >= 1 || throw(ArgumentError("tortuosity must be at least one"))
        new(
            _nonnegative_finite(pore_diameter, "pore_diameter"),
            _nonnegative_finite(permeability, "permeability"),
            τ,
            s,
        )
    end
end

function _transport_blend(blend)
    lo, hi = Float64.(blend)
    0 < lo < hi < 1 || throw(ArgumentError("blend must satisfy 0 < lo < hi < 1"))
    return (lo, hi)
end

"""
    LumpedTransfer(; endpoint=nothing, blend=(1e-4,1e-3))

Superficial partial-pressure transfer with condensed-component `λ`, plus
Darcy advection with condensed-component `κ`. The matrix occupancy is
`q = Σ_condensed γ ξ/ρ_bulk`. Coefficients blend to the explicit endpoint on
`blend`. Without an endpoint the supported domain is `q ≥ blend[2]`;
crossing it terminates integration with `:matrix_exhausted`.
"""
struct LumpedTransfer{E<:Union{Nothing,LumpedEndpoint}} <: GasTransport
    endpoint::E
    blend::Tuple{Float64,Float64}
end
LumpedTransfer(; endpoint::Union{Nothing,LumpedEndpoint} = nothing, blend = (1e-4, 1e-3)) =
    LumpedTransfer(endpoint, _transport_blend(blend))

"""
    DustyGas(; tortuosity, D_ref, T_ref, P_ref, pore_diameter=nothing,
               capillary_regularization=0, endpoint=nothing, blend=(1e-4,1e-3))

Common molecular diffusivity `D_ref*(T/T_ref)^1.75*(P_ref/P)` [m²/s], reduced
by `φ/tortuosity`, and species Knudsen coefficients from pore diameter [m].
An explicit diameter is independent of permeability. Otherwise the effective
capillary geometry is `sqrt(32τκ/φ)`. For a path crossing κ=0, supply an
explicit diameter or a positive `capillary_regularization` [m²], which replaces
`sqrt(κ)` by `κ/sqrt(κ+capillary_regularization)`. This zero-preserving closure
changes positive permeabilities near its scale and requires sensitivity checks.
The matrix domain and endpoint blending have the same meaning as `LumpedTransfer`.
"""
struct DustyGas{D,E<:Union{Nothing,DustyEndpoint}} <: GasTransport
    tortuosity::Float64
    D_ref::Float64
    T_ref::Float64
    P_ref::Float64
    pore_diameter::D
    capillary_regularization::Float64
    endpoint::E
    blend::Tuple{Float64,Float64}
end
function DustyGas(;
    tortuosity,
    D_ref,
    T_ref,
    P_ref,
    pore_diameter = nothing,
    capillary_regularization = 0.0,
    endpoint::Union{Nothing,DustyEndpoint} = nothing,
    blend = (1e-4, 1e-3),
)
    τ = _positive_finite(tortuosity, "tortuosity")
    τ >= 1 || throw(ArgumentError("tortuosity must be at least one"))
    ε = _nonnegative_finite(capillary_regularization, "capillary_regularization")
    pore_diameter !== nothing &&
        ε != 0 &&
        throw(
            ArgumentError(
                "an explicit pore_diameter does not use capillary_regularization",
            ),
        )
    d = pore_diameter === nothing ? nothing : to_property(pore_diameter)
    d !== nothing &&
        _validate_property_law(d, :DustyGas, "pore diameter", "m"; strict_lower = false)
    return DustyGas(
        τ,
        _positive_finite(D_ref, "D_ref"),
        _positive_finite(T_ref, "T_ref"),
        _positive_finite(P_ref, "P_ref"),
        d,
        ε,
        endpoint,
        _transport_blend(blend),
    )
end

has_gas(material) = any(is_gas, material.components)

function _validate_transport_inputs(components, model, mixing)
    any(is_gas, components) || return nothing
    model isa NoGasTransport && throw(
        ArgumentError("NoGasTransport is only valid for materials without gas components"),
    )
    for c in components
        is_gas(c) && continue
        if c.swelling_factor > 0 &&
           model isa LumpedTransfer &&
           c.gas_transfer_coeff === nothing
            throw(
                ArgumentError(
                    "component :$(c.name): specify λ explicitly for a participating lumped gas matrix",
                ),
            )
        end
    end
    for spec in (mixing.λ, mixing.κ)
        spec.rule == BRUGGEMAN &&
            throw(ArgumentError("BRUGGEMAN is a thermal conductivity closure"))
    end
    if model isa DustyGas &&
       model.pore_diameter === nothing &&
       model.capillary_regularization == 0
        cs = filter(c -> !is_gas(c) && c.swelling_factor > 0, components)
        laws = map(c -> c.permeability, cs)
        # Constant, uniformly positive data or an identically closed channel
        # cannot cross the square-root singularity. Other paths need explicit data.
        safe =
            mixing.κ.rule == CARMAN_KOZENY ||
            all(identically_zero, laws) ||
            (!isempty(laws) && all(p -> p isa ConstantProperty && p(300.0) > 0, laws))
        safe || throw(
            ArgumentError(
                "capillary pore diameter may cross κ=0: specify pore_diameter or capillary_regularization",
            ),
        )
    end
    return nothing
end

@inline gas_component_indices(::Tuple{}) = ()
@inline function gas_component_indices(cs::Tuple)
    rest=map(i->i+1, gas_component_indices(Base.tail(cs)))
    return is_gas(first(cs)) ? (1, rest...) : rest
end
