"""
    GasReservoir(; P, T, mole_fractions=NamedTuple())

Complete external ideal-gas state: pressure [Pa], incoming temperature [K],
and normalized mole fractions by tracked gas name. Values may be functions
of time; `mole_fractions` may also be a function returning a NamedTuple.
Only an explicit vacuum (`P=0`) permits empty composition. Unknown species
with nonzero fractions are rejected when bound to a material.
"""
struct GasReservoir{P,T,X}
    P::P
    T::T
    mole_fractions::X
end
function GasReservoir(; P, T, mole_fractions = NamedTuple())
    P isa Real &&
        (!isfinite(P) || P<0) &&
        throw(ArgumentError("reservoir pressure must be finite and nonnegative"))
    T isa Real &&
        (!isfinite(T) || T<=0) &&
        throw(ArgumentError("reservoir temperature must be finite and positive"))
    if mole_fractions isa NamedTuple
        _validate_fractions(mole_fractions, P isa Real && P==0)
    elseif !(mole_fractions isa Function)
        throw(
            ArgumentError(
                "mole_fractions must be a NamedTuple or a time function returning one",
            ),
        )
    end
    return GasReservoir(P, T, mole_fractions)
end

function _validate_fractions(x::NamedTuple, vacuum)
    all(v->isfinite(v) && v>=0, values(x)) ||
        throw(ArgumentError("reservoir fractions must be finite and nonnegative"))
    isempty(x) && vacuum && return nothing
    isapprox(sum(values(x)), 1.0; atol = 32eps(Float64), rtol = 0.0) ||
        throw(ArgumentError("reservoir mole fractions must sum to one"))
    return nothing
end
# Recursive tuple lookup keeps material names dynamic without allocating a
# NamedTuple fallback on every boundary evaluation.
@inline _named_value(x::NamedTuple, name, default) =
    _named_value(keys(x), values(x), name, default)
@inline _named_value(::Tuple{}, ::Tuple{}, name, default) = default
@inline _named_value(names::Tuple, values::Tuple, name, default) =
    first(names)===name ? first(values) :
    _named_value(Base.tail(names), Base.tail(values), name, default)
@inline _fractions(x::NamedTuple, t) = x
@inline _fractions(x::Function, t) = x(t)

"""
    MassTransferFilm(conductance)

Isobaric, common-diffusivity, nonreacting mass film with an explicit
conductance `g=ρD/δ` [kg/(m² s)], constant or a function of surface
temperature. Its integrated law is `N_j=Y_j,s*m + g*β(m/g)*(Y_j,s-Y_j,∞)`
with `β(B)=B/(e^B-1)`. The same conductance for every species is what makes
the species fluxes sum to the total flux. Zero conductance retains donor bulk
flow; infinity is the direct-contact limit. [`AnalogyFilm`](@ref) takes the
conductance from the face's convective condition instead. Thermal exchange is
specified by the thermal condition and the complete reservoir temperature.
"""
struct MassTransferFilm{G}
    conductance::G
    function MassTransferFilm(g::G) where {G}
        g isa Real &&
            (isnan(g) || g<0) &&
            throw(ArgumentError("film conductance must be nonnegative"))
        new{G}(g)
    end
end

"""
    AnalogyFilm(; Le = 1.0)

The twin of the face's convective heat film: the film of
[`MassTransferFilm`](@ref) with its conductance from the heat–mass transfer
analogy, `g = h/(c_p,∞ Le^(2/3))`, where `h` is the bare coefficient of the
`ConvectiveBC` on the same face (summed over the convective members of a
combined condition), `c_p,∞` the heat capacity of the reservoir gas at the
film temperature and `Le` the Lewis number of the exchanged species in that
gas. This is the default film of [`AmbientGasBC`](@ref) and
[`ConvectivePressureBC`](@ref). A face without a convective condition, or
facing a constant vacuum reservoir, has no external film: the reservoir is
then in direct contact with the pore gas, as with `film = nothing`. A
reservoir whose pressure falls to zero in time must declare `film = nothing`
itself.
"""
struct AnalogyFilm
    Le::Float64
    function AnalogyFilm(; Le::Real = 1.0)
        isfinite(Le) && Le>0 ||
            throw(ArgumentError("Lewis number must be finite and positive"))
        new(Float64(Le))
    end
end

const MassFilm = Union{Nothing,MassTransferFilm,AnalogyFilm}

"""
    AmbientGasBC(reservoir; film = AnalogyFilm())

Open pore exchange with a complete reservoir through a mass-transfer film:
by default the [`AnalogyFilm`](@ref) twin of the face's convective condition,
an explicit [`MassTransferFilm`](@ref), or `nothing` for direct contact, which
pins the surface composition to the reservoir's.
"""
struct AmbientGasBC{R<:GasReservoir,F<:MassFilm} <: AbstractMassBC
    reservoir::R
    film::F
end
AmbientGasBC(reservoir::GasReservoir; film::MassFilm = AnalogyFilm()) =
    AmbientGasBC(reservoir, film)

"""
    MassFluxBC(; fluxes, inflow_temperature=NamedTuple())

Species-resolved mass flux [kg/(m² s)], positive outward, as NamedTuples by
gas name. Each possibly inward flux requires its incoming temperature [K].
Both flux and temperature entries may be functions of time. Omitted tracked
species have zero flux. No scalar flux is broadcast across species.
"""
struct MassFluxBC{F<:NamedTuple,T<:NamedTuple} <: AbstractMassBC
    fluxes::F
    inflow_temperature::T
end
function MassFluxBC(; fluxes::NamedTuple, inflow_temperature::NamedTuple = NamedTuple())
    for name in keys(fluxes)
        f=fluxes[name]
        f isa Real && !isfinite(f) && throw(ArgumentError("species flux must be finite"))
        if !(f isa Real && f>=0)
            haskey(inflow_temperature, name) ||
                throw(ArgumentError("inward species :$name requires inflow_temperature"))
        end
    end
    for (name, T) in pairs(inflow_temperature)
        haskey(fluxes, name) ||
            throw(ArgumentError("inflow_temperature :$name has no prescribed flux"))
        T isa Real &&
            (!isfinite(T) || T<=0) &&
            throw(ArgumentError("inflow temperature must be finite and positive"))
    end
    return MassFluxBC(fluxes, inflow_temperature)
end

"""
    NeumannPressureBC(; gradient, T, mole_fractions)

Outward total-pressure gradient [Pa/m], normalized surface composition, and
incoming temperature [K]. The surface pressure is `P_cell+d*gradient` and
must be nonnegative. Species follow their reconstructed partial-pressure
gradients; zero total gradient can retain counterdiffusion.
"""
struct NeumannPressureBC{G,T,X} <: AbstractMassBC
    gradient::G
    T::T
    mole_fractions::X
end
function NeumannPressureBC(; gradient, T, mole_fractions)
    GasReservoir(P = 1.0, T = T, mole_fractions = mole_fractions)
    gradient isa Real &&
        !isfinite(gradient) &&
        throw(ArgumentError("pressure gradient must be finite"))
    return NeumannPressureBC(gradient, T, mole_fractions)
end

"""
    ConvectivePressureBC(reservoir; resistance, flow, film=AnalogyFilm())

Solve `P_surface-P_reservoir=resistance*flow_rate`. Select `flow` explicitly:
`:darcy_velocity` uses m/s (resistance Pa s/m) or `:mass` uses kg/(m² s)
(resistance Pa m² s/kg). Resistance is nonnegative, a constant or a function
of surface temperature like a film conductance: a Darcy skin of thickness `δ`
and permeability `K` has `R_h = μ(T) δ/(ρ_g K)`. The reservoir supplies
composition and incoming temperature; the film solves surface composition by
species-flux continuity.
"""
struct ConvectivePressureBC{R<:GasReservoir,H,F<:MassFilm} <: AbstractMassBC
    reservoir::R
    resistance::H
    flow::Symbol
    film::F
end
function ConvectivePressureBC(
    reservoir::GasReservoir;
    resistance,
    flow::Symbol,
    film::MassFilm = AnalogyFilm(),
)
    flow in (:darcy_velocity, :mass) ||
        throw(ArgumentError("flow must be :darcy_velocity or :mass"))
    resistance isa Real &&
        (!isfinite(resistance) || resistance<0) &&
        throw(ArgumentError("hydraulic resistance must be finite and nonnegative"))
    return ConvectivePressureBC(reservoir, resistance, flow, film)
end

function validate_gas_boundary(bc::Union{AmbientGasBC,ConvectivePressureBC}, material, t)
    _reservoir_values(bc.reservoir, material, t)
    return nothing
end
function validate_gas_boundary(bc::NeumannPressureBC, material, t)
    _reservoir_values(
        GasReservoir(P = 1.0, T = bc.T, mole_fractions = bc.mole_fractions),
        material,
        t,
    )
    return nothing
end
validate_gas_boundary(::ImpermeableBC, material, t) = nothing
function validate_gas_boundary(bc::MassFluxBC, material, t)
    gasnames=map(c->c.name, filter(is_gas, material.components))
    for name in keys(bc.fluxes)
        name in gasnames ||
            throw(ArgumentError("prescribed flux species :$name is not a tracked gas"))
    end
    return nothing
end

@inline function _reservoir_values(res, material, t)
    P=evaluate_param(res.P, t)
    T=evaluate_param(res.T, t)
    isfinite(P) && P>=0 ||
        throw(ArgumentError("reservoir pressure must be finite and nonnegative at t=$t"))
    isfinite(T) && T>0 ||
        throw(ArgumentError("reservoir temperature must be finite and positive at t=$t"))
    x=_fractions(res.mole_fractions, t)
    _validate_fractions(x, P==0)
    for name in keys(x)
        any(c->is_gas(c) && c.name===name, material.components) ||
            x[name]==0 ||
            throw(ArgumentError("reservoir species :$name is not a tracked gas"))
    end
    fractions=map(c->is_gas(c) ? _named_value(x, c.name, 0.0) : 0.0, material.components)
    return (P = P, T = T, x = fractions)
end

"""Integrated mass-film flux, outward, with exact signed zero-conductance limits."""
@inline function mass_film_flux(Ys, Y∞, m, g)
    if g==0
        return (ForwardDiff.value(m)>=0 ? Ys : Y∞) .* m
    end
    isinf(g) &&
        throw(ArgumentError("infinite conductance fixes surface composition directly"))
    return Ys .* m+(g*blowing_factor(m/g)) .* (Ys-Y∞)
end
