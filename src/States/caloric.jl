# ============================================================================
# Sensible enthalpy of a component
# ============================================================================
#
# hⱼ(T) = ∫_{T_REF}^{T} c_p,ⱼ(T′) dT′ with every component on the same datum.
# A property law evaluates at its argument clamped to the property bracket,
# so c_p is constant outside it and h continues linearly. Constant, linear,
# polynomial and tabulated laws are integrated in closed form, piecewise
# across the bracket ends and table knots, so dh/dT is the law itself to
# round-off. Any other law is integrated once, at setup, onto a 1 K grid and
# evaluated by cubic Hermite interpolation of the integral and the law, which
# keeps h C¹ and gives its derivative in closed form.

# Antiderivative of the clamped law c(clamp(s, lo, hi)) from `lo`, given the
# antiderivative `P` of the raw law and its end values.
@inline function _clamped_antiderivative(P, c_lo, c_hi, T)
    lo, hi = T_PROPERTY_MIN, T_PROPERTY_MAX
    T <= lo && return c_lo * (T - lo)
    T >= hi && return P(hi) - P(lo) + c_hi * (T - hi)
    return P(T) - P(lo)
end

"""
	sensible_enthalpy(law, T) -> h

Sensible enthalpy `∫_{T_REF}^{T} c(T′) dT′` [J/kg] of the enthalpy law `law`: a
constant, linear, polynomial or tabulated heat-capacity law, or a
[`TabulatedEnthalpy`](@ref) built from any other law.
"""
@inline sensible_enthalpy(c::ConstantProperty, T) = c.value * (T - T_REF)

@inline function sensible_enthalpy(c::LinearProperty, T)
    P(s) = s * (c.a + c.b * s / 2)
    G(s) = _clamped_antiderivative(P, c(T_PROPERTY_MIN), c(T_PROPERTY_MAX), s)
    return G(T) - G(T_REF)
end

@inline function sensible_enthalpy(c::PolynomialProperty{N}, T) where {N}
    coeffs = ntuple(k -> c.coeffs[k] / k, Val(N))
    P(s) = s * evalpoly(s, coeffs)
    G(s) = _clamped_antiderivative(P, c(T_PROPERTY_MIN), c(T_PROPERTY_MAX), s)
    return G(T) - G(T_REF)
end

function sensible_enthalpy(c::TableProperty, T)
    return _table_integral(c, T) - _table_integral(c, oftype(float(T), T_REF))
end

# ∫_{T_PROPERTY_MIN}^{T} c(T′) dT′ for a tabulated law: piecewise linear
# between the knots that fall inside the property bracket, constant outside
# the table and outside the bracket.
function _table_integral(c::TableProperty{N}, T) where {N}
    lo, hi = T_PROPERTY_MIN, T_PROPERTY_MAX
    T <= lo && return c(lo) * (T - lo)
    Tc = T >= hi ? oftype(float(T), hi) : T
    s = zero(float(T))
    a = oftype(Tc, lo)
    ca = c(lo)
    # Each piece between consecutive breakpoints is linear in T, so the
    # trapezoid rule is exact on it.
    for k = 1:N
        tk = c.temperatures[k]
        tk <= lo && continue
        tk >= Tc && break
        ck = c(tk)
        s += (tk - a) * (ca + ck) / 2
        a, ca = oftype(Tc, tk), ck
    end
    s += (Tc - a) * (ca + c(Tc)) / 2
    return T >= hi ? s + c(hi) * (T - hi) : s
end

"""
	TabulatedEnthalpy

The sensible enthalpy of a heat-capacity law without a closed-form integral:
the integral at the 1 K knots of the property bracket, accumulated by
five-point Gauss–Legendre quadrature on each interval, and the law at the
knots. Between knots the enthalpy is the cubic Hermite interpolant of these
values and slopes; outside the bracket it continues linearly, as the clamped
law does. [`enthalpy_slope`](@ref) is the interpolant's exact derivative.
Both are `NaN` at a `NaN` temperature, so a wild trial state gives a `NaN`
residual rather than an error.
"""
struct TabulatedEnthalpy
    h::Vector{Float64}
    c::Vector{Float64}
end

const _GL5_X =
    (-0.906179845938664, -0.5384693101056831, 0.0, 0.5384693101056831, 0.906179845938664)
const _GL5_W = (
    0.23692688505618908,
    0.47862867049936647,
    0.5688888888888889,
    0.47862867049936647,
    0.23692688505618908,
)

function TabulatedEnthalpy(c::AbstractPropertyFunction)
    knots = T_PROPERTY_MIN:1.0:T_PROPERTY_MAX
    cs = [Float64(c(T)) for T in knots]
    h = zeros(length(knots))
    for k = 2:length(knots)
        mid = knots[k] - 0.5
        h[k] = h[k-1] + 0.5 * sum(_GL5_W[q] * Float64(c(mid + 0.5 * _GL5_X[q])) for q = 1:5)
    end
    t = TabulatedEnthalpy(h, cs)
    # Shift the datum to T_REF.
    h .-= _tabulated(t, T_REF)
    return t
end

@inline function _tabulated(t::TabulatedEnthalpy, T)
    lo, hi = T_PROPERTY_MIN, T_PROPERTY_MAX
    T <= lo && return t.h[1] + t.c[1] * (T - lo)
    T >= hi && return t.h[end] + t.c[end] * (T - hi)
    isnan(T) && return T * NaN
    k = min(floor(Int, _primal(T) - lo) + 1, length(t.h) - 1)
    s = T - (lo + (k - 1))
    s2 = s * s
    s3 = s2 * s
    return (2s3 - 3s2 + 1) * t.h[k] +
           (s3 - 2s2 + s) * t.c[k] +
           (-2s3 + 3s2) * t.h[k+1] +
           (s3 - s2) * t.c[k+1]
end

@inline sensible_enthalpy(t::TabulatedEnthalpy, T) = _tabulated(t, T)

"""
	enthalpy_slope(law, T)

`dh/dT` of the enthalpy law: the heat-capacity law itself for the closed-form
integrals, the derivative of the Hermite interpolant for a
[`TabulatedEnthalpy`](@ref).
"""
@inline enthalpy_slope(c::AbstractPropertyFunction, T) = c(T)
@inline function enthalpy_slope(t::TabulatedEnthalpy, T)
    lo, hi = T_PROPERTY_MIN, T_PROPERTY_MAX
    T <= lo && return t.c[1] + zero(T)
    T >= hi && return t.c[end] + zero(T)
    isnan(T) && return T * NaN
    k = min(floor(Int, _primal(T) - lo) + 1, length(t.h) - 1)
    s = T - (lo + (k - 1))
    s2 = s * s
    return (6s2 - 6s) * t.h[k] +
           (3s2 - 4s + 1) * t.c[k] +
           (-6s2 + 6s) * t.h[k+1] +
           (3s2 - 2s) * t.c[k+1]
end

"""
	enthalpy_law(c) -> law

The enthalpy law of the heat-capacity law `c`: `c` itself when its integral
has a closed form, a [`TabulatedEnthalpy`](@ref) otherwise.
"""
enthalpy_law(c::Union{ConstantProperty,LinearProperty,PolynomialProperty,TableProperty}) = c
enthalpy_law(c::AbstractPropertyFunction) = TabulatedEnthalpy(c)

"""
	CaloricLaws{NC,H,C}

The heat-capacity laws `c` of a material's components and their enthalpy
laws `h`, built once per problem.
"""
struct CaloricLaws{NC,H<:Tuple,C<:Tuple}
    h::H
    c::C
end

function CaloricLaws(material::Material{NC}) where {NC}
    c = ntuple(j -> material.components[j].heat_capacity, Val(NC))
    h = map(enthalpy_law, c)
    return CaloricLaws{NC,typeof(h),typeof(c)}(h, c)
end

"""
	component_enthalpies(laws, T) -> NTuple{NC}

Sensible enthalpy of every component at temperature `T`.
"""
@inline component_enthalpies(laws::CaloricLaws, T) =
    map(h -> sensible_enthalpy(h, T), laws.h)

"""
	component_enthalpies(laws, T, flags) -> NTuple{NC}

Sensible enthalpy at temperature `T` of each component `flags` marks, and zero
for the others.
"""
@inline component_enthalpies(laws::CaloricLaws, T, flags::Tuple) =
    map((h, flag) -> flag ? sensible_enthalpy(h, T) : zero(T), laws.h, flags)

"""
	component_enthalpy_slopes(laws, T) -> NTuple{NC}

`dhⱼ/dT` of every component's enthalpy law at temperature `T`.
"""
@inline component_enthalpy_slopes(laws::CaloricLaws, T) =
    map(h -> enthalpy_slope(h, T), laws.h)

"""
	component_heat_capacities(laws, T) -> NTuple{NC}

Heat capacity of every component at temperature `T`.
"""
@inline component_heat_capacities(laws::CaloricLaws, T) = map(c -> c(T), laws.c)

# `init + a₁ b₁ + a₂ b₂ + …`, added from left to right. The per-component
# tuples of the laws mix constants and functions of the temperature, so under
# automatic differentiation their elements differ in type, and indexing them in
# a loop would box every element.
@inline sum_products(init, a::Tuple, b::Tuple) =
    sum_products(init + first(a) * first(b), Base.tail(a), Base.tail(b))
@inline sum_products(init, ::Tuple{}, ::Tuple{}) = init

# `init + m₁ c₁(T) + m₂ c₂(T) + …`, added from left to right, without
# evaluating the heat capacity of a component whose mass is exactly zero: its
# term adds an exact zero for every finite heat capacity, and `iszero` of a
# dual number reads its partials too, so the derivative is kept.
@inline heat_capacity_sum(init, laws::CaloricLaws, m::Tuple, T) =
    _heat_capacity_sum(init, m, laws.c, T)
@inline _heat_capacity_sum(acc, m::Tuple, c::Tuple, T) = _heat_capacity_sum(
    iszero(first(m)) ? acc : acc + first(m) * first(c)(T),
    Base.tail(m),
    Base.tail(c),
    T,
)
@inline _heat_capacity_sum(acc, ::Tuple{}, ::Tuple{}, T) = acc

# Smallest heat capacity of a law on the property bracket, sampled at the
# 1 K spacing the component constructors validate the law on. Outside the
# bracket the law is constant, so this bounds it on the whole real line up to
# the variation between samples, which the factor `1/2` in `CaloricExtension`
# absorbs.
function _minimum_heat_capacity(c::AbstractPropertyFunction)
    m = Inf
    for T = T_PROPERTY_MIN:1.0:T_PROPERTY_MAX
        m = min(m, Float64(c(T)))
    end
    return m
end

# ============================================================================
# Extended caloric map of the inventory form
# ============================================================================

"""
	CaloricExtension{NC}

The regularized relation between a cell's enthalpy inventory and its
temperature, defined for every real state:

```
H = Σⱼ mⱼ⁺ hⱼ(T) + κ(m) (T − T_REF)
```

`mⱼ⁺ = positive_part(mⱼ, δᵢ)` is exactly zero at and below zero mass and exactly
`mⱼ` at and above `δᵢ`; `κ` is a C¹ switch that is zero while the guaranteed
heat capacity `Σⱼ mⱼ⁺ cⱼ,min` is at least twice the cell's floor `Cᵢ` and rises
to `Cᵢ` as it vanishes. The relation is strictly increasing in `T` with slope
at least `Cᵢ`, so it has exactly one root for every real `H`. It agrees with
the physical caloric relation on every state in which each species is absent
or holds at least `δᵢ` and the guaranteed capacity is at least `2 Cᵢ`.

A component whose heat capacity is zero (a gas that carries no heat) adds
nothing to the guaranteed capacity, and the switch keeps the slope at least
`Cᵢ`; a negative heat capacity is refused.

Fields: `δ` and `floor` per cell, `c_min` per component (half the sampled
minimum of its heat-capacity law).
"""
struct CaloricExtension{NC}
    δ::Vector{Float64}
    floor::Vector{Float64}
    c_min::NTuple{NC,Float64}
end

function CaloricExtension(
    laws::CaloricLaws{NC},
    δ::Vector{Float64},
    floor::Vector{Float64},
) where {NC}
    length(δ) == length(floor) || throw(DimensionMismatch("δ and floor"))
    all(>(0), δ) || throw(ArgumentError("the positive-part scales must be positive"))
    all(>(0), floor) || throw(ArgumentError("the capacity floors must be positive"))
    c_min = ntuple(j -> 0.5 * _minimum_heat_capacity(laws.c[j]), Val(NC))
    all(>=(0), c_min) || throw(
        ArgumentError("the caloric relation needs every heat capacity to be non-negative"),
    )
    return CaloricExtension{NC}(δ, floor, c_min)
end

@inline function _guaranteed_capacity(ext::CaloricExtension{NC}, m, i) where {NC}
    δ = ext.δ[i]
    C = zero(eltype(m))
    for j = 1:NC
        C += positive_part(m[j], δ) * ext.c_min[j]
    end
    return C
end

"""
	regularization(ext, m, i)

The capacity switch `κ` of cell `i` at masses `m`.
"""
@inline function regularization(ext::CaloricExtension, m, i)
    f = ext.floor[i]
    C = _guaranteed_capacity(ext, m, i)
    C >= 2f && return zero(C)
    s = 1 - C / (2f)
    return f * s * s
end

"""
	extended_energy(ext, laws, T, m, i)

The right-hand side of the extended caloric relation of cell `i`.
"""
@inline function extended_energy(
    ext::CaloricExtension{NC},
    laws::CaloricLaws{NC},
    T,
    m,
    i,
) where {NC}
    h = component_enthalpies(laws, T)
    δ = ext.δ[i]
    E = regularization(ext, m, i) * (T - T_REF)
    for j = 1:NC
        E += positive_part(m[j], δ) * h[j]
    end
    return E
end

"""
	extended_capacity(ext, laws, T, m, i)

`∂/∂T` of [`extended_energy`](@ref): the cell's effective heat capacity, at
least the cell's floor.
"""
@inline function extended_capacity(
    ext::CaloricExtension{NC},
    laws::CaloricLaws{NC},
    T,
    m,
    i,
) where {NC}
    c = component_enthalpy_slopes(laws, T)
    δ = ext.δ[i]
    p = ntuple(j -> positive_part(m[j], δ), Val(NC))
    return sum_products(regularization(ext, m, i), p, c)
end

@inline _primal(x) = x
@inline _primal(x::ForwardDiff.Dual) = _primal(ForwardDiff.value(x))
@inline _depth(x) = 0
@inline _depth(x::ForwardDiff.Dual) = 1 + _depth(ForwardDiff.value(x))

"""
	invert_caloric(ext, laws, H, m, i, T_start; maxiters = 160) -> (T, iterations)

Solve the extended caloric relation of cell `i` for its temperature by a
bracketed Newton iteration started at `T_start` (the previous temperature of
the cell). The bracket follows from `E(T_REF) = 0` and a slope of at least the
guaranteed capacity; a Newton step that does not contract the bracket is
replaced by bisection. Dual-number inventories are differentiated through the
defining relation: the primal root is corrected by one Newton step per nesting
level, which carries the implicit-function derivatives exactly. A non-finite
inventory has temperature `NaN`, so the residual of a trial state carrying one
is `NaN` and the integrator rejects the step, as in the temperature form.
"""
function invert_caloric(
    ext::CaloricExtension,
    laws::CaloricLaws,
    H,
    m,
    i::Int,
    T_start::Float64;
    maxiters::Int = 160,
)
    if H isa ForwardDiff.Dual || eltype(m) <: ForwardDiff.Dual
        mp = ntuple(j -> _primal(m[j]), length(m))
        T0, its = _invert_primal(ext, laws, _primal(H), mp, i, T_start, maxiters)
        T = T0 + zero(H + sum(m))
        for _ = 1:max(_depth(H), maximum(_depth, m))
            T +=
                (H - extended_energy(ext, laws, T, m, i)) /
                extended_capacity(ext, laws, T, m, i)
        end
        return T, its
    end
    return _invert_primal(ext, laws, Float64(H), m, i, T_start, maxiters)
end

function _invert_primal(ext, laws, H::Float64, m, i, T_start, maxiters)
    isfinite(H) && all(isfinite, m) || return NaN, 0
    Cmin = _guaranteed_capacity(ext, m, i) + regularization(ext, m, i)
    H == 0 && return T_REF, 0
    endpoint = T_REF + H / Cmin
    isfinite(endpoint) ||
        throw(OverflowError("caloric bracket of cell $i is not representable"))
    lo, hi = minmax(T_REF, endpoint)
    T =
        lo <= T_start <= hi ? T_start :
        clamp(T_REF + H / extended_capacity(ext, laws, T_REF, m, i), lo, hi)
    for its = 1:maxiters
        r = extended_energy(ext, laws, T, m, i) - H
        C = extended_capacity(ext, laws, T, m, i)
        # Converged when the residual is at the round-off of the inventory, or
        # of one kelvin's worth of it near the datum.
        abs(r) <= 4eps(Float64) * max(abs(H), C * abs(T - T_REF), C) && return T, its
        if r > 0
            hi = T
        else
            lo = T
        end
        if hi - lo <= 4eps(max(abs(lo), abs(hi)))
            # A collapsed bracket holds the root unless the slope bound behind
            # its endpoint was optimistic; then the root lies beyond the
            # endpoint, and the bracket is widened away from the datum.
            abs(r) <= 1e-8 * max(abs(H), C * abs(T - T_REF), C) && return T, its
            if endpoint > T_REF
                hi = T_REF + 2 * (hi - T_REF)
            else
                lo = T_REF + 2 * (lo - T_REF)
            end
            T = lo + (hi - lo) / 2
            continue
        end
        candidate = T - r / C
        margin = 0.05 * (hi - lo)
        T = lo + margin < candidate < hi - margin ? candidate : lo + (hi - lo) / 2
    end
    error(
        "caloric inversion of cell $i did not converge in $maxiters iterations " *
        "(H = $H, m = $m, bracket = ($lo, $hi))",
    )
end
