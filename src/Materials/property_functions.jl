"""
Temperature-dependent property functions.

All property types implement the callable interface: `(p::PropertyType)(T) -> value`
"""

# ============================================================================
# Property Function Types
# ============================================================================

"""
	ConstantProperty{T<:Real}

A property that does not depend on temperature.

# Example
```julia
k = ConstantProperty(0.21)  # W/(m·K)
k(300.0)  # Returns 0.21
k(500.0)  # Returns 0.21
```
"""
struct ConstantProperty{T<:Real} <: AbstractPropertyFunction
    value::T
end

"""
	LinearProperty{T<:Real}

A property that varies linearly with temperature: `p(T) = a + b*T`

# Fields
- `a::T`: Intercept (value at T=0)
- `b::T`: Slope (change per Kelvin)

# Example
```julia
cp = LinearProperty(1420.0, 2.5)  # J/(kg·K)
cp(300.0)  # Returns 1420 + 2.5*300 = 2170
```
"""
struct LinearProperty{T<:Real} <: AbstractPropertyFunction
    a::T  # intercept
    b::T  # slope
end

"""
	PolynomialProperty{N,T<:Real}

A property defined by a polynomial: `p(T) = c₀ + c₁T + c₂T² + ...`

# Fields
- `coeffs::NTuple{N,T}`: Polynomial coefficients (c₀, c₁, ..., cₙ₋₁)

# Example
```julia
cp = PolynomialProperty((1000.0, 2.0, 0.001))  # c₀ + c₁T + c₂T²
cp(300.0)  # Returns 1000 + 2*300 + 0.001*300² = 1690
```
"""
struct PolynomialProperty{N,T<:Real} <: AbstractPropertyFunction
    coeffs::NTuple{N,T}
end

"""
	ArrheniusProperty{T<:Real}

A property following Arrhenius form: `p(T) = A * exp(-E / (R*T))`

Primarily used for reaction rate constants.

# Fields
- `A::T`: Pre-exponential factor
- `E::T`: Activation energy [J/mol]
- `T_ref::T`: Reference temperature (for normalization, optional)

# Example
```julia
k = ArrheniusProperty(8.5e12, 188e3, 298.15)
k(600.0)  # Returns 8.5e12 * exp(-188000 / (8.314 * 600))
```
"""
struct ArrheniusProperty{T<:Real} <: AbstractPropertyFunction
    A::T       # pre-exponential factor
    E::T       # activation energy [J/mol]
    T_ref::T   # reference temperature
end

# Convenience constructors: promote mixed numeric types to a common float type
# (so e.g. `ArrheniusProperty(8.5e12, 188_000)` with an integer E works).
ArrheniusProperty(A::Real, E::Real, T_ref::Real = T_REF) =
    ArrheniusProperty(promote(float(A), float(E), float(T_ref))...)

"""
	TableProperty{N,T<:Real}

A property defined by tabulated data with linear interpolation.

# Fields
- `temperatures::NTuple{N,T}`: Temperature values (must be sorted ascending)
- `values::NTuple{N,T}`: Property values at each temperature

# Example
```julia
k = TableProperty((300.0, 400.0, 500.0), (0.20, 0.22, 0.25))
k(350.0)  # Returns 0.21 (linear interpolation)
```
"""
struct TableProperty{N,T<:Real} <: AbstractPropertyFunction
    temperatures::NTuple{N,T}
    values::NTuple{N,T}

    function TableProperty(temps::NTuple{N,T}, vals::NTuple{N,T}) where {N,T<:Real}
        # Verify temperatures are sorted
        for i = 1:(N-1)
            temps[i] >= temps[i+1] && error("Temperatures must be strictly increasing")
        end
        new{N,T}(temps, vals)
    end
end

"""
	CallableProperty{F}

A property defined by an arbitrary callable (function or functor).

# Fields
- `func::F`: Any callable f(T) -> value

# Example
```julia
# Ideal gas density
ρ_gas = CallableProperty(T -> P_REF * M / (R_GAS * T))
ρ_gas(300.0)  # Returns density at 300 K
```
"""
struct CallableProperty{F} <: AbstractPropertyFunction
    func::F
end

"""
	SutherlandProperty(μ_ref, T_ref, S)
	SutherlandProperty()   # canonical air constants

Sutherland's law for the temperature dependence of a gas transport property:

	μ(T) = μ_ref × (T/T_ref)^(3/2) × (T_ref + S)/(T + S)

The zero-argument form is the canonical Sutherland **air** triple
(μ_ref = 1.716e-5 Pa·s at T_ref = 273.15 K, S = 110.4 K; e.g. White, *Viscous
Fluid Flow*), giving μ(300 K) = 1.846e-5 Pa·s. This is the default carried by
[`Material`](@ref)'s `gas_viscosity` field, so a material that says nothing
about its pore gas behaves exactly as the package did before that field
existed.

# Fields
- `μ_ref`: Reference viscosity at `T_ref` [Pa·s]
- `T_ref`: Correlation reference temperature [K] — **not** the global enthalpy
  datum `T_REF = 298.15 K`; do not conflate the two
- `S`: Sutherland constant [K]

# Note on the low-temperature clamp
This type evaluates the bare correlation and does **not** clamp `T`. The floor
that keeps `T^1.5` evaluable for diverging Newton trial states is applied once,
centrally, in `Properties.effective_gas_viscosity(T, material)`, so it protects
every property function a material might carry — not just this one.

# Example
```julia
SutherlandProperty()(300.0)                  # ≈ 1.846e-5 Pa·s (air)
SutherlandProperty(2.0e-5, 350.0, 120.0)     # a different gas
```
"""
struct SutherlandProperty{T<:Real} <: AbstractPropertyFunction
    μ_ref::T
    T_ref::T
    S::T
end

SutherlandProperty() = SutherlandProperty(1.716e-5, 273.15, 110.4)
SutherlandProperty(μ_ref::Real, T_ref::Real, S::Real) =
    SutherlandProperty(promote(float(μ_ref), float(T_ref), float(S))...)

"""
	StateDependentProperty{F}

A property that depends on both temperature and cell state (concentrations).

This is used for properties like moisture-dependent heat of sorption where
the property value depends on the current composition of the material,
not just temperature.

# Fields
- `func::F`: Any callable f(T, ξ) -> value where ξ is the concentration tuple

# Example
```julia
# Heat of sorption that increases as moisture content decreases
# HoS(MC) = 2275 + 860 * exp(-0.120 * MC) [kJ/kg]
MOISTURE_IDX = 7
heat_sorption = StateDependentProperty((T, ξ) -> begin
	ξ_total = sum(ξ[i] for i in 1:7)  # total solid + liquid mass
	MC = 100.0 * ξ[MOISTURE_IDX] / max(ξ_total, 1.0)  # moisture content %
	return (2275.0 + 860.0 * exp(-0.120 * MC)) * 1e3  # Convert to J/kg
end)
```

# Notes
- When used as a heat of reaction, the concentration ξ is automatically
  passed by the kinetics evaluation functions.
- Calling with just T (without ξ) will raise an error since state is required.
"""
struct StateDependentProperty{F} <: AbstractPropertyFunction
    func::F
end

# ============================================================================
# Callable Interface
# ============================================================================

"""
	property_temperature(T) -> T clamped to [T_PROPERTY_MIN, T_PROPERTY_MAX]

The temperature a property law is evaluated at. Inside the bracket this is
the identity (bit-for-bit, derivative 1); outside it the law sees the
nearest bracket end, so a nonphysical Newton trial temperature cannot reach
a fractional power or a fitted polynomial's far extrapolation. The bracket is
the interval on which the component constructors validate the law.
"""
@inline property_temperature(T) = clamp(T, T_PROPERTY_MIN, T_PROPERTY_MAX)

# Temperature-only evaluation
@inline (p::ConstantProperty)(T) = p.value
@inline (p::LinearProperty)(T) = (Tp = property_temperature(T); p.a + p.b * Tp)
@inline (p::PolynomialProperty)(T) = evalpoly(property_temperature(T), p.coeffs)
@inline (p::ArrheniusProperty)(T) = p.A * exp(-p.E / (R_GAS * property_temperature(T)))
@inline (p::CallableProperty)(T) = p.func(property_temperature(T))
@inline function (p::SutherlandProperty)(T)
    Tp = property_temperature(T)
    return p.μ_ref * (Tp / p.T_ref)^1.5 * (p.T_ref + p.S) / (Tp + p.S)
end
@inline (p::TableProperty{N,T})(Tq::Real) where {N,T} =
    linear_interpolate(p.temperatures, p.values, property_temperature(Tq))

# Temperature and state evaluation (state ignored for non-state-dependent properties)
@inline (p::ConstantProperty)(T, ξ) = p.value
@inline (p::LinearProperty)(T, ξ) = p(T)
@inline (p::PolynomialProperty)(T, ξ) = p(T)
@inline (p::ArrheniusProperty)(T, ξ) = p(T)
@inline (p::CallableProperty)(T, ξ) = p(T)
@inline (p::SutherlandProperty)(T, ξ) = p(T)
@inline (p::TableProperty{N,T})(Tq::Real, ξ) where {N,T} = p(Tq)

# StateDependentProperty requires both T and ξ
@inline (p::StateDependentProperty)(T, ξ) = p.func(property_temperature(T), ξ)

@inline function (p::StateDependentProperty)(T)
    error(
        "StateDependentProperty requires concentration state: use p(T, ξ) instead of p(T)",
    )
end

# ============================================================================
# Helper Functions
# ============================================================================

"""
	linear_interpolate(xs, ys, x)

Perform linear interpolation in a sorted table.
Extrapolates using the nearest segment if x is outside the table range.
"""
@inline function linear_interpolate(xs::NTuple{N,T}, ys::NTuple{N,T}, x::Real) where {N,T}
    # Handle extrapolation
    if x <= xs[1]
        return ys[1]
    elseif x >= xs[N]
        return ys[N]
    end

    # Binary search for the interval
    # For small N, linear search is faster
    if N <= 8
        @inbounds for i = 1:(N-1)
            if x < xs[i+1]
                t = (x - xs[i]) / (xs[i+1] - xs[i])
                return ys[i] + t * (ys[i+1] - ys[i])
            end
        end
    else
        # Binary search for larger tables
        lo, hi = 1, N
        while lo < hi - 1
            mid = (lo + hi) ÷ 2
            if x < xs[mid]
                hi = mid
            else
                lo = mid
            end
        end
        @inbounds begin
            t = (x - xs[lo]) / (xs[hi] - xs[lo])
            return ys[lo] + t * (ys[hi] - ys[lo])
        end
    end

    # Fallback (should not reach here)
    return ys[N]
end

# ============================================================================
# Temperature Derivatives
# ============================================================================

"""
	derivative(p::AbstractPropertyFunction, T) -> Float64

Compute the derivative of a property function with respect to temperature: dp/dT.

Used for derivative and Jacobian computations.

# Example
```julia
cp = LinearProperty(1420.0, 2.5)
derivative(cp, 300.0)  # Returns 2.5 (the slope)
```
"""
function derivative end

# Every law is flat outside the property bracket (its argument is clamped),
# so the analytic derivatives vanish there and `derivative(p, T)` equals
# `ForwardDiff.derivative(p, T)` for every T.
@inline _in_property_bracket(T) = T_PROPERTY_MIN <= T <= T_PROPERTY_MAX

# ConstantProperty: derivative is always zero
@inline derivative(p::ConstantProperty, T) = 0.0
@inline derivative(p::ConstantProperty, T, ξ) = 0.0

# LinearProperty: derivative is the slope
@inline derivative(p::LinearProperty, T) = _in_property_bracket(T) ? p.b : zero(p.b)
@inline derivative(p::LinearProperty, T, ξ) = derivative(p, T)

# PolynomialProperty: derivative of polynomial
# p(T) = c₀ + c₁T + c₂T² + c₃T³ + ...
# p'(T) = c₁ + 2c₂T + 3c₃T² + ...
@inline function derivative(p::PolynomialProperty{N,T}, Tq) where {N,T}
    N == 1 && return 0.0
    _in_property_bracket(Tq) || return 0.0
    # Compute derivative coefficients: (c₁, 2c₂, 3c₃, ...)
    result = p.coeffs[2]
    Tpow = one(Tq)
    @inbounds for i = 3:N
        Tpow *= Tq
        result += (i - 1) * p.coeffs[i] * Tpow
    end
    return result
end
@inline derivative(p::PolynomialProperty, T, ξ) = derivative(p, T)

# ArrheniusProperty: d/dT[A*exp(-E/RT)] = A*exp(-E/RT) * E/(R*T²) = p(T) * E/(R*T²)
@inline function derivative(p::ArrheniusProperty, T)
    _in_property_bracket(T) || return 0.0
    return p(T) * p.E / (R_GAS * T^2)
end
@inline derivative(p::ArrheniusProperty, T, ξ) = derivative(p, T)

# SutherlandProperty: analytic. With μ = A·T^{3/2}/(T+S),
# dμ/dT = μ·(3/(2T) − 1/(T+S)) — cheaper and better-conditioned than expanding
# the quotient rule. Kept in lockstep with the callable above by
# `test/materials/test_gas_viscosity_property.jl`, which gates it against
# `ForwardDiff.derivative`.
@inline derivative(p::SutherlandProperty, T) =
    _in_property_bracket(T) ? p(T) * (1.5 / T - 1.0 / (T + p.S)) : 0.0
@inline derivative(p::SutherlandProperty, T, ξ) = derivative(p, T)

# TableProperty: finite difference approximation of derivative
@inline function derivative(p::TableProperty{N,T}, Tq) where {N,T}
    # Use central difference with small step
    δT = 0.1  # Small temperature perturbation
    return (p(Tq + δT) - p(Tq - δT)) / (2 * δT)
end
@inline derivative(p::TableProperty, T, ξ) = derivative(p, T)

# CallableProperty: finite difference approximation
@inline function derivative(p::CallableProperty, T)
    δT = 0.1
    return (p(T + δT) - p(T - δT)) / (2 * δT)
end
@inline derivative(p::CallableProperty, T, ξ) = derivative(p, T)

# StateDependentProperty: finite difference in T (holding ξ constant)
@inline function derivative(p::StateDependentProperty, T, ξ)
    δT = 0.1
    return (p(T + δT, ξ) - p(T - δT, ξ)) / (2 * δT)
end

# ============================================================================
# Derivative with respect to concentration (∂p/∂ξⱼ)
# ============================================================================

"""
	derivative_ξ(p::AbstractPropertyFunction, T, ξ, j) -> Float64

Compute the partial derivative of property `p` with respect to concentration ξⱼ.

For most properties (temperature-dependent only), this returns 0.
For StateDependentProperty, uses finite difference approximation.
"""
@inline derivative_ξ(p::AbstractPropertyFunction, T, ξ, j::Int) = 0.0
@inline derivative_ξ(p::ConstantProperty, T, ξ, j::Int) = 0.0
@inline derivative_ξ(p::LinearProperty, T, ξ, j::Int) = 0.0
@inline derivative_ξ(p::PolynomialProperty, T, ξ, j::Int) = 0.0
@inline derivative_ξ(p::ArrheniusProperty, T, ξ, j::Int) = 0.0
@inline derivative_ξ(p::TableProperty, T, ξ, j::Int) = 0.0
@inline derivative_ξ(p::CallableProperty, T, ξ, j::Int) = 0.0
@inline derivative_ξ(p::SutherlandProperty, T, ξ, j::Int) = 0.0

# StateDependentProperty: finite difference in ξⱼ (holding T and other ξ constant)
# Note: Uses setindex on mutable copy to avoid ntuple allocation in hot paths
@inline function derivative_ξ(
    p::StateDependentProperty,
    T,
    ξ::NTuple{NC,Float64},
    j::Int,
) where {NC}
    # Use relative perturbation for better numerical accuracy
    # The perturbation should be large enough to capture meaningful changes
    # but small enough for accurate finite differences
    ξ_j = ξ[j]

    # For very small values, use absolute perturbation
    # For larger values, use relative perturbation (1e-4 is better than 1e-6 for FD accuracy)
    δξ = max(abs(ξ_j) * 1e-4, 1e-6)

    # Compute function at perturbed states
    # For efficiency, we evaluate p directly with modified ξⱼ values
    # by creating minimal perturbations using Base.setindex on tuples
    ξ_plus = Base.setindex(ξ, ξ_j + δξ, j)

    # Central difference (or one-sided if at boundary)
    if ξ_j - δξ < 0
        # Forward difference at boundary
        return (p(T, ξ_plus) - p(T, ξ)) / δξ
    else
        # ξ_j - δξ ≥ 0 here, so the perturbed entry stays non-negative
        ξ_minus = Base.setindex(ξ, ξ_j - δξ, j)
        # Central difference
        return (p(T, ξ_plus) - p(T, ξ_minus)) / (2 * δξ)
    end
end

# ============================================================================

"""
	to_property(x)

Convert various inputs to a property function.
- Numbers become ConstantProperty
- Tuples of 2 numbers become LinearProperty
- Functions become CallableProperty
- Property functions pass through unchanged
"""
to_property(x::AbstractPropertyFunction) = x
to_property(x::Real) = ConstantProperty(Float64(x))
to_property(x::Tuple{Real,Real}) = LinearProperty(Float64(x[1]), Float64(x[2]))
to_property(f::Function) = CallableProperty(f)

# ============================================================================
# Unified parameter evaluation
# ============================================================================
#
# Single, AD-compatible interface for evaluating parameters that may be
# constants, functions of time, or functions of temperature. Type-stable —
# preserves Dual number types so the kernels remain AD-clean.

"""
	evaluate_param(param, t)

Evaluate a time-dependent parameter. Handles constants (`Real`) and unary
functions. For constants the AD type of `t` is preserved.
"""
@inline evaluate_param(param::Real, t) = convert(typeof(t), param)
@inline evaluate_param(param::Function, t) = param(t)

"""
	evaluate_param_T(param, T)

Evaluate a temperature-dependent parameter. Handles constants (`Real`),
unary functions, and `AbstractPropertyFunction` callables. AD type of `T`
is preserved for constants.
"""
@inline evaluate_param_T(param::Real, T) = convert(typeof(T), param)
@inline evaluate_param_T(param::Function, T) = param(T)
@inline evaluate_param_T(param::AbstractPropertyFunction, T) = param(T)
# Fast-path: ConstantProperty bypasses the dynamic-dispatch hop
# through param(T) to the value field directly.
@inline evaluate_param_T(param::ConstantProperty, T) = param.value
