# ============================================================================
# Smooth guards for trial states
# ============================================================================
#
# Newton iterates leave the physical domain before any step can be rejected,
# so every quantity the residual evaluates is finite with a finite derivative
# at any real state. Each guard below is C¹, is exactly the identity on the
# physical states the tolerances can distinguish, and replaces a hard clamp.
# The admissibility check reads the unguarded quantities.

"""
	smooth_clamp(x, lo, hi, width)

C¹ clamp of `x` into `[lo, hi]`: constant outside the interval, the identity
on `[lo + width, hi − width]`, and a cubic join of zero slope at each end in
between.
"""
@inline function smooth_clamp(x, lo, hi, width)
    x <= lo && return zero(x) + lo
    x >= hi && return zero(x) + hi
    if x < lo + width
        t = (x - lo) / width
        return lo + width * t * t * (2 - t)
    elseif x > hi - width
        t = (hi - x) / width
        return hi - width * t * t * (2 - t)
    end
    return x
end

# Width of the transitions at the ends of the property bracket [K].
const TRIAL_TEMPERATURE_WIDTH = 1.0

"""
	trial_temperature(T)

The temperature every temperature-sensitive kernel evaluates: `T` smoothly
clamped into the property bracket `[T_PROPERTY_MIN, T_PROPERTY_MAX]` with 1 K
transitions. It is `T` on every physical state, so Arrhenius factors, gas
densities, radiation and property laws see the cell temperature, and it keeps
`exp(−E/RT)` finite with a finite derivative at zero and negative trial
temperatures.
"""
@inline trial_temperature(T) =
    smooth_clamp(T, T_PROPERTY_MIN, T_PROPERTY_MAX, TRIAL_TEMPERATURE_WIDTH)

"""
	positive_part(m, δ)

One-sided C¹ positive part: exactly zero at and below zero, exactly `m` at
and above `δ`, and the smoothstep `m (3t² − 2t³)`, `t = m/δ`, between. Values
and slopes match at both joints.
"""
@inline function positive_part(m, δ)
    m <= 0 && return zero(m)
    m >= δ && return m
    t = m / δ
    return m * t * (3t - 2t * t)
end

"""
	smooth_floor(x, floor)

C¹ floor: `floor` at and below `floor`, `x` at and above `2 floor`, and the
cubic `floor (1 + 2t² − t³)`, `t = x/floor − 1`, between.
"""
@inline function smooth_floor(x, floor)
    x <= floor && return zero(x) + floor
    x >= 2floor && return x
    t = (x - floor) / floor
    return floor * (1 + t * t * (2 - t))
end
