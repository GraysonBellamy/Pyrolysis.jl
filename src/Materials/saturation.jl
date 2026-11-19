"""
Saturation pressure of a liquid over its vapour.
"""

# Critical point of water (Wagner and Pruß 2002).
const WATER_T_CRITICAL = 647.096   # K
const WATER_P_CRITICAL = 22.064e6  # Pa

# Range of the saturation curve [K] and the width of its smooth ends.
const WATER_SATURATION_T_MIN = 255.0
const WATER_SATURATION_T_MAX = 640.0
const WATER_SATURATION_WIDTH = 1.0

"""
	WaterSaturation()

Saturation pressure of water [Pa] at temperature `T` [K]: the vapour-pressure
equation of Wagner and Pruß (2002), the auxiliary equation of the IAPWS-95
formulation,

	ln(p_sat/p_c) = (T_c/T) (a₁τ + a₂τ^1.5 + a₃τ³ + a₄τ^3.5 + a₅τ⁴ + a₆τ^7.5),   τ = 1 − T/T_c,

with `T_c = 647.096 K`, `p_c = 22.064 MPa` and
`a = (−7.85951783, 1.84408259, −11.7866497, 22.6807411, −15.9618719, 1.80122502)`.
It gives 611.66 Pa at the triple point, 3169.8 Pa at 298.15 K and 101 418 Pa
at 373.15 K.

The temperature is held in `[255, 640]` K by `smooth_clamp` with 1 K
transitions, so the curve is C¹ at every trial temperature and constant
outside that range. [`WATER_SATURATION`](@ref) is the instance; another
liquid's saturation curve is any `AbstractPropertyFunction`.
"""
struct WaterSaturation <: AbstractPropertyFunction end

"""
	WATER_SATURATION

The saturation pressure of water, [`WaterSaturation`](@ref)`()`.

```julia
WATER_SATURATION(373.15)   # ≈ 101 418 Pa
```
"""
const WATER_SATURATION = WaterSaturation()

@inline function (::WaterSaturation)(T)
    Tc = WATER_T_CRITICAL
    Tb = smooth_clamp(
        T,
        WATER_SATURATION_T_MIN,
        WATER_SATURATION_T_MAX,
        WATER_SATURATION_WIDTH,
    )
    τ = 1 - Tb / Tc
    s =
        -7.85951783τ + 1.84408259τ^1.5 - 11.7866497τ^3 + 22.6807411τ^3.5 - 15.9618719τ^4 +
        1.80122502τ^7.5
    return WATER_P_CRITICAL * exp(Tc / Tb * s)
end
