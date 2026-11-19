"""
Phase change of a liquid component and its vapour: a two-way cell source that
relaxes the pore vapour toward sorption equilibrium with the liquid.

The functions below take scalars and the law objects (activity law, saturation
law, the vapour's molar mass, the exchange rate and the latent heat at
`T_REF`), so every caller evaluates one statement of the law.
"""

using ..Materials: AbstractActivity, AbstractPropertyFunction, activity, sorption_heat

"""Floor [kg/m³] added to the sorbent concentration in [`moisture_content`](@ref)."""
const SORBENT_FLOOR = 1e-3

"""
	moisture_content(ξ_w, S) -> X

Moisture content [kg/kg] of a cell holding the liquid at bulk concentration
`ξ_w` [kg/m³] in sorbing solids of weighted bulk concentration
`S = Σₛ wₛ ξₛ` [kg/m³]:

	X = ξ_w/(S + 10⁻³ kg/m³).

Signed with `ξ_w`. The floor keeps `X` finite in a cell whose sorbing solids
are gone; what liquid such a cell holds then reads as free.
"""
@inline moisture_content(ξ_w, S) = ξ_w / (S + SORBENT_FLOOR)

"""
	equilibrium_vapour_pressure(activity, saturation, X, T) -> p_eq [Pa]

Partial pressure of the vapour in equilibrium with the liquid at moisture
content `X` [kg/kg] and temperature `T` [K]:

	p_eq = a(X, T) p_sat(T),

with `activity` an [`AbstractActivity`](@ref) and `saturation` the liquid's
saturation curve, such as [`WATER_SATURATION`](@ref). Signed with `X`.
"""
@inline equilibrium_vapour_pressure(
    law::AbstractActivity,
    saturation::AbstractPropertyFunction,
    X,
    T,
) = activity(law, X, T) * saturation(T)

"""
	phase_change_rate(activity, saturation, M, K, X, T, φ, ξ_v) -> ω_v [kg/(m³ s)]

Rate at which the liquid becomes vapour, per unit bulk volume:

	ω_v = K (ξ_v^eq − ξ_v),   ξ_v^eq = φ M p_eq(X, T)/(R T),

with `M` [kg/mol] the vapour's molar mass, `K` [1/s] the exchange rate, `φ`
the pore fraction, `ξ_v` [kg/m³] the bulk vapour concentration and `p_eq` from
[`equilibrium_vapour_pressure`](@ref). Positive for evaporation or desorption,
negative for condensation or sorption: its sign is that of `p_eq − p_v`, so
the step cannot evaporate into vapour above its equilibrium pressure.

Both ends restore: a negative trial vapour mass raises the rate, and a
negative trial liquid mass gives a negative activity, which draws vapour back.
"""
@inline function phase_change_rate(
    law::AbstractActivity,
    saturation::AbstractPropertyFunction,
    M,
    K,
    X,
    T,
    φ,
    ξ_v,
)
    ξ_eq = φ * M * equilibrium_vapour_pressure(law, saturation, X, T) / (R_GAS * T)
    return K * (ξ_eq - ξ_v)
end

"""
	phase_change_heat(activity, M, L_ref, X, T, h_v, h_w) -> h_pc [J/kg]

Heat taken per kilogram of liquid that becomes vapour:

	h_pc = L_ref + h_v(T) − h_w(T) + Q_s(X, T),

with `L_ref` [J/kg] the latent heat at `T_REF`, `h_v` and `h_w` [J/kg] the
sensible enthalpies of the vapour and the liquid on the datum `T_REF`, and
`Q_s` the [`sorption_heat`](@ref) of the activity law with `M` [kg/mol] the
vapour's molar mass. `L_ref + h_v − h_w` is Kirchhoff's latent heat: built
from the components' own enthalpies, it keeps the total enthalpy a state
function whatever path the liquid takes. Condensation returns the same heat.
"""
@inline phase_change_heat(law::AbstractActivity, M, L_ref, X, T, h_v, h_w) =
    L_ref + h_v - h_w + sorption_heat(law, X, T, M)
