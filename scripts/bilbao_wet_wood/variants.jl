# The comparison's variants: a thermal family (kinetics, shrinkage, back loss)
# whose flux factor is fitted to the dry run (calibrate_dry.jl, written to
# calibration.toml), and the moisture settings run on it.

const CALIBRATION = joinpath(@__DIR__, "calibration.toml")

"""
Thermal families, each calibrated on the dry run's 9 and 19 mm records over
0–24 min (calibrate_dry.jl):
- with the measured 1 mm record as the top boundary (`front = :tc1`), the
  factor on the solids' conductivity (wood and char) is fitted:
  - `tc1`: larch kinetics, 20 % shrinkage of the converted layer, adiabatic back;
  - `tc1_bilbao`: the authors' single-step kinetics and DSC heats instead;
  - `tc1_fixed`: no shrinkage, on a fixed mesh;
- with the measured flux on the whole slab (`front = :flux`), the flux factor
  is fitted:
  - `deck`: larch kinetics, 20 % shrinkage, adiabatic back;
  - `bilbao`, `fixed`: as above;
  - `backloss`: a back loss of 5 W/(m² K) to 300 K through the fibreglass.
"""
const FAMILIES = Dict(
    "tc1" => (front = :tc1, kinetics = :larch, shrinkage = 0.2),
    "tc1_bilbao" => (front = :tc1, kinetics = :bilbao, shrinkage = 0.2),
    "tc1_fixed" => (front = :tc1, kinetics = :larch, shrinkage = nothing),
    "deck" => (front = :flux, kinetics = :larch, shrinkage = 0.2),
    "bilbao" => (front = :flux, kinetics = :bilbao, shrinkage = 0.2),
    "fixed" => (front = :flux, kinetics = :larch, shrinkage = nothing),
    "backloss" => (front = :flux, kinetics = :larch, shrinkage = 0.2, h_back = 5.0),
)

"""The calibrated factors of `family`: (flux_factor, k_scale)."""
function calibrated(family)
    isfile(CALIBRATION) || error("no calibration.toml: run calibrate_dry.jl")
    c = TOML.parsefile(CALIBRATION)
    haskey(c, family) || error("family $family not calibrated")
    return (
        flux_factor = c[family]["flux_factor"],
        k_scale = get(c[family], "k_scale", 1.0),
    )
end

"""
A variant's keywords for `simulate`: `family/law[/tag]`, e.g. `deck/gated`,
`deck/arrhenius/kappa10`. Tags: `kappa0.1`, `kappa10` (both solids'
permeability), `sorption0`, `sorption2` (the sorption part of the heat of
desorption), `kseries` (the water adds nothing to k), `cells200`, `tight`
(reltol 1e-8, dtmax 2 s), `hfront10` (convective loss 10 W/(m² K) at the
exposed face), `flux0.9`, `flux1.1` (the calibrated factor ×0.9, ×1.1).
"""
function variant_kwargs(name::AbstractString)
    parts = split(name, '/')
    family = String(parts[1])
    law = Symbol(parts[2])
    tag = length(parts) > 2 ? String(parts[3]) : ""
    cal = calibrated(family)
    kw = Dict{Symbol,Any}(pairs(FAMILIES[family]))
    kw[:law] = law
    kw[:flux_factor] = cal.flux_factor
    kw[:k_scale] = cal.k_scale
    tag == "kappa0.1" && (kw[:κ_scale] = 0.1)
    tag == "kappa10" && (kw[:κ_scale] = 10.0)
    tag == "sorption0" && (kw[:sorption] = 0.0)
    tag == "sorption2" && (kw[:sorption] = 2.0)
    tag == "kseries" && (kw[:k_rule] = SERIES)
    tag == "cells200" && (kw[:cells] = 200)
    tag == "tight" && (kw[:reltol] = 1e-8; kw[:dtmax] = 2.0)
    tag == "hfront10" && (kw[:h_front] = 10.0)
    tag == "flux0.9" && (kw[:flux_factor] *= 0.9)
    tag == "flux1.1" && (kw[:flux_factor] *= 1.1)
    return (; kw...)
end
