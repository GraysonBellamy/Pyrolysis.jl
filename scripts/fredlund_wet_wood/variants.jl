# The comparison's variants: the moisture law and one setting changed, on the
# conductivity factors fitted to the oven-dry spruce (calibrate_dry.jl,
# calibration.toml).

const CALIBRATION = joinpath(@__DIR__, "calibration.toml")

"""The fitted factors on the solids' conductivity: (k_scale, kc_scale)."""
function calibrated()
    isfile(CALIBRATION) || error("no calibration.toml: run calibrate_dry.jl")
    c = TOML.parsefile(CALIBRATION)
    return (k_scale = Float64(c["k_scale"]), kc_scale = Float64(c["kc_scale"]))
end

"""
A variant's keywords for `simulate`: `law[/tag]`, e.g. `sorption`,
`sorption/lambda10`, `gated`. The law is `sorption` or `freewater` (the phase
change) or a one-way law of `LAWS`, which runs at the deck's transfer
coefficient whatever the moisture. Tags: `lambda1`, `lambda3`, `lambda10`
(the wood's transfer coefficient, times the deck's 5e-7 m²/s), `kappa0.1`,
`kappa10`, `kappa100` (both solids' permeability), `K1e2`, `K1e4` (the phase
change's exchange rate [1/s]), `dusty4`, `dusty12` (the dusty-gas closure at
that tortuosity), `openback` (the back open to the air), `fluxfront` (the
whole column under the flux of the case's history instead of the 10 mm
record), `uncal` (WOOD1's conductivities unscaled), `cal2` (another point
of the calibration's valley: the wood's unscaled, the char's ×3; 44 K on the
dry run), `cells2x` (0.25 mm cells), `tight` (reltol 1e-8, dtmax 5 s). Tags
join with `+`.
"""
function variant_kwargs(name::AbstractString)
    parts = split(name, '/')
    law = Symbol(parts[1])
    kw = Dict{Symbol,Any}(:law => law, pairs(calibrated())...)
    law in PHASE_CHANGE_LAWS || (kw[:λ_scale] = 1.0)
    for tag in (length(parts) > 1 ? split(parts[2], '+') : String[])
        _apply_tag!(kw, String(tag))
    end
    return kw
end

function _apply_tag!(kw, tag)
    if startswith(tag, "lambda")
        kw[:λ_scale] = parse(Float64, tag[7:end])
    elseif startswith(tag, "kappa")
        kw[:κ_scale] = parse(Float64, tag[6:end])
    elseif startswith(tag, "K1e")
        kw[:exchange_rate] = parse(Float64, tag[2:end])
    elseif startswith(tag, "dusty")
        kw[:tortuosity] = parse(Float64, tag[6:end])
    elseif tag == "openback"
        kw[:openback] = true
    elseif tag == "fluxfront"
        kw[:front] = :flux
    elseif tag == "uncal"
        kw[:k_scale] = 1.0
        kw[:kc_scale] = 1.0
    elseif tag == "cal2"
        kw[:k_scale] = 1.0
        kw[:kc_scale] = 3.0
    elseif tag == "cells2x"
        kw[:cells] = 0   # resolved per case in run_cases.jl
    elseif tag == "tight"
        kw[:reltol] = 1e-8
        kw[:abstol] = 1e-8
        kw[:dtmax] = 5.0
    else
        error("unknown tag $tag")
    end
    return kw
end

"""The keywords of `variant_kwargs` for `case` (resolves `cells2x`)."""
function case_kwargs(name, case; d = deck())
    kw = variant_kwargs(name)
    if get(kw, :cells, nothing) == 0
        kw[:cells] = 2 * round(Int, d["cases"][case]["thickness"] / 5e-4)
    end
    return (; kw...)
end
