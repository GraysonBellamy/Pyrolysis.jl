# The comparison's variants: the moisture law and one setting changed.

"""
A variant's keywords for `simulate`: `law[/tag]`, e.g. `sorption`,
`sorption/lambda1e-7`, `thermakin`. The law is `sorption` or `freewater` (the
phase change) or `thermakin` (the paper's one-way step). Tags: `lambda1e-7`,
`lambda2e-7`, `lambda1.5e-6`, `lambda2e-5` (the transfer coefficient of the
wood stages [m²/s]; the deck's is 5e-7), `kappa0.1`, `kappa10`, `kappa100`
(every solid's permeability), `K1e2`, `K1e4` (the phase change's exchange
rate [1/s]), `dusty4`, `dusty12` (the dusty-gas closure at that tortuosity),
`n2bg` (the nitrogen at the background temperature), `nofoil` (no copper
below the wood), `cells360`, `tight` (reltol 1e-8, dtmax 1 s). Tags join
with `+`: `sorption/kappa10+lambda1.5e-6`.
"""
function variant_kwargs(name::AbstractString)
    parts = split(name, '/')
    kw = Dict{Symbol,Any}(:law => Symbol(parts[1]))
    tags = length(parts) > 1 ? split(parts[2], '+') : String[]
    for tag in tags
        _apply_tag!(kw, String(tag))
    end
    return (; kw...)
end

function _apply_tag!(kw, tag)
    if startswith(tag, "lambda")
        kw[:λ_wood] = parse(Float64, tag[7:end])
    elseif startswith(tag, "kappa")
        kw[:κ_scale] = parse(Float64, tag[6:end])
    elseif startswith(tag, "K1e")
        kw[:exchange_rate] = parse(Float64, tag[2:end])
    elseif startswith(tag, "dusty")
        kw[:tortuosity] = parse(Float64, tag[6:end])
    elseif tag == "n2bg"
        kw[:n2] = :background
    elseif tag == "nofoil"
        kw[:foil] = false
    elseif tag == "cells360"
        kw[:cells] = 360
    elseif tag == "tight"
        kw[:reltol] = 1e-8
        kw[:abstol] = 1e-8
        kw[:dtmax] = 1.0
    else
        error("unknown tag $tag")
    end
    return kw
end
