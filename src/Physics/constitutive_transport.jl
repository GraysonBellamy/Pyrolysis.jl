"""A diagnosed local transport failure; no unconverged flux is accepted."""
struct TransportFailure <: Exception
    reason::Symbol
    detail::String
end
Base.showerror(io::IO, e::TransportFailure) =
    print(io, "gas transport ", e.reason, ": ", e.detail)
@inline primal(x) = x
@inline primal(x::ForwardDiff.Dual) = primal(ForwardDiff.value(x))

"""
    molar_wall_flux(K,t,p,g)

Solve the common-diffusivity dusty-gas friction equations in O(Ng). `K` is
finite nonnegative wall diffusivity, `t` molecular resistance per pressure,
`p` one nonnegative face pressure vector, and `g=-∇p/(RT)`. With `t=0` this
is the lumped wall law. Single-species molecular friction cancels exactly.
The positive-sum denominator avoids subtraction at strong coupling; the
pressure-scaled branch gives the finite vacuum continuation and derivatives.
"""
@inline function molar_wall_flux(
    K::StaticVector{N},
    t,
    p::StaticVector{N},
    g::StaticVector{N},
) where {N}
    return _molar_wall_flux(K, t, sum(p), p, g, _accurate_sum(g))
end

# P and p describe one face state. The donor-face caller supplies its mean
# total pressure directly, so all diagonal, collective and denominator terms
# use exactly the same pressure rather than a second rounded sum of P*x.
@inline function _molar_wall_flux(
    K::StaticVector{N},
    t,
    P,
    p::StaticVector{N},
    g::StaticVector{N},
    gsum = _accurate_sum(g),
) where {N}
    (N == 1 || primal(t) == 0) && return K .* g
    u=t*P .* K
    all(x->isfinite(primal(x)), u) || throw(
        TransportFailure(
            :coefficient_overflow,
            "t P K exceeds the finite coefficient envelope",
        ),
    )
    b=inv.(1 .+ u)
    a=K .* b
    if maximum(primal, u)<0.5
        c=a .* t .* p
        d=1-sum(c)
    else
        x=p ./ P
        c=(u .* b) .* x
        d=sum(x .* b)
    end
    primal(d)>0 && isfinite(primal(d)) ||
        throw(TransportFailure(:singular_friction, "nonpositive or nonfinite denominator"))
    # At strong coupling, a_j=(1-b_j)/(tP). Removing the common
    # part analytically avoids cancelling O(1/tP) terms when the total
    # pressure gradient is zero and only O(1/(tP)^2) remains.
    nt=all(v->primal(v)>=0.5, u) ? (gsum-_accurate_sum(b .* g))/(t*P*d) :
       _accurate_sum(a .* g)/d
    n=a .* g+c .* nt
    all(x->isfinite(primal(x)), n) || throw(
        TransportFailure(:nonfinite_flux, "finite-state constitutive solve overflowed"),
    )
    return n
end

"""
    pore_half_flux(pl,pr,Tl,Tr,Tpath,d,K,t,mob)

Skeleton-relative molar species flux on a half-cell path oriented left to
right. Wall resistance uses `Tpath`; Darcy density uses the donor temperature.
The collective donor is selected by the diagonal sum at mean pressure, then
normalized without floors. At two-sided vacuum the wall-law continuation has
its finite first derivative. All branch decisions use primal values.
"""
@inline function pore_half_flux(
    pl::StaticVector{N},
    pr::StaticVector{N},
    Tl,
    Tr,
    Tpath,
    d,
    K,
    t,
    mob,
) where {N}
    return _pore_half_flux(
        pl,
        pr,
        Tl,
        Tr,
        Tpath,
        d,
        K,
        t,
        mob,
        pl-pr,
        _pressure_difference(pl, pr),
    )
end

# Pressure increments can be represented more accurately than their sum with
# an absolute interface pressure. Local interface solves pass these differences
# directly; the donor state and the mean pressure still use the same endpoints.
@inline function _pore_half_flux(
    pl::StaticVector{N},
    pr::StaticVector{N},
    Tl,
    Tr,
    Tpath,
    d,
    K,
    t,
    mob,
    Δp,
    ΔP,
) where {N}
    Pl, Pr=sum(pl), sum(pr)
    Pbar=(Pl+Pr)/2
    g=Δp/(R_GAS*Tpath*d)
    gsum=ΔP/(R_GAS*Tpath*d)
    if primal(Pbar)==0 || primal(t)==0 || N==1
        nd=K .* g
    else
        a=K ./ (1 .+ t*Pbar .* K)
        u=t*Pbar .* K
        drive=all(v->primal(v)>=0.5, u) ? gsum-_accurate_sum(g ./ (1 .+ u)) :
              _accurate_sum(a .* g)
        nd =
            primal(drive)>=0 ? _collective_donor(Pbar, pl, Pl, a, g, K, t, gsum) :
            _collective_donor(Pbar, pr, Pr, a, g, K, t, gsum)
    end
    v=mob*ΔP/d
    up=primal(v)>=0
    adv=v .* (up ? pl/(R_GAS*Tl) : pr/(R_GAS*Tr))
    return (n = nd+adv, v = v)
end

@inline function _collective_donor(Pbar, donor, Pup, a, g, K, t, gsum)
    # A nonnegative empty donor implies a zero driving sum and inactive
    # driven channels. This branch preserves the finite vacuum derivative.
    primal(Pup)==0 && return a .* g
    pf=Pbar .* (donor ./ Pup)
    return _molar_wall_flux(K, t, Pbar, pf, g, gsum)
end

# Error-free addition followed by one final rounding. AD follows the same
# algebra; no pressure or composition floor enters these identities.
@inline function _accurate_sum(values)
    total=zero(first(values)+last(values))
    correction=zero(total)
    for value in values
        next=total+value
        correction+=abs(primal(total))>=abs(primal(value)) ? (total-next)+value :
                    (value-next)+total
        total=next
    end
    return total+correction
end
@inline function _pressure_difference(
    left::StaticVector{N},
    right::StaticVector{N},
) where {N}
    z=zero(left[1]+right[1])
    values=SVector{2N}(ntuple(j->j<=N ? left[j]+z : -right[j-N]+z, Val(2N)))
    return _accurate_sum(values)
end
