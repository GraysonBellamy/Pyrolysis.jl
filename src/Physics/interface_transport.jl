# Local systems use only gas species, with statically sized vectors/matrices.
# The returned flux is shared by both cells. No mechanism is averaged separately.

@inline function _lumped_interface(pl, pr, Tl, Tr, Ts, dl, dr, λl, λr, ml, mr)
    if primal(_pressure_difference(pl, pr))>=0
        return _lumped_forward(pl, pr, Tl, Tr, Ts, dl, dr, λl, λr, ml, mr)
    end
    reverse=_lumped_forward(pr, pl, Tr, Tl, Ts, dr, dl, λr, λl, mr, ml)
    return (n = -reverse.n, p = reverse.p, iterations = 0, residual = reverse.residual)
end

@inline function _lumped_forward(pl, pr, Tl, Tr, Ts, dl, dr, λl, λr, ml, mr)
    Pl, Pr=sum(pl), sum(pr)
    Δ=_pressure_difference(pl, pr)
    A=(λl/Tl+ml*Pl/Tl)/dl
    B=(λr/Tr+mr*Pr/Ts)/dr
    C=mr/(Ts*dr)
    # δ=P*-Pr is the nonnegative root Cδ²+(A+B)δ-AΔ=0.
    # Scaling and the conjugate root avoid overflow and cancellation.
    scale=max(primal(A), primal(B), primal(C*Δ))
    if scale==0
        return (
            n = zero(pl+pr),
            p = pl+(dl/(dl+dr))*(pr-pl),
            iterations = 0,
            residual = 0.0,
        )
    end
    a, b, c=A/scale, B/scale, C/scale
    denom=a+b+sqrt((a+b)^2+4c*a*Δ)
    δ=primal(denom)==0 ? zero(Δ) : 2a*Δ/denom
    drop=Δ*(b+c*δ)/(a+b+c*δ)
    vl=ml*drop/dl
    vr=mr*δ/dr
    al=λl/(Tl*dl)
    ar=λr/(Tr*dr)
    ul=vl/Tl
    ur=vr/Ts
    if primal(λl)==0 && primal(λr)==0
        q=A*drop/R_GAS
        return (n = q .* (pl/Pl), p = (Pr+δ) .* (pl/Pl), iterations = 0, residual = 0.0)
    end
    den=al+ar+ur
    if primal(den)==0
        # A closed downstream channel fixes P*=Pl; its undetermined
        # composition is selected by the upstream state. State flux is zero.
        return (n = zero(pl+pr), p = pl, iterations = 0, residual = 0.0)
    end
    ps=((al+ul) .* pl+ar .* pr)/den
    # Elimination evaluates transmitted species flux without subtracting two
    # nearly equal interface pressures across a highly conducting half.
    n=(al*ar .* (pl-pr)+(al*ur+ul*ar+ul*ur) .* pl)/(R_GAS*den)
    return (n = n, p = ps, iterations = 0, residual = 0.0)
end

# Fixed local seed type avoids inference through ForwardDiff's run-time tag
# construction in these nested face solves. Inputs here are primal; outer state
# derivatives are restored separately by implicit_transport_root.
struct LocalTransportTag end
@inline @generated function _local_jacobian(f::F, x::SVector{N,T}) where {F,N,T}
    seeds=Expr[]
    for j = 1:N
        parts=Expr(:tuple, (j==k ? :(one(T)) : :(zero(T)) for k = 1:N)...)
        push!(seeds, :(D(x[$j], ForwardDiff.Partials($parts))))
    end
    entries=[:(ForwardDiff.partials(y[$i])[$j]) for j = 1:N for i = 1:N]
    return quote
        D=ForwardDiff.Dual{LocalTransportTag,T,N}
        z=SVector{N,D}(($(seeds...),))
        y=SVector{N,D}(@inline f(z))
        SMatrix{N,N,T,$(N*N)}(($(entries...),))
    end
end

"""
    bounded_transport_root(f, x0; tolerance=5e-13, maxiters=60)

Newton solve of a dimensionless local residual on the nonnegative orthant.
A merit-function line search and a relative trust bound keep trials in the
physical root domain. A converged residual and nonsingular local Jacobian
are required; exhaustion raises `TransportFailure`. This is a local root
qualification, not a theorem of global uniqueness for arbitrary material data.
"""
# Equilibrate equations and unknowns before solving or estimating condition.
# A trace-species flux can be many orders smaller than its derivative with
# respect to pressure. Dividing its equation by that flux is a useful residual
# test, but the resulting row magnitude is not physical ill-conditioning.
# J δ = r becomes A y = r/rows, with δ = y/columns. No state or composition
# is changed, and the same transformation is used for implicit derivatives.
@inline function _equilibrated_jacobian(J::StaticMatrix{N,N}) where {N}
    rows=SVector{N}(ntuple(i->maximum(abs, J[i, :]), Val(N)))
    all(x->isfinite(x) && x>0, rows) || throw(
        TransportFailure(
            :singular_interface,
            "local Jacobian has a zero or nonfinite equation row",
        ),
    )
    row_scaled=J ./ rows
    columns=SVector{N}(ntuple(j->maximum(abs, row_scaled[:, j]), Val(N)))
    all(x->isfinite(x) && x>0, columns) || throw(
        TransportFailure(
            :singular_interface,
            "local Jacobian has a zero or nonfinite variable column",
        ),
    )
    return (matrix = row_scaled ./ transpose(columns), rows = rows, columns = columns)
end

function bounded_transport_root(
    f::F,
    x0::S;
    tolerance = 5e-13,
    maxiters = 60,
    lower = zero(x0),
    admissible::A = x->true,
) where {F,S<:StaticVector,A}
    x::S=x0
    for iter = 0:maxiters
        r=f(x)
        e=maximum(abs, r)
        isfinite(e) || throw(
            TransportFailure(:nonfinite_local_residual, "local residual is not finite"),
        )
        linear=_equilibrated_jacobian(_local_jacobian(f, x))
        step=try
            (linear.matrix\(r ./ linear.rows)) ./ linear.columns
        catch
            throw(
                TransportFailure(
                    :singular_interface,
                    "local physical root has a singular Jacobian",
                ),
            )
        end
        all(isfinite, step) ||
            throw(TransportFailure(:singular_interface, "local correction is not finite"))
        if e<=tolerance
            inverse=try
                inv(linear.matrix)
            catch
                throw(
                    TransportFailure(
                        :singular_interface,
                        "local physical root has no finite inverse",
                    ),
                )
            end
            condition=opnorm(linear.matrix, Inf)*opnorm(inverse, Inf)
            isfinite(condition) && condition*eps(Float64)<=1/64 || throw(
                TransportFailure(
                    :ill_conditioned_local_root,
                    "equilibrated local Jacobian cannot resolve its implicit derivative in Float64",
                ),
            )
            maximum(abs, step)<=sqrt(tolerance)*max(1.0, maximum(abs, x)) || throw(
                TransportFailure(
                    :unresolved_local_root,
                    "small residual does not resolve the interface state",
                ),
            )
            return (x = x, linear = linear, iterations = iter, residual = e)
        end
        α=one(eltype(x))
        for j in eachindex(x)
            step[j]>0 && (α=min(α, (x[j]-lower[j])/step[j]))
        end
        α=min(α, 2max(maximum(abs, x), 1.0)/max(maximum(abs, step), eps(Float64)))
        accepted=false
        for backtrack = 1:40
            candidate=x-α*step
            if all(candidate .>= lower) && admissible(candidate)
                er=maximum(abs, f(candidate))
                if er<=tolerance || er<(1-1e-4α)*e
                    x=candidate
                    accepted=true
                    break
                end
            end
            α/=2
        end
        accepted || throw(
            TransportFailure(
                :interface_line_search,
                "no decreasing nonnegative Newton step; residual=$e",
            ),
        )
    end
    throw(
        TransportFailure(
            :interface_iterations,
            "local solve exceeded $maxiters iterations",
        ),
    )
end

# Differentiate the converged equations, never the stopping criterion or the
# sequence of damped iterates. This also handles an exact initial root.
@inline function implicit_transport_root(f::F, root::Q) where {F,Q}
    residual=f(root.x)
    linear=root.linear
    correction=-(linear.matrix\(residual ./ linear.rows)) ./ linear.columns
    return root.x+map(z->z-primal(z), correction)
end

@inline _primal_coefficients(c) = map(x->x isa StaticVector ? map(primal, x) : primal(x), c)

@inline function _interface_halves(offset, base, pl, pr, Tl, Tr, Ts, dl, dr, cl, cr)
    ps=base+offset
    left=_pore_half_flux(
        pl,
        ps,
        Tl,
        Ts,
        Tl,
        dl,
        cl.K,
        cl.t,
        cl.mob,
        (pl-base)-offset,
        _pressure_difference(pl, base)-_accurate_sum(offset),
    ).n
    right=_pore_half_flux(
        ps,
        pr,
        Ts,
        Tr,
        Tr,
        dr,
        cr.K,
        cr.t,
        cr.mob,
        (base-pr)+offset,
        _pressure_difference(base, pr)+_accurate_sum(offset),
    ).n
    return left, right
end

function _interface_residual(
    offset::StaticVector{N},
    base,
    pl,
    pr,
    Tl,
    Tr,
    Ts,
    dl,
    dr,
    cl,
    cr,
    scale,
    pressure_scale,
) where {N}
    left, right=_interface_halves(offset, base, pl, pr, Tl, Tr, Ts, dl, dr, cl, cr)
    # A blocked species has no continuity equation. The base pressure is the
    # distance-weighted state; zero offset pins its unobservable value.
    return SVector{N}(
        ntuple(Val(N)) do j
            if primal(cl.K[j])==0 &&
               primal(cr.K[j])==0 &&
               primal(cl.mob)==0 &&
               primal(cr.mob)==0
                offset[j]/pressure_scale
            else
                (left[j]-right[j])/scale[j]
            end
        end,
    )
end

"""
    interface_flux(pl,pr,Tl,Tr,dl,dr,cl,cr)

Match total molar species flux through two half cells. Coefficients contain
`K` (gas species only), `t` and Darcy `mob`. The common interface temperature
is `(dr*Tl+dl*Tr)/(dl+dr)`. The lumped and pure-gas cases have an analytic
positive quadratic root. Other cases use a scaled, nonnegative local solve
and the implicit-function derivative; failures are explicit.
"""
function interface_flux(
    pl::StaticVector{N},
    pr::StaticVector{N},
    Tl,
    Tr,
    dl,
    dr,
    cl,
    cr,
) where {N}
    Ts=(dr*Tl+dl*Tr)/(dl+dr)
    if N==1 || (
        primal(cl.t)==0 &&
        primal(cr.t)==0 &&
        all(==(cl.K[1]), cl.K) &&
        all(==(cr.K[1]), cr.K)
    )
        return _lumped_interface(
            pl,
            pr,
            Tl,
            Tr,
            Ts,
            dl,
            dr,
            cl.K[1],
            cr.K[1],
            cl.mob,
            cr.mob,
        )
    end
    # At complete vacuum the collective and viscous terms are second order.
    # The linear wall continuation carries the exact finite state Jacobian.
    if primal(sum(pl)+sum(pr))==0
        linearL=cl.K/(R_GAS*Tl*dl)
        linearR=cr.K/(R_GAS*Tr*dr)
        n=map(
            (a, b, l, r)->primal(a+b)==0 ? zero(l+r) : a*b/(a+b)*(l-r),
            linearL,
            linearR,
            pl,
            pr,
        )
        p=map(
            (a, b, l, r)->primal(a+b)==0 ? (dr*l+dl*r)/(dl+dr) : (a*l+b*r)/(a+b),
            linearL,
            linearR,
            pl,
            pr,
        )
        return (n = n, p = p, iterations = 0, residual = 0.0)
    end
    l, r=map(primal, pl), map(primal, pr)
    tl, tr, dL, dR=primal(Tl), primal(Tr), primal(dl), primal(dr)
    ts=(dR*tl+dL*tr)/(dL+dR)
    cL, cR=_primal_coefficients(cl), _primal_coefficients(cr)
    pscale=max(sum(l), sum(r))
    # This interpolation preserves an identical state exactly; multiplying
    # both absolute pressures by distances can manufacture an ulp-sized jump.
    base=l+(dL/(dL+dR))*(r-l)
    guess=zero(base)
    initialL, initialR=_interface_halves(guess, base, l, r, tl, tr, ts, dL, dR, cL, cR)
    # Scaling tracks transmitted flux, rather than K*P, at strong isobaric
    # coupling. Offsets retain small pressure jumps without adding them to P.
    initialscales=_interface_flux_scales(initialL, initialR)
    firstroot=_solve_interface_at_scale(
        base,
        l,
        r,
        tl,
        tr,
        ts,
        dL,
        dR,
        cL,
        cR,
        initialscales,
        pscale,
        guess,
    )
    matchedL, matchedR=_interface_halves(
        firstroot.x*pscale,
        base,
        l,
        r,
        tl,
        tr,
        ts,
        dL,
        dR,
        cL,
        cR,
    )
    scales=_interface_flux_scales(matchedL, matchedR)
    root=_solve_interface_at_scale(
        base,
        l,
        r,
        tl,
        tr,
        ts,
        dL,
        dR,
        cL,
        cR,
        scales,
        pscale,
        firstroot.x,
    )
    f=z->_interface_residual(
        z*pscale,
        base,
        pl,
        pr,
        Tl,
        Tr,
        Ts,
        dl,
        dr,
        cl,
        cr,
        scales,
        pscale,
    )
    offset=implicit_transport_root(f, root)*pscale
    ps=base+offset
    nl, nr=_interface_halves(offset, base, pl, pr, Tl, Tr, Ts, dl, dr, cl, cr)
    mismatch=maximum(abs, map(primal, nl-nr))
    fluxscale=max(maximum(abs, map(primal, nl)), maximum(abs, map(primal, nr)))
    mismatch<=1e-8*fluxscale || throw(
        TransportFailure(
            :unresolved_flux_continuity,
            "the two half paths do not match within 1e-8 of the transmitted flux; the interface is unresolved in Float64",
        ),
    )
    return (n = (nl+nr)/2, p = ps, iterations = root.iterations, residual = root.residual)
end

@inline function _interface_flux_scales(left, right)
    raw=max.(abs.(left), abs.(right))
    common=maximum(raw)
    return map(x->x==0 ? (common==0 ? 1.0 : common) : x, raw)
end
@inline function _solve_interface_at_scale(
    base,
    l,
    r,
    tl,
    tr,
    ts,
    dl,
    dr,
    cl,
    cr,
    scales,
    pscale,
    guess,
)
    f=z->_interface_residual(
        z*pscale,
        base,
        l,
        r,
        tl,
        tr,
        ts,
        dl,
        dr,
        cl,
        cr,
        scales,
        pscale,
    )
    # Division in the normalized lower bound and reconstruction can round in
    # opposite directions. Qualify the physical pressures themselves as well.
    admissible=z->all(base+z*pscale .>= 0)
    return bounded_transport_root(
        f,
        guess;
        lower = -base/pscale,
        admissible,
        tolerance = 2e-10,
    )
end
