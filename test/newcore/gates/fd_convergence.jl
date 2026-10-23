# W10's convergence check of a central-difference derivative (Phase 3 plan
# §4 W10, D144), as pure functions of an objective's values: the gate
# `fd_gradient_qualification.jl` reads its verdict from them and
# `test/newcore/test_fd_convergence.jl` tests them on a synthetic objective.
#
# A parameter p is perturbed relatively, p = p₀ (1 + ε), and the central
# difference at the step h is D(h) = (J(h) − J(−h))/(2h), an estimate of
# dJ/d ln p. Truncation makes D(h) = D + c h² + …, so it settles as h falls,
# and the solve's error ϵ adds about ϵ/h, so it departs again below some
# step. At one tolerance:
#
#   - the agreeing steps: the three largest when their differences lie within
#     `AGREE` of one another (the spread (max − min)/min|D| at most `AGREE`, one
#     sign), extended down through each smaller step whose difference keeps
#     the spread there; none when the three largest disagree;
#   - the converged value: Richardson's (r² D(h) − D(rh))/(r² − 1) from the
#     two smallest agreeing steps h and rh ((4 D(h) − D(2h))/3 on doubling
#     steps), its error estimate |D(h) − D(rh)|/(r² − 1);
#   - the smallest usable step against a reference value: the smallest step
#     whose difference lies within `AGREE` of the reference.
#
# Across two tolerances: at each of the three largest steps, which define the
# converged value, the looser tolerance's difference within `AGREE` of the
# tighter's (a smaller step may be noisier at the looser tolerance; its
# smallest usable step records that). A derivative is qualified when the
# three largest steps agree at the tightest tolerance and the two tightest
# tolerances agree.

const AGREE = 0.05

# The central difference at `h` of the values `J` (a map from ε to J(ε)).
central(J, h) = (J[h] - J[-h]) / (2h)

# The relative spread of the values `d`: (max − min)/min |d|, infinite when
# their signs differ or one is zero.
function spread(d)
    (all(>(0), d) || all(<(0), d)) || return Inf
    return (maximum(d) - minimum(d)) / minimum(abs, d)
end

"""
    convergence(J, steps) -> NamedTuple

At one tolerance: the differences `D` at `steps` (ascending), the agreeing
steps, the converged value and its error estimate (`NaN` when the three
largest disagree).
"""
function convergence(J, steps)
    h = sort(collect(steps))
    D = [central(J, x) for x in h]
    n = length(h)
    n >= 3 || throw(ArgumentError("the check needs at least three steps"))
    agree = spread(D[(n-2):n]) <= AGREE ? collect((n-2):n) : Int[]
    if !isempty(agree)
        k = n - 3
        while k >= 1 && spread(D[vcat(k, agree)]) <= AGREE
            pushfirst!(agree, k)
            k -= 1
        end
    end
    isempty(agree) &&
        return (; steps = h, D, agreeing = Float64[], value = NaN, error = NaN)
    i = agree[1]
    a, b, r2 = D[i], D[i+1], (h[i+1] / h[i])^2
    return (;
        steps = h,
        D,
        agreeing = h[agree],
        value = (r2 * a - b) / (r2 - 1),
        error = abs(a - b) / (r2 - 1),
    )
end

"""
    smallest_usable(c, reference) -> Float64

The smallest step of the check `c` whose difference lies within `AGREE` of
`reference`; `NaN` when none does.
"""
function smallest_usable(c, reference)
    isfinite(reference) || return NaN
    k = findfirst(d -> abs(d - reference) <= AGREE * abs(reference), c.D)
    return k === nothing ? NaN : c.steps[k]
end

"""
    tolerances_agree(loose, tight) -> Bool

Whether, at each of the three largest steps, the looser tolerance's
difference lies within `AGREE` of the tighter's.
"""
function tolerances_agree(loose, tight)
    loose.steps == tight.steps || throw(ArgumentError("the checks use different steps"))
    n = length(tight.steps)
    return all(abs(loose.D[k] - tight.D[k]) <= AGREE * abs(tight.D[k]) for k = (n-2):n)
end

"""
    qualification(checks) -> NamedTuple

The verdict over the checks of one derivative at each tolerance, loosest
first: the three largest steps agreeing at the tightest tolerance and the
two tightest tolerances agreeing; the tightest's converged value; and per
tolerance the smallest usable step against it.
"""
function qualification(checks)
    tight = checks[end]
    agrees = !isempty(tight.agreeing)
    across = length(checks) >= 2 && tolerances_agree(checks[end-1], tight)
    value = tight.value
    return (;
        pass = agrees && across,
        agrees,
        across,
        value,
        error = tight.error,
        usable = [smallest_usable(c, value) for c in checks],
    )
end
