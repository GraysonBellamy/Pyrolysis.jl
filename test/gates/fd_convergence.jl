# The convergence check of a central-difference derivative, as pure
# functions of an objective's values: `fd_gradient_qualification.jl` reads
# its verdict from them and `test/sensitivity/test_fd_convergence.jl` tests
# them on a synthetic objective.
#
# A parameter p is perturbed relatively, p = p₀ (1 + ε), and the central
# difference at the step h is D(h) = (J(h) − J(−h))/(2h), an estimate of
# dJ/d ln p. Truncation makes D(h) = D + c h² + …, so it settles as h falls,
# and the solve's error ϵ adds about ϵ/h, so it departs again below some
# step. Richardson's extrapolant of two steps h and rh,
# R = (r² D(h) − D(rh))/(r² − 1) ((4 D(h) − D(2h))/3 on doubling steps),
# removes the h² term. At one tolerance:
#
#   - the converged window: of the windows of three consecutive steps whose
#     two extrapolants (of the smaller pair and of the larger) lie within
#     `AGREE` of each other (the spread (max − min)/min |R| at most `AGREE`, one
#     sign), the one where they agree best; none when no window agrees;
#   - the converged value: the extrapolant of the window's smaller pair, its
#     error estimate the difference of the window's two extrapolants;
#   - the usable steps against a reference value: those whose difference lies
#     within `AGREE` of the reference, a step too large for its truncation or
#     too small for the solve's error falling outside.
#
# Across two tolerances: at each step of the tighter tolerance's window, the
# looser tolerance's difference within `AGREE` of the tighter's. A
# derivative is qualified when the tightest tolerance has a converged window
# and the two tightest tolerances agree on it.
#
# The literal reading, reported beside (the `literal` fields): the three
# largest steps' differences within `AGREE` of one another, the set extended
# down through each smaller step that keeps that spread, its value the
# extrapolant of the two smallest steps of the set, and the tolerances
# compared at the three largest steps. It presumes the truncation at the
# largest step is within `AGREE`, which the converged window does not.

const AGREE = 0.05

# The central difference at `h` of the values `J` (a map from ε to J(ε)).
central(J, h) = (J[h] - J[-h]) / (2h)

# The relative spread of the values `d`: (max − min)/min |d|, infinite when
# their signs differ or one is zero.
function spread(d)
    (all(>(0), d) || all(<(0), d)) || return Inf
    return (maximum(d) - minimum(d)) / minimum(abs, d)
end

# Richardson's extrapolant of the differences at steps `h[k]` and `h[k+1]`.
function extrapolant(D, h, k)
    r2 = (h[k+1] / h[k])^2
    return (r2 * D[k] - D[k+1]) / (r2 - 1)
end

"""
    convergence(J, steps) -> NamedTuple

At one tolerance: the differences `D` at `steps` (ascending), the
extrapolants of successive pairs, the converged window, the converged value
and its error estimate (`NaN` without a window), and the literal reading's
agreeing steps and value.
"""
function convergence(J, steps)
    h = sort(collect(steps))
    D = [central(J, x) for x in h]
    n = length(h)
    n >= 3 || throw(ArgumentError("the check needs at least three steps"))
    R = [extrapolant(D, h, k) for k = 1:(n-1)]
    spreads = [spread(R[k:(k+1)]) for k = 1:(n-2)]
    best = argmin(spreads)
    window = spreads[best] <= AGREE ? h[best:(best+2)] : Float64[]
    value, error = isempty(window) ? (NaN, NaN) : (R[best], abs(R[best] - R[best+1]))
    agree = spread(D[(n-2):n]) <= AGREE ? collect((n-2):n) : Int[]
    if !isempty(agree)
        k = n - 3
        while k >= 1 && spread(D[vcat(k, agree)]) <= AGREE
            pushfirst!(agree, k)
            k -= 1
        end
    end
    return (;
        steps = h,
        D,
        extrapolants = R,
        window,
        value,
        error,
        agreeing = h[agree],
        literal_value = isempty(agree) ? NaN : extrapolant(D, h, agree[1]),
    )
end

"""
    usable_steps(c, reference) -> Vector{Float64}

The steps of the check `c` whose difference lies within `AGREE` of
`reference`; none when the reference is not finite.
"""
usable_steps(c, reference) =
    isfinite(reference) ?
    [
        c.steps[k] for
        k in eachindex(c.D) if abs(c.D[k] - reference) <= AGREE * abs(reference)
    ] : Float64[]

"""
    tolerances_agree(loose, tight, steps) -> Bool

Whether, at each of `steps`, the looser tolerance's difference lies within
`AGREE` of the tighter's; false for no steps.
"""
function tolerances_agree(loose, tight, steps)
    loose.steps == tight.steps || throw(ArgumentError("the checks use different steps"))
    isempty(steps) && return false
    return all(steps) do x
        k = findfirst(==(x), tight.steps)
        abs(loose.D[k] - tight.D[k]) <= AGREE * abs(tight.D[k])
    end
end

"""
    qualification(checks) -> NamedTuple

The verdict over the checks of one derivative at each tolerance, loosest
first: a converged window at the tightest tolerance and the two tightest
tolerances agreeing on it; the tightest's converged value, its error and its
window; per tolerance the usable steps against that value; and the literal
reading's verdict and value.
"""
function qualification(checks)
    tight = checks[end]
    loose = length(checks) >= 2 ? checks[end-1] : nothing
    converges = !isempty(tight.window)
    across = loose !== nothing && tolerances_agree(loose, tight, tight.window)
    literal_agrees = !isempty(tight.agreeing)
    literal_across =
        loose !== nothing && tolerances_agree(loose, tight, tight.steps[(end-2):end])
    return (;
        pass = converges && across,
        converges,
        across,
        value = tight.value,
        error = tight.error,
        window = tight.window,
        usable = [usable_steps(c, tight.value) for c in checks],
        literal = (;
            pass = literal_agrees && literal_across,
            agrees = literal_agrees,
            across = literal_across,
            value = tight.literal_value,
        ),
    )
end
