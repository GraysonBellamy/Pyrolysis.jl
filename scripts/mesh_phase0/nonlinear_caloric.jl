module NonlinearCaloric

using ForwardDiff
using ..MeshPhase0: positive_mass

export PolynomialLaw,
    PeakLaw,
    CaloricModel,
    heat_capacity,
    enthalpy,
    extended_energy,
    extended_capacity,
    regularization,
    invert_temperature,
    inverse_report

# Manufactured laws, not fits to production materials. Both have a known,
# strictly positive global lower bound and an analytic enthalpy on the real line.
struct PolynomialLaw
    a::Float64
    b::Float64
    scale::Float64
    function PolynomialLaw(a, b, scale)
        a > 0 && b >= 0 && scale > 0 || throw(ArgumentError("invalid cp law"))
        new(a, b, scale)
    end
end
struct PeakLaw
    a::Float64
    b::Float64
    center::Float64
    width::Float64
    function PeakLaw(a, b, center, width)
        a > 0 && b >= 0 && width > 0 || throw(ArgumentError("invalid cp law"))
        new(a, b, center, width)
    end
end
heat_capacity(l::PolynomialLaw, T) = l.a + l.b * (T / l.scale)^2
heat_capacity(l::PeakLaw, T) = l.a + l.b / (1 + ((T - l.center) / l.width)^2)
enthalpy(l::PolynomialLaw, T, Tr) =
    (T - Tr) * (l.a + l.b * (T^2 + T * Tr + Tr^2) / (3l.scale^2))
enthalpy(l::PeakLaw, T, Tr) =
    l.a * (T - Tr) +
    l.b * l.width * (atan((T - l.center) / l.width) - atan((Tr - l.center) / l.width))

struct CaloricModel{L}
    laws::L
    delta::Float64
    floor::Float64
    Tref::Float64
    function CaloricModel(laws; delta = 1e-10, floor = 1e-12, Tref = 300.0)
        !isempty(laws) && delta > 0 && floor > 0 || throw(ArgumentError("invalid model"))
        new{typeof(laws)}(laws, delta, floor, Tref)
    end
end

# kappa is chosen from a GLOBAL lower capacity bound, not cp at a guessed T.
# Exact physical agreement requires both identity masses and inactive kappa.
# A pointwise floor in a denominator alone does not define this caloric map.
function regularization(m, model)
    length(m) == length(model.laws) || throw(DimensionMismatch("masses and laws"))
    Cmin = sum(positive_mass(m[i], model.delta) * l.a for (i, l) in enumerate(model.laws))
    # Nonnegative C1 switch, inactive at 2*floor. Unlike subtracting Cmin from
    # an arbitrary denominator floor, this cannot yield a negative kappa.
    return Cmin < 2model.floor ? model.floor * (1 - Cmin / (2model.floor))^2 : zero(Cmin)
end
function extended_energy(T, m, model)
    kappa = regularization(m, model)
    return sum(
        positive_mass(m[i], model.delta) * enthalpy(l, T, model.Tref) for
        (i, l) in enumerate(model.laws)
    ) + kappa * (T - model.Tref)
end
function extended_capacity(T, m, model)
    return sum(
        positive_mass(m[i], model.delta) * heat_capacity(l, T) for
        (i, l) in enumerate(model.laws)
    ) + regularization(m, model)
end

# Safeguarded Newton; no persistent guesses, so rejected steps and AD calls
# cannot mutate the map. The bracket follows from E(Tref)=0 and E' >= Cmin.
# Finite floating-point inputs must still have representable intermediates;
# "all real states" is a mathematical contract, not immunity to overflow.
function inverse_report(H, m, model; maxiters = 160)
    isfinite(H) && all(isfinite, m) || throw(ArgumentError("nonfinite inventory"))
    kappa = regularization(m, model)
    Cmin =
        sum(positive_mass(m[i], model.delta) * l.a for (i, l) in enumerate(model.laws)) +
        kappa
    Tr = model.Tref
    H == 0 && return (; T = Tr, iterations = 0, residual = zero(H))
    endpoint = Tr + H / Cmin
    isfinite(endpoint) || throw(OverflowError("caloric bracket is not representable"))
    lo, hi = minmax(Tr, endpoint)
    T = clamp(Tr + H / extended_capacity(Tr, m, model), lo, hi)
    for iterations = 1:maxiters
        r = extended_energy(T, m, model) - H
        C = extended_capacity(T, m, model)
        isfinite(r) && isfinite(C) || throw(OverflowError("caloric evaluation overflow"))
        # Stop on temperature accuracy, including the attainable spacing at Tref.
        abs(r) <= 16eps(Float64) * max(abs(H), C * max(abs(T), abs(Tr), 1.0)) &&
            return (; T, iterations, residual = r)
        if r > 0
            hi = T
        else
            lo = T
        end
        candidate = T - r / C
        # Enforce bracket contraction: an in-bracket Newton step alone can
        # oscillate across a narrow heat-capacity peak almost indefinitely.
        margin = 0.05 * (hi - lo)
        T = lo + margin < candidate < hi - margin ? candidate : lo + (hi - lo) / 2
    end
    error(
        "nonlinear caloric inverse did not converge after $maxiters iterations: H=$H m=$m T=$T bracket=($lo,$hi)",
    )
end

# Differentiate the defining equation, not the safeguarded root search. Starting
# at the primal root, one Newton correction recovers its first-order Dual jet;
# a second recovers the second order. One correction per nesting level suffices.
# This avoids tag confusion from combining partially stripped nested Duals.
# Coefficients are fixed: property-parameter AD remains a separate gate.
# Regularization is C1; Hessians at its joints are branch derivatives, not a
# claim of global C2 smoothness.
primal(x) = x
primal(x::ForwardDiff.Dual) = primal(ForwardDiff.value(x))
depth(x) = 0
depth(x::ForwardDiff.Dual) = 1 + depth(ForwardDiff.value(x))
function invert_temperature(H, m, model)
    if H isa ForwardDiff.Dual || any(x -> x isa ForwardDiff.Dual, m)
        T = inverse_report(primal(H), primal.(m), model).T + zero(H + sum(m))
        for _ = 1:max(depth(H), maximum(depth, m))
            T += (H - extended_energy(T, m, model)) / extended_capacity(T, m, model)
        end
        return T
    end
    return inverse_report(H, m, model).T
end

end
