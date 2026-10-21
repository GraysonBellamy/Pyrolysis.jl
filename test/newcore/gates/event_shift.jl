# The event-time shift of a solve against its reference, shared by the gate
# scripts that time deactivation events (`case_accuracy.jl`,
# `burnout_timing.jl`). An invocation is one firing of the deactivation
# callback: on the finite-volume core the events that share an invocation
# number, at the time of its first event; on the old core one park time.

using Statistics: median

# The time of each deactivation invocation [s].
function invocation_times(sol)
    hasproperty(sol, :events) && return [
        first(e.t for e in sol.events if e.invocation == k) for
        k in unique(e.invocation for e in sol.events)
    ]
    return unique(sol.extras.park_times)
end

# The largest shift of the invocation times [s] and the median interval
# between the reference's invocations, or, with a single invocation, between
# its saved times [s] (NaN without events or with different counts).
function event_shift(sol, ref)
    a, b = invocation_times(sol), invocation_times(ref)
    (isempty(b) || length(a) != length(b)) && return NaN, NaN
    return maximum(abs, a .- b), median(diff(length(b) > 1 ? b : ref.t))
end
