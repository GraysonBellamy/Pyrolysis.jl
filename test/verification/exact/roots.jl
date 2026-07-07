# =============================================================================
# Root finding for transcendental eigenconditions
# =============================================================================

"""
    bisect_root(f, a, b; iters = 200)

Bisection to ~1e-15 relative bracket width. `f(a)` and `f(b)` must have
opposite signs. Deterministic and immune to the derivative pathologies of the
eigencondition functions (tan poles), which is why the verification suite
uses bisection over Newton for reference roots.
"""
function bisect_root(f, a::Float64, b::Float64; iters::Int = 200)
    fa = f(a)
    fb = f(b)
    fa * fb <= 0 || throw(ArgumentError("root not bracketed on [$a, $b]"))
    for _ = 1:iters
        m = 0.5 * (a + b)
        (m == a || m == b) && break  # bracket at floating-point resolution
        fm = f(m)
        if fa * fm <= 0
            b = m
        else
            a, fa = m, fm
        end
    end
    return 0.5 * (a + b)
end
