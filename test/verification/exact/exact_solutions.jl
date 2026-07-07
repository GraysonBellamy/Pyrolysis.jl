# =============================================================================
# Exact-solution library for the verification suite
# =============================================================================
#
# Reference implementations of the closed-form solutions the verification
# cases compare against. Conventions (see docs/Verification_Suite_Plan):
#
# - Every formula is transcribed from a cited source and carries a permanent
#   self-test in `selftests.jl` (independent check: algebraic identity, limit
#   case, quadrature, or a published value). Never trust a transcription that
#   has no self-test.
# - Series are summed to machine tolerance, transcendental roots bracketed
#   and bisected to ~1e-14 — no hard-coded floats without provenance.
# - `x` always denotes depth below the exposed (top, z = L) surface, i.e.
#   x = L − z in mesh coordinates.

include("roots.jl")
include("semi_infinite.jl")
include("slab_series.jl")
