# Aqua.jl package-hygiene checks.
#
# Run standalone:  julia --project=test/quality test/quality/aqua.jl
#
# Policy: a high-signal subset gates the build (the script exits non-zero on a
# regression). Gating: deps/weakdeps compat, undefined exports, method
# ambiguities, and unbound type parameters (regression gate — see ALLOWED_UNBOUND).
using Aqua
using Pyrolysis
using Test

# Methods whose only unbound type parameter is the element type of an `NTuple{N,T}`
# argument — flagged solely because the empty tuple `NTuple{0,T}` leaves T free
# (impossible here: `Material{0}` / a 0-length table never occurs). T is genuinely
# needed in each body (`zero(T)`, `new{N,T}`, type-uniform AD output tuples), so
# dropping it would *hurt* inference. We accept this documented set and gate that
# no NEW unbound-arg method appears. Identified by function-name symbol.
const ALLOWED_UNBOUND = Set([
    :_gt_face_fluxes,
    :_gt_S_conv_unrolled,
    :_gt_species_div_unrolled!,
    :_dt_face_fluxes,
    :_dt_S_conv_unrolled,
    :_dt_species_div_unrolled!,
    :_bc_mass_contribution!,
    :_bc_boundary_face_fluxes,
    :_bc_thermal_contribution,
    :_bc_convection_contribution,
    :_bcc_interior_face_fluxes,
    :apply_cell!,
    :TableProperty,
    :linear_interpolate,
])

# ---------------------------------------------------------------------------------
# Gating checks. A failure here fails CI (throws on testset finish).
# ---------------------------------------------------------------------------------
@testset "Aqua — gating checks" begin
    # Dependencies declared but unused, or used but undeclared.
    # RecursiveFactorization is intentionally a dep with no direct `import`: it is
    # loaded so LinearSolve's RFLUFactorization backend is available at runtime.
    Aqua.test_stale_deps(Pyrolysis; ignore = [:RecursiveFactorization])
    # Every dependency (incl. stdlibs and weakdeps) and `julia` has a [compat] entry.
    Aqua.test_deps_compat(Pyrolysis)
    # [extras]/test setup is consistent.
    Aqua.test_project_extras(Pyrolysis)
    # Every exported name resolves to a binding (Internal/Experimental included).
    Aqua.test_undefined_exports(Pyrolysis)
    # No method ambiguities (the CombinedThermalBC/CombinedMassBC `+` pairs are
    # disambiguated in src/BoundaryConditions/boundary_conditions.jl).
    Aqua.test_ambiguities(Pyrolysis; recursive = true)
    # No unbound type parameters beyond the documented NTuple-element-type set.
    let unbound = Aqua.detect_unbound_args_recursively(Pyrolysis)
        unexpected = filter(m -> !(m.name in ALLOWED_UNBOUND), unbound)
        if !isempty(unexpected)
            @warn "Aqua: $(length(unexpected)) NEW method(s) with unbound type parameters"
            for m in unexpected
                @warn "  unbound args: $m"
            end
        end
        @test isempty(unexpected)
    end
end

println("Aqua gating checks passed.")
