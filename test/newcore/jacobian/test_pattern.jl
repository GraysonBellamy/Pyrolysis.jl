using Test
using SparseArrays
using Pyrolysis
using Pyrolysis.Materials: DepletionLimiter, BEER_LAMBERT, SURFACE_ABSORPTION
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.States: StateLayout, state_length
using Pyrolysis.Kernels: RelativeFluxKernel, find_kernel
using Pyrolysis.Closures: EquidistributionCorrection, PrescribedCorrection
using Pyrolysis.ColoredJacobian
using Pyrolysis.Integration: setup

# ============================================================================
# The per-entry comparison
# ============================================================================

@testset "the comparison is per entry, with a row-magnitude floor" begin
    # A stiff row (1e6) above a small one (1e-2), as an Arrhenius mass row
    # above a velocity row.
    J_ref = [
        1.0e6 2.0e3
        1.0e-2 3.0e-2
    ]
    # The small row wholly missing: a global criterion accepts it
    # (3e-2 / 1e6 = 3e-8); the per-entry one must not.
    J_missing = copy(J_ref)
    J_missing[2, :] .= 0.0
    @assert maximum(abs, J_missing .- J_ref) / maximum(abs, J_ref) <= 1e-6
    cmp = compare_entries(J_missing, J_ref; rtol = 1e-6, atol = 1e-10)
    @test !cmp.passed
    @test cmp.max_abs_dev ≈ 3.0e-2
    @test cmp.worst_index[1] == 2
    # Round-off at a true zero of the stiff row passes through the row floor.
    J_round = copy(J_ref)
    J_round[1, 2] += 1.0e-9
    @test compare_entries(J_round, J_ref; rtol = 1e-6, atol = 1e-10).passed
    # A disagreement inside the small row fails, invisible at the global scale.
    J_wrong = copy(J_ref)
    J_wrong[2, 2] *= 2
    @test !compare_entries(J_wrong, J_ref; rtol = 1e-6, atol = 1e-10).passed
    cmp = compare_entries(copy(J_ref), J_ref; rtol = 1e-6, atol = 1e-10)
    @test cmp.passed
    @test cmp.max_abs_dev == 0.0
end

# ============================================================================
# Color counts
# ============================================================================

function _material(NC; law = nothing, reacting = true)
    a = SolidComponent(:a, ρ = 1000.0, c = 1000.0, k = 0.5, ε = 0.9, α = 2.0)
    b = SolidComponent(:b, ρ = 400.0, c = 1200.0, k = 0.2, ε = 0.95, α = 20.0)
    g = GaseousComponent(:g, M = 0.03, c = 1100.0, k = 0.03, λ = 1e-5)
    components = NC == 2 ? (a, g) : (a, b, g)
    reactions =
        !reacting ? () :
        NC == 2 ? (Reaction(:r, 1 => 2; A = 5e9, E = 1.4e5, h = 8e5),) :
        (Reaction(:r, 1 => (2, 3); A = 5e9, E = 1.4e5, h = 8e5, yields = (0.4, 0.6)),)
    return Material(; name = :colors, components, reactions, lateral_shrinkage_law = law)
end

function _problem(NC, n; law = nothing, reacting = true, coordinates = Planar())
    mat = _material(NC; law, reacting)
    L = 0.01
    return PyrolysisProblem(
        mesh = create_stretched_mesh(L, n, 1.05; coordinates),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = RadiativeFluxBC(flux = 3e4) + RadiativeBC(T_inf = 300.0),
                bottom = AdiabaticBC(),
            ),
            mass = (top = DiffusiveMassBC(0.02), bottom = ImpermeableBC()),
        ),
        T_initial = z -> 500.0 + 2e4 * z,
        tspan = (0.0, 10.0),
    )
end

# The new core's color count against the design's formula (§3.7),
# `3NC + 9 + 2·bl + 7·χ` with `c = 0` and `5NC + 13 + 2·bl + 11·χ` with the
# reconstruction of the relative fluxes. That model is a moving mesh whose
# cell rows read their two node velocities; at `c = 0` the cell rows of this
# core read none (the velocities enter the position rows only), two colors
# fewer, as under the equidistribution correction, whose node corrections
# read positions and masses only and whose balance rows reach no further
# than the cell rows; only a prescribed velocity, the correction itself,
# reaches them. A fixed mesh has no node block (four positions and the two
# velocities), six fewer.
_phase0(NC, bl, χ, reconstruction) =
    reconstruction ? 5NC + 13 + 2bl + 11χ : 3NC + 9 + 2bl + 7χ
_expected(NC, moving, bl, χ, reconstruction, prescribed) =
    _phase0(NC, bl, χ, reconstruction) - (!moving ? 6 : prescribed ? 0 : 2)

const _color_cases = (
    fixed = (options = (;), law = nothing, reacting = true),
    fixed_beer_lambert = (
        options = (radiation_model = BEER_LAMBERT,),
        law = nothing,
        reacting = true,
    ),
    moving = (options = (use_ale = true,), law = nothing, reacting = true),
    moving_beer_lambert = (
        options = (use_ale = true, radiation_model = BEER_LAMBERT),
        law = nothing,
        reacting = true,
    ),
    cross_section = (options = (use_ale = true,), law = χ -> 1 - 0.3χ, reacting = true),
    cross_section_beer_lambert = (
        options = (use_ale = true, radiation_model = BEER_LAMBERT),
        law = χ -> 1 - 0.3χ,
        reacting = true,
    ),
    reconstruction = (
        options = (use_ale = true, w_override = (x, t) -> -0.01x),
        law = nothing,
        reacting = false,
    ),
    reconstruction_cross_section = (
        options = (use_ale = true, w_override = (x, t) -> -0.01x),
        law = χ -> 1 - 0.3χ,
        reacting = false,
    ),
    equidistribution = (
        options = (use_ale = true, correction = EquidistributionCorrection()),
        law = nothing,
        reacting = true,
    ),
    equidistribution_beer_lambert = (
        options = (
            use_ale = true,
            radiation_model = BEER_LAMBERT,
            correction = EquidistributionCorrection(),
        ),
        law = nothing,
        reacting = true,
    ),
    equidistribution_monitor_cross_section = (
        options = (
            use_ale = true,
            radiation_model = BEER_LAMBERT,
            correction = EquidistributionCorrection(monitor = (1.0, 1.0)),
        ),
        law = χ -> 1 - 0.3χ,
        reacting = true,
    ),
)

# The largest number of stored entries in a row: the lower bound of any
# column coloring.
_largest_row(S) = maximum(diff(sparse(transpose(S)).colptr))

@testset "colors: $(name), $(NC) components" for name in keys(_color_cases), NC in (2, 3)
    case = _color_cases[name]
    counts = Int[]
    for n in (6, 8, 12, 25)
        cp = setup(_problem(NC, n; case.law, case.reacting); case.options...)
        cj = cp.jacobian
        layout = cp.def.layout
        reconstruction = find_kernel(RelativeFluxKernel, cp.def.kernels) !== nothing
        prescribed =
            cp.def.closure !== nothing && cp.def.closure.correction isa PrescribedCorrection
        @test all(cj.pattern[r, r] for r in axes(cj.pattern, 1))
        # The coloring is optimal: it reaches the largest row of the colored
        # pattern.
        @test cj.n_colors == _largest_row(cj.colored)
        @test cj.n_colors == _expected(
            NC,
            is_moving(layout),
            has_intensities(layout),
            has_chi(layout),
            reconstruction,
            prescribed,
        )
        push!(counts, cj.n_colors)
    end
    # Constant in the cell count.
    @test all(==(counts[1]), counts)
end

# ============================================================================
# Radial meshes
# ============================================================================

# The half-bandwidth of the pattern without the rows and columns of the
# global scalars, which are dense by construction.
function _bandwidth(pattern, layout)
    globals = has_chi(layout) ? (area_index(layout), lateral_rate_index(layout)) : ()
    w = 0
    rows = rowvals(pattern)
    for c in axes(pattern, 2), k in nzrange(pattern, c)
        r = rows[k]
        (r in globals || c in globals) && continue
        w = max(w, abs(r - c))
    end
    return w
end

# In-depth radiation is planar only: setup refuses it on a cylinder or a sphere.
@testset "radial meshes: the planar pattern, bandwidth and colors" begin
    for options in ((;), (use_ale = true,))
        widths = Int[]
        for n in (6, 8, 12, 25)
            cps = map(
                c -> setup(_problem(3, n; coordinates = c); options...),
                (Planar(), Cylindrical(), Spherical()),
            )
            patterns = map(cp -> cp.jacobian.pattern, cps)
            @test patterns[1] == patterns[2] == patterns[3]
            colors = map(cp -> cp.jacobian.n_colors, cps)
            @test colors[1] == colors[2] == colors[3]
            w = map(cp -> _bandwidth(cp.jacobian.pattern, cp.def.layout), cps)
            @test w[1] == w[2] == w[3]
            push!(widths, w[1])
        end
        @test all(==(widths[1]), widths)
    end
end
