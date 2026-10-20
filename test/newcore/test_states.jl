using Test
using ForwardDiff
using LinearAlgebra
using Pyrolysis
using Pyrolysis.Materials: T_REF, T_PROPERTY_MIN, T_PROPERTY_MAX, BEER_LAMBERT
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.Meshes: MeshGeometry
using Pyrolysis.States: StateLayout, state_length
using Pyrolysis.Integration: setup

@testset "state layout: interleaved blocks and row kinds" begin
    for (moving, intensities, chi) in (
        (false, false, false),
        (true, false, false),
        (false, true, false),
        (true, true, false),
        (true, true, true),
        (true, false, true),
    )
        NC, n = 3, 5
        l = StateLayout(NC, n, TemperatureState(); moving, intensities, chi)
        rows = Int[]
        for i = 1:n
            append!(rows, (mass_index(l, i, j) for j = 1:NC))
            push!(rows, energy_index(l, i))
            chi && append!(rows, (history_index(l, i, 1), history_index(l, i, 2)))
        end
        for nd = 1:(n+1)
            moving && append!(rows, (position_index(l, nd), velocity_index(l, nd)))
            intensities && push!(rows, intensity_index(l, nd))
        end
        chi && append!(rows, (area_index(l), lateral_rate_index(l)))
        @test sort(rows) == collect(1:state_length(l))
        @test state_length(l) ==
              n * (NC + 1 + 2chi) + (n + 1) * (2moving + intensities) + 2chi
        # Interleaving: node n's block, then cell n's block.
        for i = 1:n
            @test first(cell_rows(l, i)) == last(node_rows(l, i)) + 1
            @test first(node_rows(l, i + 1)) == last(cell_rows(l, i)) + 1
        end
        d = is_differential(l)
        @test count(!, d) == (n + 1) * (moving + intensities) + 2chi
        moving && @test !d[velocity_index(l, 2)] && d[position_index(l, 2)]
        intensities && @test !d[intensity_index(l, 1)]
        @test mass_matrix(l) == Diagonal(Float64.(d))
        @test all(l.mask)
    end
    @test_throws ArgumentError StateLayout(2, 3, TemperatureState(); chi = true)
    l = StateLayout(2, 3, EnthalpyState())
    @test energy_state(l) isa EnthalpyState
    @test_throws ArgumentError position_index(l, 1)
end

@testset "smooth guards: C¹, exact on physical states, finite everywhere" begin
    for T in (-1e6, -1.0, 0.0, 199.0, 200.0, 200.5, 300.0, 2999.5, 3000.0, 1e6)
        Tt = trial_temperature(T)
        @test T_PROPERTY_MIN <= Tt <= T_PROPERTY_MAX
        kernels(y) = (
            t = trial_temperature(y[1]);
            [exp(-1e5 / (8.314 * t)), sqrt(t), log(t), t^1.75, t^4]
        )
        @test all(isfinite, kernels([T]))
        @test all(isfinite, ForwardDiff.jacobian(kernels, [T]))
    end
    @test trial_temperature(700.0) === 700.0
    @test ForwardDiff.derivative(trial_temperature, 700.0) == 1
    for edge in (T_PROPERTY_MIN, T_PROPERTY_MAX)
        @test ForwardDiff.derivative(trial_temperature, edge) == 0
    end
    for x in (T_PROPERTY_MIN + 1.0, T_PROPERTY_MAX - 1.0)
        @test ForwardDiff.derivative(trial_temperature, x) ≈ 1
    end
    # Negative control: the raw Arrhenius factor has a NaN derivative at zero.
    @test !isfinite(ForwardDiff.derivative(T -> exp(-1000 / T), 0.0))

    δ = 1e-8
    for m in (-1.0, 0.0, δ / 10, δ / 2, δ, 10δ, 1.0)
        p = positive_part(m, δ)
        @test p >= 0
        @test isfinite(ForwardDiff.derivative(x -> positive_part(x, δ), m))
        m <= 0 && @test p == 0
        m >= δ && @test p == m
        0 < m < δ && @test abs(p - m) <= δ
    end
    @test ForwardDiff.derivative(x -> positive_part(x, δ), 0.0) == 0
    @test ForwardDiff.derivative(x -> positive_part(x, δ), δ) == 1
    for x in (-1.0, 0.5, 1.0, 1.5, 2.0, 3.0)
        @test smooth_floor(x, 1.0) >= 1.0
        x >= 2 && @test smooth_floor(x, 1.0) == x
    end
    @test ForwardDiff.derivative(x -> smooth_floor(x, 1.0), 1.0) == 0
    @test ForwardDiff.derivative(x -> smooth_floor(x, 1.0), 2.0) == 1
end

# Heat-capacity laws of every kind the enthalpy must integrate.
const _cp_laws = (
    ConstantProperty(1500.0),
    LinearProperty(900.0, 1.2),
    PolynomialProperty((700.0, 0.9, -2e-4, 3e-8)),
    TableProperty((150.0, 400.0, 900.0, 3500.0), (800.0, 1200.0, 1900.0, 2400.0)),
    CallableProperty(T -> 800.0 + 6000.0 / (1 + ((T - 800.0) / 12.0)^2)),
)

# Independent integral: composite Simpson's rule, fine enough to resolve the
# peaked law, with each law evaluated through its own (clamped) callable.
function _simpson(c, a, b; n = 400_000)
    a == b && return 0.0
    h = (b - a) / n
    s = c(a) + c(b)
    for k = 1:(n-1)
        s += (isodd(k) ? 4 : 2) * c(a + k * h)
    end
    return s * h / 3
end

@testset "sensible enthalpy of every heat-capacity law" begin
    temps =
        (-500.0, 150.0, 200.0, 250.0, T_REF, 401.0, 799.0, 805.0, 2999.0, 3000.0, 4000.0)
    # Closed-form laws: the slope is the law to round-off.
    for c in _cp_laws[1:4]
        h = enthalpy_law(c)
        @test h === c
        @test sensible_enthalpy(h, T_REF) == 0
        for T in temps
            @test ForwardDiff.derivative(t -> sensible_enthalpy(h, t), T) ≈ c(T) rtol =
                1e-12
            @test enthalpy_slope(h, T) == c(T)
            @test sensible_enthalpy(h, T) ≈ _simpson(c, T_REF, T) rtol = 1e-10 atol = 1e-6
        end
    end
    # Tabulated laws: C¹, the slope is the interpolant's exact derivative, and
    # a quadratic law (cubic enthalpy) is reproduced exactly.
    quadratic = CallableProperty(T -> 1000.0 + 0.5T - 1e-4T^2)
    for (c, value_tol, slope_tol) in ((quadratic, 1e-12, 1e-12), (_cp_laws[5], 1e-7, 1e-4))
        h = enthalpy_law(c)
        @test h isa TabulatedEnthalpy
        @test abs(sensible_enthalpy(h, T_REF)) < 1e-9
        cmax = maximum(c(T) for T = 200.0:0.5:3000.0)
        for T in (temps..., 233.37, 812.5, 1777.77)
            @test ForwardDiff.derivative(t -> sensible_enthalpy(h, t), T) ≈
                  enthalpy_slope(h, T) rtol = 1e-12
            @test abs(enthalpy_slope(h, T) - c(T)) <= slope_tol * cmax
            @test sensible_enthalpy(h, T) ≈ _simpson(c, T_REF, T) rtol = value_tol atol =
                1e-5
        end
        # Continuous slope at a knot.
        @test enthalpy_slope(h, 800.0 - 1e-9) ≈ enthalpy_slope(h, 800.0 + 1e-9) rtol = 1e-7
        # A NaN trial temperature gives NaN, value and derivative, not an error.
        @test isnan(sensible_enthalpy(h, NaN))
        @test isnan(enthalpy_slope(h, NaN))
        @test isnan(ForwardDiff.derivative(t -> sensible_enthalpy(h, t), NaN))
        @test isnan(ForwardDiff.derivative(t -> enthalpy_slope(h, t), NaN))
    end
end

# A material whose components carry the peaked and polynomial laws.
_caloric_material() = Material(
    name = :caloric,
    components = (
        SolidComponent(:poly, ρ = 1000.0, c = PolynomialProperty((700.0, 0.9)), k = 0.3),
        SolidComponent(:peak, ρ = 800.0, c = _cp_laws[5], k = 0.2),
        GaseousComponent(:gas, M = 0.03, c = 1100.0, k = 0.03, λ = 0.0),
    ),
    reactions = (),
)

# Independent 256-bit bisection on the extended relation itself.
function _oracle(ext, laws, H, m, i)
    setprecision(256) do
        lo, hi = big"-1e9", big"1e9"
        for _ = 1:300
            mid = (lo + hi) / 2
            if extended_energy(ext, laws, mid, BigFloat.(m), i) < BigFloat(H)
                lo = mid
            else
                hi = mid
            end
        end
        Float64((lo + hi) / 2)
    end
end

@testset "extended caloric map: inversion, implicit derivatives, exact region" begin
    material = _caloric_material()
    laws = CaloricLaws(material)
    δ, floor = 1e-8, 1e-3
    ext = CaloricExtension(laws, [δ], [floor])
    masses = (
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.7, 0.3, 1e-3),
        (δ, 0.0, 0.0),
        (δ / 2, -1.0, 0.0),
        (0.0, 0.0, 0.0),
        (-1.0, -2.0, -0.5),
        (1e-18, 2e-18, 0.0),
    )
    worst = 0.0
    for m in masses, target in (-10000.0, -1.0, 200.0, T_REF, 800.0, 3000.0, 1e5)
        H = extended_energy(ext, laws, target, m, 1)
        for start in (T_REF, 750.0, 1e4)
            T, _ = invert_caloric(ext, laws, H, m, 1, start)
            ref = _oracle(ext, laws, H, m, 1)
            worst = max(worst, abs(T - ref) / max(1, abs(ref)))
            @test T ≈ ref atol = 1e-9 rtol = 4e-15
        end
        @test extended_capacity(ext, laws, target, m, 1) >= floor
        f(v) = first(invert_caloric(ext, laws, v[1], (v[2], v[3], v[4]), 1, T_REF))
        v = [H; collect(m)]
        g = ForwardDiff.gradient(f, v)
        @test all(isfinite, g)
        @test g[1] ≈ 1 / extended_capacity(ext, laws, f(v), m, 1) rtol = 1e-12
        @test all(isfinite, ForwardDiff.hessian(f, v))
        identity(v) = extended_energy(ext, laws, f(v), (v[2], v[3], v[4]), 1) - v[1]
        # The identity cancels terms of the size of the component enthalpies.
        Tv = f(v)
        scale = max(
            abs(H),
            1.0,
            extended_capacity(ext, laws, Tv, m, 1) * abs(Tv),
            maximum(abs, component_enthalpies(laws, Tv)),
        )
        @test norm(ForwardDiff.gradient(identity, v), Inf) <= 1e-14 * scale
    end
    @test worst < 1e-13

    # Exact agreement with the physical relation: every species absent or at
    # least δ, and the guaranteed capacity above twice the floor.
    for m in ((1.0, 0.0, 0.0), (0.0, δ * 1e6, 0.0), (0.3, 0.7, 1e-6)),
        T in (250.0, 812.0, 2500.0)

        H = sum(m[j] * sensible_enthalpy(laws.h[j], T) for j = 1:3)
        @test regularization(ext, m, 1) == 0
        @test first(invert_caloric(ext, laws, H, m, 1, T_REF)) ≈ T atol = 1e-10
    end
    # Inside the onset the regularization error is bounded by δ c ΔT.
    for m1 in (δ / 4, δ / 2, 3δ / 4)
        m = (m1, 1.0, 0.0)
        Hphys = sum(m[j] * sensible_enthalpy(laws.h[j], 900.0) for j = 1:3)
        @test abs(extended_energy(ext, laws, 900.0, m, 1) - Hphys) <=
              δ * 900.0 * (900.0 - T_REF)
    end
    # The capacity switch keeps the relation invertible as the cell empties.
    for scale in (0.0, 1e-9, 1e-6, 1e-5)
        m = (scale, 0.0, 0.0)
        @test extended_capacity(ext, laws, 500.0, m, 1) >= floor
        @test regularization(ext, m, 1) >= 0
    end
    # A gas that carries no heat: the relation stays invertible with a slope of
    # at least the floor, a cell of that gas alone included.
    cold = Material(;
        name = :cold_gas,
        components = (
            SolidComponent(:s, ρ = 1000.0, c = 900.0, k = 0.3),
            GaseousComponent(:g, M = 0.1, c = 0.0, k = 0.02, λ = 0.0),
        ),
        reactions = (),
    )
    cold_laws = CaloricLaws(cold)
    cold_ext = CaloricExtension(cold_laws, [δ], [floor])
    for m in ((1.0, 0.5), (0.0, 0.5), (0.0, 0.0)), target in (250.0, 900.0)
        H = extended_energy(cold_ext, cold_laws, target, m, 1)
        @test first(invert_caloric(cold_ext, cold_laws, H, m, 1, T_REF)) ≈ target rtol =
            1e-12
        @test extended_capacity(cold_ext, cold_laws, target, m, 1) >= floor
    end
    @test isnan(first(invert_caloric(ext, laws, Inf, (1.0, 0.0, 0.0), 1, T_REF)))
    @test isnan(first(invert_caloric(ext, laws, NaN, (1.0, 0.0, 0.0), 1, T_REF)))
    @test isnan(first(invert_caloric(ext, laws, 1e3, (NaN, 0.0, 0.0), 1, T_REF)))
    @test_throws ErrorException invert_caloric(
        ext,
        laws,
        1e9,
        (1.0, 0.0, 0.0),
        1,
        T_REF;
        maxiters = 1,
    )
end

@testset "state mapping: masses, volumes, temperatures, trial states" begin
    material = _caloric_material()
    mesh = line_mesh(collect(range(0.0, 0.01, length = 7)); cross_section = 2.0)
    topo = mesh.topology
    n = 6
    m0 = [j == 1 ? 1000.0 * 2.0 * 0.01 / n : 0.0 for j = 1:3, i = 1:n]
    ref = CellReference(mesh, m0, material; reltol = 1e-6, threshold = 1e-2)
    @test ref.volume ≈ fill(2.0 * 0.01 / n, n)
    @test ref.solid_mass ≈ m0[1, :]
    @test ref.mass_scale ≈ 1e-8 .* m0[1, :]
    for energy in (TemperatureState(), EnthalpyState())
        ext =
            energy isa EnthalpyState ?
            CaloricExtension(CaloricLaws(material), ref.mass_scale, ref.capacity_floor) :
            nothing
        setup = StateSetup(topo, material, ref, ext)
        layout = StateLayout(3, n, energy)
        T0 = collect(range(300.0, 900.0, length = n))
        u = zeros(state_length(layout))
        for i = 1:n
            mi = (m0[1, i] * (1 - 0.1i), m0[1, i] * 0.1i, 1e-6)
            for j = 1:3
                u[mass_index(layout, i, j)] = mi[j]
            end
            u[energy_index(layout, i)] = energy_state_value(energy, setup, T0[i], mi, i)
        end
        ds = DerivedState{Float64}(layout, setup, copy(T0))
        derive!(ds, u, layout, setup)
        @test ds.T ≈ T0 rtol = 1e-14
        @test ds.T_trial == ds.T
        @test ds.V ≈ ref.volume
        @test ds.ξ[1, :] ≈ [u[mass_index(layout, i, 1)] for i = 1:n] ./ ref.volume
        @test all(ds.k .> 0)
        # Trial states off the physical manifold: finite values and derivatives.
        for bad in (zeros(length(u)), -abs.(u), u .* (1 .+ 0.3 .* sin.(1:length(u))))
            if energy isa TemperatureState
                bad[energy_index(layout, 1)] = -50.0
                bad[energy_index(layout, 2)] = 0.0
            end
            derive!(ds, bad, layout, setup)
            @test all(isfinite, ds.T) && all(isfinite, ds.ξ) && all(isfinite, ds.cp_sum)
            @test all(>(0), ds.cp_sum)
            J = ForwardDiff.jacobian(bad) do v
                dsd = DerivedState{eltype(v)}(layout, setup, copy(T0))
                derive!(dsd, v, layout, setup)
                [dsd.T; vec(dsd.ξ⁺); dsd.cp_sum; dsd.k; vec(dsd.h)]
            end
            @test all(isfinite, J)
        end
    end
end

@testset "carried fields: axis priority and physical tensors" begin
    reference = [1 -1 0; 1 1 0; 0 0 sqrt(2)] / sqrt(2)
    F = Diagonal([2.0, 1.0, 1.0])
    axes = material_axes(Matrix(F), reference)
    @test axes ≈ [2 -1 0; 1 2 0; 0 0 sqrt(5)] / sqrt(5) atol = 1e-14
    @test axes' * axes ≈ I atol = 1e-14
    @test det(axes) ≈ 1 atol = 1e-14
    @test atan(axes[2, 1], axes[1, 1]) ≈ atan(1 / 2) atol = 1e-14
    k = physical_tensor([4.0, 1.0, 1.0], axes)
    @test k[1:2, 1:2] ≈ [3.4 1.2; 1.2 1.6] atol = 1e-14
    swapped = material_axes(Matrix(F), reference[:, [2, 1, 3]])
    @test physical_tensor([1.0, 4.0, 1.0], swapped)[1:2, 1:2] ≈ [1.6 1.2; 1.2 3.4] atol =
        1e-14
    rotation = [0.0 -1 0; 1 0 0; 0 0 1]
    rotated = material_axes(rotation, Matrix{Float64}(I, 3, 3))
    @test physical_tensor([4.0, 1.0, 1.0], rotated) ≈ Diagonal([1.0, 4.0, 1.0]) atol = 1e-14
    @test material_axes(rotation * F, reference) ≈ rotation * axes atol = 1e-14
    @test all(
        isfinite,
        ForwardDiff.jacobian(
            x -> vec(material_axes(reshape(x, 3, 3), reference)),
            vec(Matrix(F)),
        ),
    )
    @test material_axes(fill(3.0, 1, 1), fill(1.0, 1, 1)) == fill(1.0, 1, 1)
    @test [f.name for f in CONVERSION_HISTORY] == [:G, :M]
    @test_throws ArgumentError CarriedField(:x, :mesh, false, :none)
end

@testset "the caloric extension refuses a non-positive scale, floor or heat capacity" begin
    laws = CaloricLaws(_caloric_material())
    @test CaloricExtension(laws, [1e-8], [1e-3]) isa CaloricExtension
    @test_throws ArgumentError CaloricExtension(laws, [0.0], [1e-3])
    @test_throws ArgumentError CaloricExtension(laws, [1e-8], [-1e-3])
    # A heat capacity that turns negative inside the property bracket, which a
    # component already refuses: laws built directly.
    c = (LinearProperty(1000.0, -1.0),)
    h = map(Pyrolysis.States.enthalpy_law, c)
    negative = CaloricLaws{1,typeof(h),typeof(c)}(h, c)
    @test_throws ArgumentError CaloricExtension(negative, [1e-8], [1e-3])
end

# A slab that thins (the char is lighter), absorbs in depth and shrinks
# laterally, on a stretched mesh: every row kind of the tolerance vector.
function _thinning_slab(; law = χ -> 1 - 0.3χ)
    virgin =
        SolidComponent(:virgin, ρ = 1100.0, c = (900.0, 1.5), k = 0.2, ε = 0.8, α = 1.5)
    char = SolidComponent(:char, ρ = 300.0, c = 1100.0, k = 0.1, ε = 0.95, α = 20.0)
    gas = GaseousComponent(:gas, M = 0.03, c = 1100.0, k = 0.03, λ = 1e-5)
    mat = Material(
        name = :thinning,
        components = (virgin, char, gas),
        reactions = (
            Reaction(:pyro, 1 => (2, 3); A = 5e9, E = 1.4e5, h = 8e5, yields = (0.3, 0.7)),
        ),
        lateral_shrinkage_law = law,
    )
    return PyrolysisProblem(
        mesh = create_stretched_mesh(0.01, 6, Val(3), 1.2),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = RadiativeFluxBC(
                    flux = t -> 3e4 * (1 + 0.5 * sin(t)),
                    absorptivity = 0.9,
                ),
                bottom = AdiabaticBC(),
            ),
        ),
        T_initial = z -> 400.0 + 1e4 * z,
        ξ_initial = [z -> 1100.0 * (1 - 20z), z -> 300.0 * 20z, z -> 0.2],
        tspan = (0.0, 10.0),
    )
end

@testset "absolute tolerances of every row kind: $(energy) form" for energy in (
    :temperature,
    :enthalpy,
)
    reltol, abstol, threshold, dt = 1e-5, 1e-7, 2e-2, 1e-3
    cp = setup(
        _thinning_slab();
        energy_state = energy,
        use_ale = true,
        radiation_model = BEER_LAMBERT,
        reltol,
        abstol,
        threshold,
        time_scale = dt,
    )
    layout, ref = cp.def.layout, cp.def.setup.reference
    @test is_moving(layout) && has_chi(layout) && has_intensities(layout)
    tol = cp.abstol
    n = n_cells(layout)
    x = ref.positions
    m0 = cp.m0
    # Cell rows: the masses and the history pair at the mass tolerance of the
    # cell's initial condensed mass; the energy row at the temperature
    # tolerance or, in the enthalpy form, at the enthalpy of that tolerance
    # over the threshold's share of the condensed mass.
    for i = 1:n
        condensed = m0[1, i] + m0[2, i]
        @test ref.solid_mass[i] ≈ condensed rtol = 1e-14
        mass_tol = reltol * threshold * condensed
        for j = 1:3
            @test tol[mass_index(layout, i, j)] ≈ mass_tol rtol = 1e-14
        end
        for k = 1:2
            @test tol[history_index(layout, i, k)] ≈ mass_tol rtol = 1e-14
        end
        T = cp.T0[i]
        c_ref = (m0[1, i] * (900.0 + 1.5T) + (m0[2, i] + m0[3, i]) * 1100.0) / sum(m0[:, i])
        expected = energy === :temperature ? abstol : threshold * condensed * c_ref * abstol
        @test tol[energy_index(layout, i)] ≈ expected rtol = 1e-12
    end
    # Node rows: positions at the mass tolerance's share of the thinner
    # neighbouring cell, velocities that over the first step, the face
    # intensities at the relative tolerance of the largest incident intensity.
    I_scale = maximum(t -> 0.9 * 3e4 * (1 + 0.5 * sin(t)), range(0.0, 10.0; length = 101))
    w_tol = Float64[]
    for nd = 1:(n+1)
        h = min(nd > 1 ? x[nd] - x[nd-1] : Inf, nd <= n ? x[nd+1] - x[nd] : Inf)
        @test tol[position_index(layout, nd)] ≈ reltol * threshold * h rtol = 1e-14
        @test tol[velocity_index(layout, nd)] ≈ reltol * threshold * h / dt rtol = 1e-14
        @test tol[intensity_index(layout, nd)] ≈ reltol * I_scale rtol = 1e-14
        push!(w_tol, reltol * threshold * h / dt)
    end
    # On the stretched mesh a node's two cells differ, so the rule picks one.
    @test maximum(diff(x)) > 1.5 * minimum(diff(x))
    # The cross section at the relative tolerance of its initial measure, its
    # lateral rate at the smallest velocity tolerance over the thickness.
    @test tol[area_index(layout)] ≈ reltol * ref.cross_section rtol = 1e-14
    @test tol[lateral_rate_index(layout)] ≈ minimum(w_tol) / (x[end] - x[1]) rtol = 1e-14
    @test length(tol) == state_length(layout)
end

@testset "derived volumes on a moving mesh are the shells between the nodes: $(nameof(typeof(coords)))" for coords in
                                                                                                            (
    Planar(),
    Cylindrical(),
    Spherical(),
)
    cp = setup(_thinning_slab(law = nothing); use_ale = true, coordinates = coords)
    layout, s = cp.def.layout, cp.def.setup
    topo = s.topology
    @test topo.coordinates isa typeof(coords)
    n = n_cells(layout)
    u = copy(cp.u0)
    # Move the nodes off their initial places, keeping their order.
    for nd = 2:(n+1)
        u[position_index(layout, nd)] *= 1 + 0.05 * sin(nd)
    end
    x = [u[position_index(layout, nd)] for nd = 1:(n+1)]
    @test issorted(x) && x != s.reference.positions
    ds = cp.ws.derived
    derive!(ds, u, layout, s)
    @test ds.x == x
    for i = 1:n
        @test ds.V[i] ≈ shell_volume(coords, x[i], x[i+1], topo.cross_section) rtol = 1e-14
    end
end
