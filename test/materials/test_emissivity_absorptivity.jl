"""
Per-component temperature-dependent emissivity and per-component surface
absorptivity.

- `ε` accepts the same forms as `k` (number, `(p0, p1)` linear tuple,
  `Function`, property function) and is evaluated at the cell temperature in
  `ε_eff = Σⱼ vⱼ εⱼ(T)`.
- `absorptivity` is an optional constant per condensed component. The
  mixture `a_eff = Σⱼ vⱼ aⱼ` uses `aⱼ = εⱼ(T)` (Kirchhoff) for components
  that leave it unset, so with none set `a_eff == ε_eff` bit-for-bit — the
  invariant for materials that set no absorptivity.
- The derived state of the finite-volume core carries both per cell
  (`ds.ε`, `ds.a`), the unclamped mixtures `mixture_emissivity` and
  `mixture_absorptivity` at the cell's composition and temperature.
- ForwardDiff through the mixtures yields the component-emissivity slope.
"""

using Test
using ForwardDiff
using Pyrolysis
using Pyrolysis.Materials:
    get_emissivity,
    get_absorptivity,
    surface_absorptivity,
    ConstantProperty,
    LinearProperty,
    CallableProperty,
    AbstractPropertyFunction
using Pyrolysis.Properties:
    mixture_emissivity, mixture_absorptivity, mixing_volume_fractions

const _EA_GAS = GaseousComponent(:g, M = 0.03, c = 1000.0, k = 0.03, λ = 1e-5)

# Two condensed components with constant bulk densities, so the mixing
# weights v_j = γ_j ξ_j / ρ_j (normalised over the matrix) are independent of
# T and the only temperature dependence is the component emissivity.
_ε_lin(T) = 0.95 - 1e-4 * T
_ε_fun(T) = 0.80 + 5e-5 * (T - 300.0)

function _ea_material(; a_virgin = nothing, a_char = nothing, ε_virgin = (0.95, -1e-4))
    # `ε_virgin = 1.3` below probes the ε_eff clamp with an unphysical law, so
    # the constructor's own range check is bypassed for this fixture.
    virgin = SolidComponent(
        :virgin,
        ρ = 500.0,
        c = 1500.0,
        k = 0.15,
        ε = ε_virgin,
        absorptivity = a_virgin,
        validate_properties = false,
    )
    char = SolidComponent(
        :char,
        ρ = 150.0,
        c = 1700.0,
        k = 0.08,
        ε = _ε_fun,
        absorptivity = a_char,
    )
    return Material(name = :ea, components = (virgin, char, _EA_GAS), reactions = ())
end

# Documented weights: v_j = γ_j ξ_j / ρ_j normalised over the contributing
# (γ > 0) components; the gas has γ = 0.
function _hand_weights(ξ)
    w = (ξ[1] / 500.0, ξ[2] / 150.0)
    s = w[1] + w[2]
    return (w[1] / s, w[2] / s)
end

@testset "component ε accepts every property-function form" begin
    s_const = SolidComponent(:s, ρ = 1.0, c = 1.0, k = 1.0, ε = 0.9)
    s_lin = SolidComponent(:s, ρ = 1.0, c = 1.0, k = 1.0, ε = (0.95, -1e-4))
    s_fun = SolidComponent(:s, ρ = 1.0, c = 1.0, k = 1.0, ε = _ε_fun)
    s_pf = SolidComponent(:s, ρ = 1.0, c = 1.0, k = 1.0, ε = LinearProperty(0.95, -1e-4))
    @test s_const.emissivity isa ConstantProperty
    @test s_lin.emissivity isa LinearProperty
    @test s_fun.emissivity isa CallableProperty
    @test s_pf.emissivity isa LinearProperty
    @test get_emissivity(s_const, 300.0) == 0.9
    @test get_emissivity(s_const, 900.0) == 0.9
    @test get_emissivity(s_lin, 500.0) ≈ _ε_lin(500.0)
    @test get_emissivity(s_fun, 700.0) ≈ _ε_fun(700.0)
    @test get_emissivity(s_pf, 500.0) ≈ _ε_lin(500.0)

    # Defaults are unchanged and are property functions.
    @test get_emissivity(SolidComponent(:s, ρ = 1.0, c = 1.0, k = 1.0), 400.0) == 0.9
    @test get_emissivity(LiquidComponent(:l, ρ = 1.0, c = 1.0, k = 1.0), 400.0) == 0.95
    @test get_emissivity(_EA_GAS, 400.0) == 0.0
    @test _EA_GAS.emissivity isa AbstractPropertyFunction
    g_T = GaseousComponent(:g, M = 0.03, c = 1.0, k = 1.0, λ = 1e-5, ε = (0.1, 1e-4))
    @test get_emissivity(g_T, 500.0) ≈ 0.15

    # Duals propagate through the component evaluation.
    @test ForwardDiff.derivative(T -> get_emissivity(s_lin, T), 500.0) ≈ -1e-4
    @test ForwardDiff.derivative(T -> get_emissivity(s_const, T), 500.0) == 0.0
end

@testset "surface absorptivity: construction, validation, accessors" begin
    s_none = SolidComponent(:s, ρ = 1.0, c = 1.0, k = 1.0, ε = (0.95, -1e-4))
    s_set =
        SolidComponent(:s, ρ = 1.0, c = 1.0, k = 1.0, ε = (0.95, -1e-4), absorptivity = 0.7)
    l_set = LiquidComponent(:l, ρ = 1.0, c = 1.0, k = 1.0, absorptivity = 1)
    @test surface_absorptivity(s_none) === nothing
    @test surface_absorptivity(s_set) === 0.7
    @test surface_absorptivity(l_set) === 1.0
    @test surface_absorptivity(_EA_GAS) === nothing
    # Nothing / Float64 are struct type parameters (static Kirchhoff branch).
    @test s_none.surface_absorptivity === nothing
    @test typeof(s_set).parameters[7] === Float64
    @test typeof(s_none).parameters[7] === Nothing

    # Kirchhoff fallback tracks ε(T); an explicit value is T-independent.
    @test get_absorptivity(s_none, 500.0) == get_emissivity(s_none, 500.0)
    @test get_absorptivity(s_none, 900.0) == get_emissivity(s_none, 900.0)
    @test get_absorptivity(s_set, 500.0) == 0.7
    @test get_absorptivity(s_set, 900.0) == 0.7
    @test get_absorptivity(_EA_GAS, 500.0) == 0.0
    @test ForwardDiff.derivative(T -> get_absorptivity(s_set, T), 500.0) == 0.0
    @test ForwardDiff.derivative(T -> get_absorptivity(s_none, T), 500.0) ≈ -1e-4

    @test_throws ArgumentError SolidComponent(
        :bad,
        ρ = 1.0,
        c = 1.0,
        k = 1.0,
        absorptivity = 1.2,
    )
    @test_throws ArgumentError LiquidComponent(
        :bad,
        ρ = 1.0,
        c = 1.0,
        k = 1.0,
        absorptivity = -0.1,
    )
end

@testset "ε_eff = Σ v_j ε_j(T) at the cell temperature" begin
    mat = _ea_material()
    ξ = (300.0, 45.0, 0.5)
    v = _hand_weights(ξ)
    @test collect(mixing_volume_fractions(ξ, 500.0, mat)[1:2]) ≈ collect(v)
    @test mixing_volume_fractions(ξ, 500.0, mat)[3] == 0.0

    for T in (500.0, 1000.0)
        ε_hand = v[1] * _ε_lin(T) + v[2] * _ε_fun(T)
        @test mixture_emissivity(ξ, T, mat) ≈ ε_hand rtol = 1e-14
    end
    @test mixture_emissivity(ξ, 500.0, mat) != mixture_emissivity(ξ, 1000.0, mat)

    # Scalar ε on every component: the constant volume-weighted sum.
    mat_c = Material(
        name = :ea_const,
        components = (
            SolidComponent(:virgin, ρ = 500.0, c = 1500.0, k = 0.15, ε = 0.9),
            SolidComponent(:char, ρ = 150.0, c = 1700.0, k = 0.08, ε = 0.85),
            _EA_GAS,
        ),
        reactions = (),
    )
    ε_old = v[1] * 0.9 + v[2] * 0.85
    @test mixture_emissivity(ξ, 500.0, mat_c) ≈ ε_old rtol = 1e-14
    @test mixture_emissivity(ξ, 1000.0, mat_c) == mixture_emissivity(ξ, 500.0, mat_c)

    # AD: dε_eff/dT is the weighted component slope (v is T-independent here).
    dε_mix = ForwardDiff.derivative(T -> mixture_emissivity(ξ, T, mat), 700.0)
    @test dε_mix ≈ v[1] * (-1e-4) + v[2] * 5e-5 rtol = 1e-12
end

@testset "a_eff: Kirchhoff invariant and per-component override" begin
    ξ = (300.0, 45.0, 0.5)
    v = _hand_weights(ξ)

    # No component sets an absorptivity → a_eff is ε_eff, bit-for-bit, for
    # T-dependent and constant emissivities alike.
    for mat in (_ea_material(), _ea_material(ε_virgin = 0.9)), T in (400.0, 800.0)
        @test mixture_absorptivity(ξ, T, mat) === mixture_emissivity(ξ, T, mat)
    end

    # One component overrides: its constant replaces ε_j(T) in the same sum.
    mat_a = _ea_material(a_char = 0.6)
    for T in (500.0, 1000.0)
        a_hand = v[1] * _ε_lin(T) + v[2] * 0.6
        @test mixture_absorptivity(ξ, T, mat_a) ≈ a_hand rtol = 1e-14
        @test mixture_absorptivity(ξ, T, mat_a) != mixture_emissivity(ξ, T, mat_a)
        # Emission is untouched by the override.
        @test mixture_emissivity(ξ, T, mat_a) == mixture_emissivity(ξ, T, _ea_material())
    end
    # Only the Kirchhoff component contributes a temperature slope.
    da = ForwardDiff.derivative(T -> mixture_absorptivity(ξ, T, mat_a), 700.0)
    @test da ≈ v[1] * (-1e-4) rtol = 1e-12

    # Both override: a_eff is T-independent.
    mat_ab = _ea_material(a_virgin = 0.75, a_char = 0.6)
    @test mixture_absorptivity(ξ, 500.0, mat_ab) ≈ v[1] * 0.75 + v[2] * 0.6 rtol = 1e-14
    @test mixture_absorptivity(ξ, 500.0, mat_ab) == mixture_absorptivity(ξ, 1000.0, mat_ab)

    # Depleted cell (gas only): both radiative mixtures are 0.
    @test mixture_emissivity((0.0, 0.0, 1.0), 500.0, mat_a) == 0.0
    @test mixture_absorptivity((0.0, 0.0, 1.0), 500.0, mat_a) == 0.0
end

# The per-cell values the residual reads: one evaluation at a graded state.
@testset "the derived state fills ε and a per cell at the cell's composition and temperature" begin
    _composition(i) = (300.0 - 20.0 * i, 10.0 * i, 0.5)
    _temperature(i) = 400.0 + 100.0 * i
    function _derived_for(mat)
        problem = PyrolysisProblem(
            mesh = create_uniform_mesh(0.01, 5),
            material = mat,
            bc_set = BoundaryConditionSet(
                thermal = (top = HeatFluxBC(0.0), bottom = AdiabaticBC()),
            ),
            tspan = (0.0, 1.0),
        )
        cp = Pyrolysis.Integration.setup(problem)
        layout = cp.def.layout
        V = Pyrolysis.Problem.cell_volumes(problem.mesh)
        u = copy(cp.u0)
        for i = 1:5
            u[Pyrolysis.States.energy_index(layout, i)] = _temperature(i)
            for (j, ξ) in enumerate(_composition(i))
                u[Pyrolysis.States.mass_index(layout, i, j)] = ξ * V[i]
            end
        end
        du = zeros(length(u))
        Pyrolysis.Assembly.residual!(du, u, cp.def, 0.0, cp.ws)
        ds, _ = Pyrolysis.Assembly.buffers(cp.ws, cp.def, Float64)
        return ds
    end

    mat_k = _ea_material()
    ds_k = _derived_for(mat_k)
    @test length(ds_k.a) == 5
    @test ds_k.a == ds_k.ε
    for i = 1:5
        @test ds_k.ε[i] ≈ mixture_emissivity(_composition(i), _temperature(i), mat_k) rtol =
            1e-14
    end
    # ε varies with the cell temperature across the graded state.
    @test length(unique(round.(ds_k.ε; digits = 6))) == 5

    mat_a = _ea_material(a_char = 0.6)
    ds_a = _derived_for(mat_a)
    @test ds_a.a != ds_a.ε
    for i = 1:5
        ξ, T = _composition(i), _temperature(i)
        @test ds_a.a[i] ≈ mixture_absorptivity(ξ, T, mat_a) rtol = 1e-14
        @test ds_a.ε[i] ≈ mixture_emissivity(ξ, T, mat_a) rtol = 1e-14
    end
end
