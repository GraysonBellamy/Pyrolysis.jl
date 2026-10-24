"""
Blowing models against their reference values.

The property-ratio and pooling forms were fitted to resolved simulations of
the CAPA chamber and are implemented independently in the Python calculator
`blowing_factor.py` of the oxidative_pyrolysis repository.
The reference numbers below come from that calculator with the same gas data
(NASA-7 heat capacities, CFD-case conductivity and viscosity polynomials,
Fuller diffusion volumes). The one deliberate difference is the mean free-
stream heat capacity between the gas and the surface temperature, which the
calculator takes as the mean of 41 equally spaced samples and this package
integrates exactly (5-point Gauss–Legendre), so the comparison is to 1e-4.

1. Gas data: `GasSpecies`, `NITROGEN`, `nasa7_heat_capacity`, Eucken and
   Fuller helpers against tabulated values.
2. `PropertyRatioFilm`: limits, reduction to the Couette film, reference ratios.
3. `PoolingFilm`: reference ratios, pooling factor, face fraction and buoyancy
   for a light gas (surrogate volatiles), water vapour, nitrogen and methyl
   methacrylate on both sides of the switch; limits; AD through the face solve.
4. Per-species mixture path: a single component reproduces the lumped
   description; Wilke mixing lies between the pure-species conductivities.
"""

using Test
using ForwardDiff
using Pyrolysis
using Pyrolysis.Materials
using Pyrolysis.BoundaryConditions
using Pyrolysis.BoundaryConditions:
    BlowingFace,
    blowing_state,
    blowing_ratio,
    blowing_factor,
    solve_boundary_temperature!,
    BoundaryState
using Pyrolysis.States: trial_temperature

# ----------------------------------------------------------------------------
# Gas data of the calculator (blowing_gases.json): NASA-7 cp/R low set,
# conductivity and viscosity polynomials in SI, Fuller diffusion volumes.
# ----------------------------------------------------------------------------
_poly(c) = PolynomialProperty(c)

# The calculator's binary diffusivity into nitrogen (its Fuller constant is
# 1e-7 with M in g/mol and no pressure factor; the package's helper carries the
# textbook 0.00143/P), given explicitly so the comparison isolates the model.
function _calc_diffusivity(M, V)
    Mg, Mn = 1e3 * M, 28.014
    Vn = 18.5
    return CallableProperty(
        T -> 1e-7 * T^1.75 * sqrt(1 / Mg + 1 / Mn) / (cbrt(V) + cbrt(Vn))^2,
    )
end

const SURROGATE = GasSpecies(
    :pyrolysate;
    M = 0.042,
    c = nasa7_heat_capacity(
        (
            3.2530406371,
            0.0096897643924,
            3.4722749667e-06,
            -7.8787625275e-09,
            2.8176572383e-12,
        ),
        0.042,
    ),
    k = _poly((
        0.0019742358149,
        2.1293609878e-05,
        1.6502619984e-07,
        -1.1903863518e-10,
        3.1187193584e-14,
    )),
    μ = _poly((
        -2.5907503777e-07,
        3.6642563501e-08,
        -6.5443252412e-14,
        -6.3521134277e-15,
        2.3401486488e-18,
    )),
    V_D = 65.4,
    D = _calc_diffusivity(0.042, 65.4),
    T_range = (200.0, 1250.0),
)
const MMA = GasSpecies(
    :MMA;
    M = 0.1001158,
    c = nasa7_heat_capacity(
        (
            -1.7019146976,
            0.063342153415,
            -3.8278617755e-05,
            8.2820920523e-09,
            -1.8859789348e-13,
        ),
        0.1001158,
    ),
    k = _poly((
        -0.0038579814074,
        1.4610157687e-05,
        1.4629652404e-07,
        -1.0909786588e-10,
        2.8445996362e-14,
    )),
    μ = _poly((
        -1.6197096892e-06,
        3.1780836895e-08,
        -5.725870742e-12,
        -3.8291292585e-16,
        5.0658457784e-19,
    )),
    V_D = 110.2,
    D = _calc_diffusivity(0.1001158, 110.2),
    T_range = (200.0, 1250.0),
)
const WATER = GasSpecies(
    :H2O;
    M = 0.0180153,
    c = nasa7_heat_capacity(
        (
            3.9677872019,
            -0.00045567811047,
            2.674206457e-06,
            -1.5241823097e-09,
            3.0155465519e-13,
        ),
        0.0180153,
    ),
    k = _poly((
        -0.010726753419,
        0.00010493612203,
        -6.3041297194e-08,
        1.011413797e-10,
        -3.6459810842e-14,
    )),
    μ = _poly((
        7.606206286e-06,
        -1.8920862849e-08,
        1.1670130357e-10,
        -9.4505713056e-14,
        2.6753982374e-17,
    )),
    V_D = 13.1,
    D = _calc_diffusivity(0.0180153, 13.1),
    T_range = (200.0, 1250.0),
)
# Nitrogen blown into nitrogen: the free stream itself, with the calculator's D.
const N2_BLOWN = GasSpecies(
    :N2_blown;
    M = NITROGEN.molar_mass,
    c = NITROGEN.heat_capacity,
    k = NITROGEN.thermal_conductivity,
    μ = NITROGEN.viscosity,
    V_D = 18.5,
    D = _calc_diffusivity(0.028014, 18.5),
    T_range = (200.0, 1250.0),
)

# A material whose one gas component stands for the lumped outflow; its own
# properties are never read when `volatiles` is given.
function _carrier_material(; cp = 1000.0, k = 0.02, M = 0.050)
    solid = SolidComponent(:solid; ρ = 1000.0, c = 1000.0, k = 0.2, ε = 0.9)
    gas = GaseousComponent(:gas; M = M, c = cp, k = k, λ = 1e-5)
    rxn = Reaction(:decomp, 1 => 2, A = 1e5, E = 5e4, h = 1e5)
    return Material(name = :carrier, components = (solid, gas), reactions = (rxn,))
end

_face(mat, ṁ_g_m2s, T_s) = BlowingFace((0.0, ṁ_g_m2s * 1e-3), T_s, mat)

# Reference rows: gas, ṁ″ [g/(m²·s)], T_s, T_in, h₀ → ratio, light-gas form,
# φ, ψ, face fraction, buoyancy (calculator output, 10 digits).
const REFERENCE = (
    (
        SURROGATE,
        8.0,
        870.0,
        360.0,
        11.2,
        0.5948057161,
        0.5986304188,
        0.9513589228,
        0.9958767108,
        0.6396778633,
        0.8089398248,
    ),
    (
        SURROGATE,
        3.0,
        650.0,
        336.0,
        11.2,
        0.8311061527,
        0.8332662806,
        0.3543512404,
        0.9978171245,
        0.3185579802,
        0.8730175786,
    ),
    (
        MMA,
        8.0,
        870.0,
        360.0,
        11.2,
        0.1752582883,
        0.5068607936,
        1.2337616673,
        0.5919057109,
        0.8706383645,
        -0.1866643156,
    ),
    (
        MMA,
        16.0,
        650.0,
        336.0,
        9.591,
        0.0215950821,
        0.1695576891,
        2.8942062748,
        0.5900000007,
        0.9918584978,
        -1.6756172004,
    ),
    (
        MMA,
        3.0,
        870.0,
        360.0,
        11.2,
        0.7789335449,
        0.7864442763,
        0.4626606252,
        0.9923532860,
        0.3671029426,
        0.7463054509,
    ),
    (
        WATER,
        8.0,
        870.0,
        360.0,
        11.2,
        0.5760944066,
        0.5762177221,
        1.0170682035,
        0.9998656795,
        0.4958439603,
        1.1523355910,
    ),
    (
        N2_BLOWN,
        8.0,
        870.0,
        360.0,
        11.2,
        0.6470656999,
        0.6476482841,
        0.8139121718,
        0.9993845153,
        0.5371363753,
        1.0000000000,
    ),
)

@testset "Blowing models" begin

    # ========================================================================
    @testset "gas data" begin
        # NASA-7 → cp per unit mass: nitrogen at 870 K from the calculator.
        @test NITROGEN.heat_capacity(870.0) ≈ 1139.477041 rtol = 1e-8
        @test NITROGEN.thermal_conductivity(870.0) ≈ 0.05902199 rtol = 1e-6
        @test SURROGATE.heat_capacity(870.0) ≈ 2125.600469 rtol = 1e-8
        @test SURROGATE.thermal_conductivity(870.0) ≈ 0.08488778 rtol = 1e-6
        @test MMA.heat_capacity(870.0) ≈ 2473.050823 rtol = 1e-8
        @test NITROGEN.name === :N2 && NITROGEN.molar_mass == 0.028014
        @test NITROGEN.diffusion_volume == 18.5

        # Eucken's relation is exact for a gas defined through it.
        k, cp, M = 0.05, 1200.0, 0.03
        μ = eucken_viscosity(k, cp, M)
        @test μ * (cp + 1.25 * Pyrolysis.Materials.R_GAS / M) ≈ k rtol = 1e-14
        # Nitrogen's tabulated viscosity and Eucken agree within the relation's
        # usual accuracy for a diatomic gas.
        μ_e = eucken_viscosity(
            NITROGEN.thermal_conductivity(500.0),
            NITROGEN.heat_capacity(500.0),
            NITROGEN.molar_mass,
        )
        @test μ_e ≈ NITROGEN.viscosity(500.0) rtol = 0.15

        # Fuller: CO₂ in N₂ at 298 K, measured 0.165 cm²/s (Marrero and Mason).
        D = fuller_diffusivity(298.15, 0.04401, 26.9, 0.028014, 18.5)
        @test D ≈ 1.62e-5 rtol = 0.02
        @test D ≈ 1.65e-5 rtol = 0.10
        # Inverse in pressure, T^1.75.
        @test fuller_diffusivity(298.15, 0.04401, 26.9, 0.028014, 18.5; P = 2 * 101325.0) ≈
              D / 2 rtol = 1e-12
        @test fuller_diffusivity(596.3, 0.04401, 26.9, 0.028014, 18.5) ≈ D * 2^1.75 rtol =
            1e-12

        @test_throws ArgumentError GasSpecies(:bad; M = 0.0, c = 1.0, k = 1.0)
        @test_throws ArgumentError GasSpecies(:bad; M = 0.03, c = 1.0, k = 1.0, V_D = -1.0)
        @test GasSpecies(:plain; M = 0.03, c = 1000.0, k = 0.03).viscosity === nothing
        @test GasSpecies(:plain; M = 0.03, c = 1000.0, k = 0.03).T_range == (200.0, 3000.0)
        # A fit that goes negative inside the stated range is rejected; held
        # inside its range it is accepted and evaluated at the range edge beyond it.
        kpoly = NITROGEN.thermal_conductivity
        @test_throws ArgumentError GasSpecies(:n2; M = 0.028, c = 1000.0, k = kpoly)
        gs = GasSpecies(:n2; M = 0.028, c = 1000.0, k = kpoly, T_range = (200.0, 1250.0))
        @test gas_property_temperature(gs, 1600.0) == 1250.0
        @test gas_property_temperature(gs, 150.0) == 200.0
        @test gas_property_temperature(gs, 700.0) == 700.0
        @test NITROGEN.T_range == (200.0, 1250.0)
        # Two-set NASA-7: continuous to the sets' agreement at the common temperature.
        @test NITROGEN.heat_capacity(999.999) ≈ NITROGEN.heat_capacity(1000.001) rtol = 1e-3
        @test NITROGEN.heat_capacity(1200.0) > NITROGEN.heat_capacity(1000.0)
    end

    # ========================================================================
    @testset "PropertyRatioFilm" begin
        mat = _carrier_material()
        h, T_inf = 11.2, 360.0
        @test_throws ArgumentError PropertyRatioFilm(NITROGEN; n_Γ = -1.0)
        @test PropertyRatioFilm(NITROGEN).n_Γ == 0.6

        # Zero outflow: exactly 1. Net inflow: the free-stream Couette film.
        m = PropertyRatioFilm(NITROGEN; volatiles = SURROGATE)
        @test blowing_ratio(m, _face(mat, 0.0, 870.0), h, T_inf) === 1.0
        s_in = blowing_state(m, _face(mat, -4.0, 870.0), h, T_inf)
        @test s_in.φ ≈ -4e-3 * NITROGEN.heat_capacity(870.0) / h rtol = 1e-14
        @test s_in.χ > 1

        # n_Γ = 0 with the free stream's own heat capacity on the gas is the
        # Couette film (λ = 0.5) in that heat capacity.
        mat_n2 = _carrier_material(cp = NITROGEN.heat_capacity(870.0))
        m0 = PropertyRatioFilm(NITROGEN; n_Γ = 0.0)
        f = _face(mat_n2, 8.0, 870.0)
        @test blowing_ratio(m0, f, h, T_inf) ≈ blowing_ratio(CouetteFilm(), f, h, T_inf) rtol =
            1e-14

        # Reference light-gas form (ψ = 1) for every row.
        for (gas, ṁ, T_s, T_in, h0, _, form, φ0, _, _, _) in REFERENCE
            s = blowing_state(
                PropertyRatioFilm(NITROGEN; volatiles = gas),
                _face(mat, ṁ, T_s),
                h0,
                T_in,
            )
            @test s.φ ≈ φ0 rtol = 1e-7
            @test s.χ ≈ form rtol = 1e-7
            @test s.ψ == 1 && s.Y == 0 && s.states == 1
        end

        # Monotone in the outflow.
        χs = [
            blowing_ratio(m, _face(mat, ṁ, 870.0), h, T_inf) for
            ṁ in (1.0, 3.0, 8.0, 16.0, 24.0)
        ]
        @test all(diff(χs) .< 0)
    end

    # ========================================================================
    @testset "PoolingFilm: reference values" begin
        mat = _carrier_material()
        @test_throws ArgumentError PoolingFilm(NITROGEN; ψ0 = 1.0)
        @test_throws ArgumentError PoolingFilm(NITROGEN; Λw = 0.0)
        @test_throws ArgumentError PoolingFilm(NITROGEN; diffusion_volume = (gas = -1.0,))
        d = PoolingFilm(NITROGEN)
        @test (d.n_Γ, d.n_D, d.ψ0, d.Λc, d.Λw) == (0.6, 0.5, 0.59, 0.35, 0.10)

        for (gas, ṁ, T_s, T_in, h0, ratio, _, φ0, ψ, Y, Λ) in REFERENCE
            s = blowing_state(
                PoolingFilm(NITROGEN; volatiles = gas),
                _face(mat, ṁ, T_s),
                h0,
                T_in,
            )
            @test s.φ ≈ φ0 rtol = 1e-7
            @test s.χ ≈ ratio rtol = 1e-4
            @test s.ψ ≈ ψ rtol = 1e-4
            @test s.Y ≈ Y rtol = 5e-4      # the face balance carries the c̄_p quadrature difference
            @test s.Λ ≈ Λ atol = 5e-4      # and Λ amplifies it by ∂Λ/∂Y ≈ −4 for a heavy gas
            @test s.states == 1
            # The ratio is the film factor at φ/ψ, scaled by ψ.
            @test s.χ ≈ s.ψ * blowing_factor(s.φ / s.ψ) rtol = 1e-12
        end

        # Limits: zero outflow is exactly 1; a face no hotter than the gas
        # keeps ψ = 1; suction is the free-stream film.
        m = PoolingFilm(NITROGEN; volatiles = MMA)
        @test blowing_ratio(m, _face(mat, 0.0, 870.0), 11.2, 360.0) === 1.0
        s_cold = blowing_state(m, _face(mat, 8.0, 350.0), 11.2, 360.0)
        @test s_cold.ψ == 1
        @test s_cold.χ ≈ blowing_factor(s_cold.φ) rtol = 1e-14
        s_in = blowing_state(m, _face(mat, -8.0, 870.0), 11.2, 360.0)
        @test s_in.ψ == 1 && s_in.χ > 1

        # Pooled gas keeps less exchange than the light-gas form at the same φ,
        # and the floor is approached for a heavy gas at high outflow.
        s_pool = blowing_state(m, _face(mat, 16.0, 650.0), 9.591, 336.0)
        @test s_pool.χ < blowing_factor(s_pool.φ)
        @test s_pool.ψ ≈ 0.59 atol = 1e-6
    end

    # ========================================================================
    @testset "PoolingFilm: AD through the face solve" begin
        mat = _carrier_material()
        m = PoolingFilm(NITROGEN; volatiles = MMA)
        h, T_inf = 11.2, 360.0
        # On the switch: MMA at 5 g/(m²·s) and 870 K sits between the lifted
        # and pooled states, where ∂ψ/∂ṁ and ∂ψ/∂T_s are largest.
        χ_of_ṁ =
            ṁ -> blowing_ratio(m, BlowingFace((zero(ṁ), ṁ), 870.0 + zero(ṁ), mat), h, T_inf)
        χ_of_T =
            T -> blowing_ratio(m, BlowingFace((zero(T), 5e-3 + zero(T)), T, mat), h, T_inf)
        ṁ0 = 5e-3
        d_ad = ForwardDiff.derivative(χ_of_ṁ, ṁ0)
        δ = 1e-7
        d_fd = (χ_of_ṁ(ṁ0 + δ) - χ_of_ṁ(ṁ0 - δ)) / (2δ)
        @test d_ad ≈ d_fd rtol = 1e-5
        dT_ad = ForwardDiff.derivative(χ_of_T, 870.0)
        δT = 1e-3
        dT_fd = (χ_of_T(870.0 + δT) - χ_of_T(870.0 - δT)) / (2δT)
        @test dT_ad ≈ dT_fd rtol = 1e-5
        # And the pooling factor itself is on the switch here.
        s = blowing_state(m, _face(mat, 5.0, 870.0), h, T_inf)
        @test 0.6 < s.ψ < 0.99

        # Through the surface Newton: ∂T_s/∂ṁ″ matches finite differences.
        bc = HeatFluxBC(3e4) + ConvectiveBC(h = h, T_inf = T_inf, blowing = m)
        k_eff, Δz_half, t, ε = 0.2, 5e-4, 0.0, 0.9
        function _Ts(J::Real)
            st = BoundaryState()
            r = solve_boundary_temperature!(
                st,
                800.0,
                k_eff,
                Δz_half,
                bc,
                mat,
                (zero(J), J),
                t,
                ε;
                trial = trial_temperature,
            )
            return r.T_s
        end
        dJ_ad = ForwardDiff.derivative(_Ts, 5e-3)
        δJ = 1e-7
        dJ_fd = (_Ts(5e-3 + δJ) - _Ts(5e-3 - δJ)) / (2δJ)
        @test dJ_ad ≈ dJ_fd rtol = 1e-5
    end

    # ========================================================================
    @testset "per-species mixture path" begin
        h, T_inf = 11.2, 360.0
        T_s = 870.0
        # One gas component with the surrogate's properties reproduces the
        # lumped description: same Γ ratio, same molar mass, same diffusivity
        # (its Fuller volume on the model, the calculator's D on the lumped gas
        # differ by the Fuller constant only).
        solid = SolidComponent(:solid; ρ = 1000.0, c = 1000.0, k = 0.2, ε = 0.9)
        g1 = GaseousComponent(
            :pyro;
            M = 0.042,
            c = SURROGATE.heat_capacity,
            k = SURROGATE.thermal_conductivity,
            λ = 1e-5,
        )
        mat1 = Material(
            name = :one,
            components = (solid, g1),
            reactions = (Reaction(:d, 1 => 2, A = 1e5, E = 5e4, h = 1e5),),
        )
        face1 = BlowingFace((0.0, 8e-3), T_s, mat1)
        s_mix = blowing_state(PropertyRatioFilm(NITROGEN), face1, h, T_inf)
        s_lump = blowing_state(
            PropertyRatioFilm(NITROGEN; volatiles = SURROGATE),
            face1,
            h,
            T_inf,
        )
        @test s_mix.φ ≈ s_lump.φ rtol = 1e-12
        p_mix = blowing_state(
            PoolingFilm(NITROGEN; diffusion_volume = (pyro = 65.4,)),
            face1,
            h,
            T_inf,
        )
        p_lump =
            blowing_state(PoolingFilm(NITROGEN; volatiles = SURROGATE), face1, h, T_inf)
        @test p_mix.χ ≈ p_lump.χ rtol = 2e-3
        @test p_mix.Y ≈ p_lump.Y rtol = 2e-3

        # Free-stream gas venting from the pores is purge at the face: it
        # counts in the outflow but not in the face composition.
        # The CFD-case polynomials are fitted below 1200 K; the component validator
        # samples the whole property bracket, so skip it for this fixture.
        n2 = GaseousComponent(
            :N2;
            M = 0.028014,
            c = NITROGEN.heat_capacity,
            k = NITROGEN.thermal_conductivity,
            λ = 1e-5,
            validate_properties = false,
        )
        mat2 = Material(
            name = :two,
            components = (solid, g1, n2),
            reactions = (Reaction(:d, 1 => 2, A = 1e5, E = 5e4, h = 1e5),),
        )
        m2 = PoolingFilm(NITROGEN; diffusion_volume = (pyro = 65.4,))
        s_n2only = blowing_state(m2, BlowingFace((0.0, 0.0, 8e-3), T_s, mat2), h, T_inf)
        @test s_n2only.Y == 0
        @test s_n2only.Λ ≈ 1 atol = 1e-12
        # The logistic switch does not reach 1 at Λ = 1: ψ(1) = ψ0 + (1 − ψ0)/(1 + e^{−6.5}),
        # the equation's own value for pure purge at the face.
        ψ1 = 0.59 + 0.41 / (1 + exp(-(1 - 0.35) / 0.10))
        @test s_n2only.ψ ≈ ψ1 rtol = 1e-10
        φ_n2 = 8e-3 * NITROGEN.heat_capacity(T_s) / h
        @test s_n2only.φ ≈ φ_n2 rtol = 1e-12
        @test s_n2only.χ ≈ ψ1 * blowing_factor(φ_n2 / ψ1) rtol = 1e-10
        s_both = blowing_state(m2, BlowingFace((0.0, 8e-3, 2e-3), T_s, mat2), h, T_inf)
        @test 0 < s_both.Y < p_mix.Y          # the nitrogen share dilutes the face
        @test s_both.Λ > p_mix.Λ

        # Wilke mixing: two identical species give the species' conductivity;
        # two different ones lie between them.
        mat_same = Material(
            name = :same,
            components = (
                solid,
                g1,
                GaseousComponent(
                    :pyro2;
                    M = 0.042,
                    c = SURROGATE.heat_capacity,
                    k = SURROGATE.thermal_conductivity,
                    λ = 1e-5,
                ),
            ),
            reactions = (Reaction(:d, 1 => 2, A = 1e5, E = 5e4, h = 1e5),),
        )
        s_same = blowing_state(
            PropertyRatioFilm(NITROGEN),
            BlowingFace((0.0, 5e-3, 3e-3), T_s, mat_same),
            h,
            T_inf,
        )
        @test s_same.φ ≈ s_lump.φ rtol = 1e-12
        cp_m, k_m = Pyrolysis.BoundaryConditions._leaving_mixture(
            nothing,
            BlowingFace((0.0, 5e-3, 3e-3), T_s, mat2),
        )
        @test min(SURROGATE.thermal_conductivity(T_s), NITROGEN.thermal_conductivity(T_s)) <
              k_m <
              max(SURROGATE.thermal_conductivity(T_s), NITROGEN.thermal_conductivity(T_s))
        @test cp_m ≈
              (5 * SURROGATE.heat_capacity(T_s) + 3 * NITROGEN.heat_capacity(T_s)) / 8 rtol =
            1e-12

        # Components flowing inward do not enter the leaving mixture.
        s_in = blowing_state(
            PropertyRatioFilm(NITROGEN),
            BlowingFace((0.0, 8e-3, -1e-3), T_s, mat2),
            h,
            T_inf,
        )
        @test s_in.φ ≈ 7e-3 / 8e-3 * s_mix.φ rtol = 1e-12
    end
end
