"""Activity laws and the saturation pressure of water: published values, the
Clausius–Clapeyron identity of the heats, the limits, dual numbers, inference,
allocations and the constructors' refusals."""

using Test
using ForwardDiff
using Pyrolysis
using Pyrolysis.Materials: AbstractActivity, AbstractPropertyFunction, WaterSaturation
using Pyrolysis.Materials: activity, sorption_heat, R_GAS

const M_V = 0.018015          # molar mass of water vapour [kg/mol]
const P_ATM = 101325.0
const X_GRID = range(-0.01, 0.4; length = 83)
const T_GRID = range(260.0, 620.0; length = 37)
const X_FAR = (1.0, 5.0, 50.0)  # far above every cap

# CAST's differential Gibbs energy of sorption [kJ/mol], Dietenberger, Glass and
# Boardman 2024, eq. 18, written out here from the paper.
cast_gibbs(β, u, T) = (u * (β[1] * exp(-β[2] / T) + β[3]) + β[4])^(-β[5]) - β[6]
# The moisture content of the paper's law at activity h, eq. 18 inverted.
cast_moisture(β, h, T) =
    ((-1e-3 * R_GAS * T * log(h) + β[6])^(-1 / β[5]) - β[4]) /
    (β[1] * exp(-β[2] / T) + β[3])
# The moisture content at which Δg = 0, eq. 25.
cast_cap(β, T) = (β[6]^(-1 / β[5]) - β[4]) / (β[1] * exp(-β[2] / T) + β[3])

# Activity functions: CAST spruce as a function of (X, T) (its activity is
# CAST_SPRUCE's and its heat is eq. 19, here by forward differentiation); a law
# linear in X that does not depend on temperature; and one linear in X with a
# constant heat of sorption, 500 R/M.
const AF_CAST =
    ActivityFunction((X, T) -> exp(-1e3 * cast_gibbs(CAST_SPRUCE.β, X, T) / (R_GAS * T)))
const AF_LINEAR = ActivityFunction((X, T) -> X / 0.2)
const AF_HEAT = ActivityFunction((X, T) -> X / 0.25 * exp(500 * (1 / 298.15 - 1 / T)))

const LAWS = (
    ("CAST_SPRUCE", CAST_SPRUCE),
    ("CAST_BEECH", CAST_BEECH),
    ("FreeLiquid", FreeLiquid()),
    ("ActivityFunction, CAST spruce", AF_CAST),
    ("ActivityFunction, linear", AF_LINEAR),
    ("ActivityFunction, constant heat", AF_HEAT),
)

# The activity before the cap: what the law multiplies by X/(|X| + X0).
raw_activity(law::CASTIsotherm, X, T) =
    exp(-1e3 * cast_gibbs(law.β, max(X, 0.0), T) / (R_GAS * T)) * X / (abs(X) + law.X0)
raw_activity(law::ActivityFunction, X, T) = law.f(max(X, 0.0), T) * X / (abs(X) + law.X0)
raw_activity(law::FreeLiquid, X, T) = X / (abs(X) + law.X0)

# The cap at one, 1 − δ ln(1 + e^{(1 − x)/δ}), on the raw activity.
function capped(law, X, T)
    x = raw_activity(law, X, T)
    z = (1 - x) / law.δ
    return z > 40 ? x : 1 - law.δ * log1p(exp(z))
end

# Equilibrium moisture content of a law at activity h, by bisection.
function equilibrium_moisture(law, h, T)
    lo, hi = 0.0, 2.0
    for _ = 1:200
        mid = 0.5 * (lo + hi)
        activity(law, mid, T) < h ? (lo = mid) : (hi = mid)
    end
    return 0.5 * (lo + hi)
end

function allocations(law, X, T)
    activity(law, X, T)
    sorption_heat(law, X, T, M_V)
    return (@allocated activity(law, X, T)) + (@allocated sorption_heat(law, X, T, M_V))
end

function saturation_allocations(T)
    WATER_SATURATION(T)
    return @allocated WATER_SATURATION(T)
end

@testset "Activity laws" begin
    @testset "CAST constants (Table 3 of the paper)" begin
        @test CAST_SPRUCE.β == (304.2, 1661.0, 3.520, 0.3895, 2.886, 0.2006)
        @test CAST_BEECH.β == (46.33, 777.6, 2.028, 0.2210, 1.725, 0.3664)
        for law in (CAST_SPRUCE, CAST_BEECH, AF_CAST)
            @test law.δ == 0.01
            @test law.X0 == 1e-4
        end
        @test FreeLiquid().X0 == 1e-4
        @test CAST_SPRUCE isa AbstractActivity
        @test FreeLiquid() isa AbstractActivity
        @test AF_CAST isa AbstractActivity
    end

    @testset "the activity is the paper's law, multiplied and capped" begin
        for law in (CAST_SPRUCE, CAST_BEECH)
            @test all(
                isapprox(activity(law, X, T), capped(law, X, T); rtol = 1e-13) for
                X in X_GRID, T in T_GRID
            )
        end
        # the function form of CAST spruce gives CAST_SPRUCE's activity
        @test all(
            isapprox(activity(AF_CAST, X, T), activity(CAST_SPRUCE, X, T); rtol = 1e-14) for
            X in X_GRID, T in T_GRID
        )
    end

    @testset "published values" begin
        # Moisture content (% dry basis) of wood in steam at 1 atm, h = p/p_sat(T): the
        # paper's law, as printed to 0.1 % in the model's evidence.
        isobar = (
            (105, 10.4, 9.8),
            (110, 7.1, 7.0),
            (120, 4.2, 4.4),
            (130, 2.9, 3.1),
            (140, 2.0, 2.2),
            (150, 1.5, 1.7),
        )
        for (Tc, beech, spruce) in isobar
            T = Tc + 273.15
            h = P_ATM / WATER_SATURATION(T)
            for (law, printed) in ((CAST_BEECH, beech), (CAST_SPRUCE, spruce))
                # the paper's law, to the printed digits
                X_paper = cast_moisture(law.β, h, T)
                @test abs(100X_paper - printed) <= 0.05
                # the law as built: 10⁻³ in X, and above the paper's by what the
                # vanishing factor takes (X0/X of the activity)
                X_law = equilibrium_moisture(law, h, T)
                @test abs(X_law - printed / 100) <= 1e-3
                @test 0 < X_law - X_paper < 3e-4
            end
        end
        # beech at a = 0.5 and 25 °C
        @test abs(equilibrium_moisture(CAST_BEECH, 0.5, 298.15) - 0.0795) <= 5e-5
        # the free water of spruce sets in at 29.0, 19.2 and 14.2 % (Δg = 0, eq. 25)
        for (Tc, printed) in ((25, 29.0), (100, 19.2), (150, 14.2))
            T = Tc + 273.15
            X_f = cast_cap(CAST_SPRUCE.β, T)
            @test abs(100X_f - printed) <= 0.05
            @test abs(Pyrolysis.Materials._cast_terms(CAST_SPRUCE, X_f, T)[4]) < 1e-9
        end
    end

    @testset "sorption heat = (R T²/M) ∂ln a/∂T" begin
        for (name, law) in LAWS
            worst = 0.0
            for X in X_GRID, T in T_GRID
                activity(law, X, T) > 0 || continue
                ref =
                    R_GAS * T^2 / M_V *
                    ForwardDiff.derivative(t -> log(activity(law, X, t)), T)
                q = sorption_heat(law, X, T, M_V)
                worst = max(worst, ref == q ? 0.0 : abs(q - ref) / abs(ref))
            end
            @test worst <= 1e-10
        end
        @test all(sorption_heat(FreeLiquid(), X, T, M_V) == 0 for X in X_GRID, T in T_GRID)
        @test all(sorption_heat(AF_LINEAR, X, T, M_V) == 0 for X in X_GRID, T in T_GRID)
        # a law with a constant heat of sorption below its cap
        @test sorption_heat(AF_HEAT, 0.05, 300.0, M_V) ≈ 500 * R_GAS / M_V rtol = 1e-12
        # CAST's heat at 25 °C, as printed in the model's evidence (kJ/kg, to its last digit)
        for (law, printed) in (
            (CAST_SPRUCE, (791, 530, 263, 142, 83)),
            (CAST_BEECH, (706, 460, 245, 148, 96)),
        )
            for (X, Q) in zip((0.02, 0.05, 0.10, 0.15, 0.20), printed)
                @test abs(sorption_heat(law, X, 298.15, M_V) / 1e3 - Q) <= 1
            end
        end
    end

    @testset "limits" begin
        for (name, law) in LAWS, T in T_GRID
            @test activity(law, 0.0, T) == 0
            @test all(activity(law, X, T) <= 1 for X in (X_GRID..., X_FAR...))
        end
        # below zero the activity is negative where the law is positive at X = 0,
        # and zero where it vanishes there
        for law in (CAST_SPRUCE, CAST_BEECH, FreeLiquid(), AF_CAST), T in T_GRID
            @test all(activity(law, X, T) < 0 for X in X_GRID if X < 0)
        end
        for law in (AF_LINEAR, AF_HEAT), T in T_GRID
            @test all(activity(law, X, T) == 0 for X in X_GRID if X < 0)
        end
        # the cap: free water where the raw activity is well above one
        for law in (CAST_SPRUCE, CAST_BEECH, AF_CAST, AF_LINEAR, AF_HEAT),
            X in (X_GRID..., X_FAR...),
            T in T_GRID

            raw_activity(law, X, T) > 1.05 || continue
            @test 1 - activity(law, X, T) < 1e-4
        end
        @test activity(AF_LINEAR, 0.3, 300.0) == 1
        # the free liquid: one to within X0/X, no cap
        fl = FreeLiquid()
        @test activity(fl, 0.25, 373.0) ≈ 0.25 / 0.2501 rtol = 1e-15
        @test activity(fl, 0.01, 373.0) ≈ 0.01 / 0.0101 rtol = 1e-15
        @test all(activity(fl, X, 300.0) < 1 for X in (0.4, 1.0, 1e3, 1e9))
        @test activity(fl, -0.25, 373.0) ≈ -0.25 / 0.2501 rtol = 1e-15
        @test activity(FreeLiquid(; X0 = 1e-6), 0.25, 300.0) ≈ 0.25 / (0.25 + 1e-6) rtol =
            1e-15
        # C¹ through X = 0: the slope there is the raw activity of the bracket over X0
        for law in (CAST_SPRUCE, FreeLiquid()), T in (300.0, 450.0)
            left = ForwardDiff.derivative(X -> activity(law, X, T), -1e-12)
            right = ForwardDiff.derivative(X -> activity(law, X, T), 1e-12)
            @test left ≈ right rtol = 1e-6
        end
    end

    @testset "dual numbers" begin
        for (name, law) in LAWS
            bad = Tuple{Float64,Float64}[]
            for X in (X_GRID..., X_FAR...), T in T_GRID
                ga = ForwardDiff.gradient(v -> activity(law, v[1], v[2]), [X, T])
                gq = ForwardDiff.gradient(v -> sorption_heat(law, v[1], v[2], M_V), [X, T])
                finite =
                    all(isfinite, ga) &&
                    all(isfinite, gq) &&
                    isfinite(activity(law, X, T)) &&
                    isfinite(sorption_heat(law, X, T, M_V))
                finite || push!(bad, (X, T))
            end
            isempty(bad) || @info "non-finite value or derivative" name first(bad)
            @test isempty(bad)
        end
        # The heat of an activity function differentiates inside the outer dual
        # numbers: against central differences in X and in T, on the scales Q_s/X
        # and Q_s/T of the derivatives.
        for law in (AF_CAST, AF_HEAT, CAST_SPRUCE), (X, T) in ((0.05, 300.0), (0.12, 400.0))
            Q = sorption_heat(law, X, T, M_V)
            dX = ForwardDiff.derivative(x -> sorption_heat(law, x, T, M_V), X)
            dT = ForwardDiff.derivative(t -> sorption_heat(law, X, t, M_V), T)
            hX, hT = 1e-6, 1e-3
            fdX =
                (sorption_heat(law, X + hX, T, M_V) - sorption_heat(law, X - hX, T, M_V)) /
                (2hX)
            fdT =
                (sorption_heat(law, X, T + hT, M_V) - sorption_heat(law, X, T - hT, M_V)) /
                (2hT)
            @test isapprox(dX, fdX; atol = 1e-6 * abs(Q) / X)
            @test isapprox(dT, fdT; atol = 1e-6 * abs(Q) / T)
        end
        # Under the solver's Jacobian tag, the function form of CAST spruce has the
        # closed form's value and derivatives.
        J = Pyrolysis.ColoredJacobian.JacobianTag
        for (X, T) in ((0.05, 300.0), (0.12, 400.0), (0.25, 330.0), (-0.002, 350.0))
            Xd = ForwardDiff.Dual{J}(X, 1.0, 0.0)
            Td = ForwardDiff.Dual{J}(T, 0.0, 1.0)
            for f in (activity, (law, X, T) -> sorption_heat(law, X, T, M_V))
                a, b = f(AF_CAST, Xd, Td), f(CAST_SPRUCE, Xd, Td)
                @test ForwardDiff.value(a) ≈ ForwardDiff.value(b) rtol = 1e-12
                # (below zero CAST's heat is Δg(0)/M, constant in T: its T-partial is
                # zero in closed form and round-off in the function form)
                @test all(
                    isapprox.(
                        ForwardDiff.partials(a).values,
                        ForwardDiff.partials(b).values;
                        rtol = 1e-9,
                        atol = 1e-12 * abs(ForwardDiff.value(b)),
                    ),
                )
            end
        end
        # ... and agree with the closed form where the law is CAST spruce
        for (X, T) in ((0.05, 300.0), (0.12, 400.0), (0.25, 330.0))
            for i = 1:2
                f(law) = v -> sorption_heat(law, v[1], v[2], M_V)
                g_fun = ForwardDiff.gradient(f(AF_CAST), [X, T])
                g_cast = ForwardDiff.gradient(f(CAST_SPRUCE), [X, T])
                @test g_fun[i] ≈ g_cast[i] rtol = 1e-9
            end
        end
    end

    @testset "inference and allocations" begin
        J = Pyrolysis.ColoredJacobian.JacobianTag
        dX = ForwardDiff.Dual{J}(0.1, 1.0)
        dT = ForwardDiff.Dual{J}(350.0, 1.0)
        for (name, law) in LAWS
            for (X, T) in ((0.1, 350.0), (dX, 350.0), (0.1, dT))
                @inferred activity(law, X, T)
                @inferred sorption_heat(law, X, T, M_V)
            end
            for (X, T) in
                ((0.1, 350.0), (-0.01, 300.0), (5.0, 450.0), (dX, 350.0), (0.1, dT))
                @test allocations(law, X, T) == 0
            end
        end
        @inferred WATER_SATURATION(350.0)
        @inferred WATER_SATURATION(dT)
        @test saturation_allocations(350.0) == 0
        @test saturation_allocations(dT) == 0
    end

    @testset "constructor refusals" begin
        β = (304.2, 1661.0, 3.520, 0.3895, 2.886, 0.2006)
        for i = 1:6, bad in (NaN, Inf, 0.0, -1.0)
            b = Base.setindex(β, bad, i)
            @test_throws ArgumentError CASTIsotherm(b...)
        end
        for bad in (NaN, Inf, 0.0, -0.01)
            @test_throws ArgumentError CASTIsotherm(β...; δ = bad)
            @test_throws ArgumentError CASTIsotherm(β...; X0 = bad)
            @test_throws ArgumentError FreeLiquid(; X0 = bad)
            @test_throws ArgumentError ActivityFunction((X, T) -> X; δ = bad)
            @test_throws ArgumentError ActivityFunction((X, T) -> X; X0 = bad)
        end
        @test CASTIsotherm(β...; δ = 0.02, X0 = 1e-5).δ == 0.02
    end
end

@testset "Saturation pressure of water" begin
    @test WATER_SATURATION isa WaterSaturation
    @test WATER_SATURATION isa AbstractPropertyFunction
    @test abs(WATER_SATURATION(273.16) - 611.66) <= 0.005
    @test abs(WATER_SATURATION(298.15) - 3169.8) <= 0.05
    @test abs(WATER_SATURATION(373.15) - 101418) <= 0.5
    # held constant outside [255, 640] K, with C¹ joins 1 K wide
    low, high = WATER_SATURATION(255.0), WATER_SATURATION(640.0)
    @test all(WATER_SATURATION(T) == low for T in (-50.0, 0.0, 200.0, 254.0, 254.99))
    @test all(WATER_SATURATION(T) == high for T in (640.01, 641.0, 1000.0, 3000.0))
    @test ForwardDiff.derivative(WATER_SATURATION, 254.0) == 0
    @test ForwardDiff.derivative(WATER_SATURATION, 641.0) == 0
    for T in (255.0, 255.5, 256.0, 639.0, 639.5, 640.0)
        @test isfinite(ForwardDiff.derivative(WATER_SATURATION, T))
    end
    # increasing on the curve, and the identity in temperature away from its ends
    @test all(
        ForwardDiff.derivative(WATER_SATURATION, T) > 0 for
        T in range(256.0, 639.0; length = 50)
    )
    # evaluated on the state like every property law
    @test WATER_SATURATION(350.0, (1.0, 2.0)) == WATER_SATURATION(350.0)
end
