"""The declaration of a phase change: resolution against the components, the
refusals, its display, and the rules that read what each component leaves."""

using Test
using Pyrolysis
using Pyrolysis.Materials: n_phase_changes, dry_gas_flags
using Pyrolysis.Closures: residues, retained_volumes

const WOOD =
    SolidComponent(:wood; ρ = 400.0, ρ_intrinsic = 1500.0, c = 1500.0, k = 0.12, λ = 5e-7)
const CHAR =
    SolidComponent(:char; ρ = 100.0, ρ_intrinsic = 1800.0, c = 1100.0, k = 0.08, λ = 1e-6)
const WATER = LiquidComponent(:water; ρ = 1000.0, c = 4180.0, k = 0.6)
const VAPOUR = GaseousComponent(:vapour; M = 0.018015, c = 1900.0, k = 0.025)
const GAS = GaseousComponent(:gas; M = 0.03, c = 1100.0, k = 0.03)
const COMPONENTS = (WOOD, CHAR, WATER, VAPOUR, GAS)
const PYROLYSIS = Reaction(
    :pyro,
    :wood => (:char, :gas);
    A = 1e8,
    E = 1.2e5,
    h = 0.0,
    yields = (0.25, 0.75),
)

_wet(phase_changes; reactions = (PYROLYSIS,), components = COMPONENTS) =
    Material(; name = :wet, components, reactions, phase_changes)

_moisture(pair = :water => :vapour; kwargs...) =
    PhaseChange(:moisture, pair; activity = CAST_SPRUCE, sorbents = (:wood,), kwargs...)

_refusal(f) =
    try
        f()
        nothing
    catch err
        err
    end

@testset "Phase change: declaration" begin
    @testset "by name and by index the same resolved value" begin
        by_name = _wet((
            PhaseChange(
                :moisture,
                :water => :vapour;
                activity = CAST_SPRUCE,
                sorbents = (:wood, :char => 0.5),
                exchange_rate = 2e3,
            ),
        ))
        by_index = _wet((
            PhaseChange(
                :moisture,
                3 => 4;
                activity = CAST_SPRUCE,
                sorbents = (1 => 1.0, 2 => 0.5),
                exchange_rate = 2e3,
            ),
        ))
        pc = only(by_name.phase_changes)
        @test pc isa PhaseChange
        @test pc == only(by_index.phase_changes)
        @test n_phase_changes(by_name) == 1
        @test n_phase_changes(_wet(())) == 0
        @test (pc.name, pc.liquid, pc.vapour) == (:moisture, 3, 4)
        @test pc.sorbents == ((1, 1.0), (2, 0.5))
        @test pc.activity === CAST_SPRUCE
        @test pc.saturation === WATER_SATURATION
        @test pc.latent_heat == 2.4417e6
        @test pc.exchange_rate == 2e3
        # the vapour's molar mass is copied in from the component
        @test pc.molar_mass == 0.018015
        # defaults
        d = only(_wet((_moisture(),)).phase_changes)
        @test (d.latent_heat, d.exchange_rate, d.sorbents) == (2.4417e6, 1e3, ((1, 1.0),))
        # a single sorbent outside a tuple, and a saturation curve given as a function
        f = only(
            _wet((
                PhaseChange(
                    :moisture,
                    :water => :vapour;
                    activity = FreeLiquid(),
                    sorbents = :wood,
                    saturation = T -> 1e3 * exp(0.05 * (T - 300)),
                ),
            )).phase_changes,
        )
        @test f.sorbents == ((1, 1.0),)
        @test f.saturation isa CallableProperty
        @test f.saturation(300.0) == 1e3
        # a resolved phase change is resolved again against the new components
        again = _wet(by_name.phase_changes)
        @test only(again.phase_changes) == pc
        heavier = GaseousComponent(:vapour; M = 0.018, c = 1900.0, k = 0.025)
        swapped =
            _wet(by_name.phase_changes; components = (WOOD, CHAR, WATER, heavier, GAS))
        @test only(swapped.phase_changes).molar_mass == 0.018
        @test _refusal(
            () ->
                _wet(by_name.phase_changes; components = (WOOD, CHAR, VAPOUR, WATER, GAS)),
        ) isa ArgumentError
        # two phase changes of two liquids
        oil = LiquidComponent(:oil; ρ = 900.0, c = 2000.0, k = 0.15)
        oil_vapour = GaseousComponent(:oil_vapour; M = 0.15, c = 1800.0, k = 0.02)
        two = Material(
            name = :two,
            components = (COMPONENTS..., oil, oil_vapour),
            reactions = (PYROLYSIS,),
            phase_changes = (
                _moisture(),
                PhaseChange(
                    :oil,
                    :oil => :oil_vapour;
                    activity = FreeLiquid(),
                    sorbents = (:wood,),
                    saturation = T -> 10.0,
                    latent_heat = 3e5,
                ),
            ),
        )
        @test n_phase_changes(two) == 2
        @test map(pc -> (pc.liquid, pc.vapour, pc.molar_mass), two.phase_changes) ==
              ((3, 4, 0.018015), (6, 7, 0.15))
        # a liquid may also be the reactant or a product of a reaction
        hydrolysis = Reaction(:hydrolysis, :water => :gas; A = 1.0, E = 5e4, h = 0.0)
        wetting = Reaction(
            :wetting,
            :wood => (:char, :water);
            A = 1.0,
            E = 1e5,
            h = 0.0,
            yields = (0.5, 0.5),
        )
        @test n_phase_changes(_wet((_moisture(),); reactions = (PYROLYSIS, hydrolysis))) ==
              1
        @test n_phase_changes(_wet((_moisture(),); reactions = (wetting,))) == 1
    end

    @testset "refusals" begin
        refused(f, pattern) =
            (err = _refusal(f); err isa ArgumentError && occursin(pattern, err.msg))
        # the liquid is not a liquid, the vapour not a gas
        @test refused(
            () -> _wet((_moisture(:wood => :vapour),)),
            "phase change :moisture: component :wood is not a LiquidComponent",
        )
        @test refused(
            () -> _wet((_moisture(:water => :char),)),
            "phase change :moisture: component :char is not a GaseousComponent with a positive molar mass",
        )
        @test refused(
            () -> _wet((_moisture(3 => 1),)),
            "component :wood is not a GaseousComponent",
        )
        # unknown components
        @test refused(
            () -> _wet((_moisture(:steam => :vapour),)),
            "phase change :moisture references unknown component :steam",
        )
        @test refused(
            () -> _wet((_moisture(3 => 9),)),
            "phase change :moisture references component 9, outside 1:5",
        )
        @test refused(
            () -> _wet((_moisture(:water => "vapour"),)),
            "give a component name or index",
        )
        # sorbents
        sorbing(s) =
            () -> _wet((
                PhaseChange(
                    :moisture,
                    :water => :vapour;
                    activity = CAST_SPRUCE,
                    sorbents = s,
                ),
            ))
        @test refused(
            sorbing((:wood, :water)),
            "phase change :moisture: the sorbent :water is not a SolidComponent",
        )
        @test refused(sorbing((:wood, :gas)), "the sorbent :gas is not a SolidComponent")
        @test refused(sorbing((:wood => -1.0,)), "sorbent weight -1.0")
        @test refused(sorbing((:wood => Inf,)), "sorbent weight Inf")
        @test refused(sorbing((:wood => NaN,)), "sorbent weight NaN")
        @test refused(
            sorbing((:wood => 0.0, :char => 0.0)),
            "no sorbing solid has a positive weight",
        )
        @test refused(sorbing(()), "declares no sorbing solid")
        @test refused(sorbing((:wood, :wood => 0.5)), "names a sorbing solid twice")
        @test refused(sorbing((:bark,)), "references unknown component :bark")
        @test_throws UndefKeywordError PhaseChange(
            :moisture,
            :water => :vapour;
            activity = CAST_SPRUCE,
        )
        @test_throws UndefKeywordError PhaseChange(
            :moisture,
            :water => :vapour;
            sorbents = (:wood,),
        )
        @test_throws TypeError PhaseChange(
            :moisture,
            :water => :vapour;
            activity = 0.5,
            sorbents = (:wood,),
        )
        # two phase changes of one liquid
        @test refused(
            () -> _wet((
                _moisture(),
                PhaseChange(
                    :again,
                    :water => :vapour;
                    activity = FreeLiquid(),
                    sorbents = (:wood,),
                ),
            )),
            "the phase changes :moisture, :again share a liquid",
        )
        # the exchange rate and the latent heat
        for K in (0.0, -1.0, Inf, NaN)
            @test refused(
                () -> _wet((_moisture(; exchange_rate = K),)),
                "phase change :moisture has exchange_rate = $K",
            )
        end
        for L in (Inf, -Inf, NaN)
            @test refused(
                () -> _wet((_moisture(; latent_heat = L),)),
                "phase change :moisture has latent_heat = $L",
            )
        end
        # a negative or zero latent heat is a law, not a malformed declaration
        @test _wet((_moisture(; latent_heat = 0.0),)) isa Material
        # anything else in the tuple
        @test refused(
            () -> _wet((PYROLYSIS,)),
            "phase_changes must hold the declarations of PhaseChange",
        )
    end

    @testset "show" begin
        mat = _wet((_moisture(; exchange_rate = 2e3),))
        text = sprint(show, MIME"text/plain"(), mat)
        @test occursin("Reactions (1):\n    1. pyro", text)
        @test occursin(
            "Phase changes (1):\n    1. moisture: water ⇄ vapour, CASTIsotherm, sorbents wood, K=2000.0 1/s",
            text,
        )
        @test findfirst("Reactions", text).start < findfirst("Phase changes", text).start
        @test sprint(show, mat) ==
              "Material(:wet, 5 components, 1 reactions, 1 phase changes)"
        # a material without phase changes shows as before
        dry = _wet(())
        @test sprint(show, dry) == "Material(:wet, 5 components, 1 reactions)"
        @test !occursin("Phase changes", sprint(show, MIME"text/plain"(), dry))
    end
end

@testset "Phase change: what each component leaves" begin
    drying = Reaction(:drying, :water => :vapour; A = 1e5, E = 5e4, h = 2.3e6)
    wet = _wet((_moisture(),))
    twin = _wet((); reactions = (PYROLYSIS, drying))
    # the liquid leaves nothing, as under a one-way drying reaction
    @test residues(wet) == residues(twin)
    @test residues(wet)[3] == 0
    for T in (300.0, 450.0)
        @test retained_volumes(wet, T) == retained_volumes(twin, T)
        @test retained_volumes(wet, T)[3] == 0
    end
    # without the rule the water would be inert and leave all of itself
    inert = _wet(())
    @test residues(inert)[3] == 1
    @test retained_volumes(inert, 300.0)[3] == 1 / 1000
    # whatever reactions also consume it: a reaction of the water into char
    # leaves char, but the water can leave as vapour
    charring = Reaction(:charring, :water => :char; A = 1.0, E = 5e4, h = 0.0)
    both = _wet((_moisture(),); reactions = (PYROLYSIS, charring))
    @test residues(both)[3] == 0
    @test residues(_wet((); reactions = (PYROLYSIS, charring)))[3] == 1
    # a material that gasifies completely leaves nothing from a wet cell
    ablator = Material(
        name = :ablator,
        components = (WOOD, WATER, VAPOUR, GAS),
        reactions = (Reaction(:ablate, :wood => :gas; A = 1e8, E = 1.2e5, h = 0.0),),
        phase_changes = (_moisture(),),
    )
    @test residues(ablator) == (0.0, 0.0, 0.0, 0.0)
    @test retained_volumes(ablator, 300.0) == (0.0, 0.0, 0.0, 0.0)
end

@testset "Phase change: dry pyrolysis gas" begin
    # the vapour of a phase change is not dry gas; the pyrolysis gas is
    @test dry_gas_flags(_wet((_moisture(),))) == (false, false, false, false, true)
    # a solid's reaction producing the vapour as well: two origins
    wet_pyrolysis = Reaction(
        :pyro,
        :wood => (:char, :gas, :vapour);
        A = 1e8,
        E = 1.2e5,
        h = 0.0,
        yields = (0.25, 0.65, 0.1),
    )
    mat = _wet((_moisture(),); reactions = (wet_pyrolysis,))
    err = _refusal(() -> dry_gas_flags(mat))
    @test err isa ArgumentError
    @test occursin("gas :vapour is produced from the solids (reactions :pyro)", err.msg)
    @test occursin("vapour of phase change :moisture", err.msg)
    @test occursin("split it into one species per origin", err.msg)
    # without the phase change the same reaction is dry gas
    @test dry_gas_flags(_wet((); reactions = (wet_pyrolysis,)))[4]
    # the refusal reaches a solve through the lateral shrinkage law
    shrinking = Material(
        name = :shrinking,
        components = COMPONENTS,
        reactions = (wet_pyrolysis,),
        phase_changes = (_moisture(),),
        lateral_shrinkage_law = χ -> 1 - 0.2χ,
    )
    @test_throws ArgumentError Pyrolysis.Closures.LateralShrinkage(shrinking, 1.0, 1.0)
end
