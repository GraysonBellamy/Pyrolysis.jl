"""
A solve ends where a live cell reaches the end of its declared transport
domain: located at the crossing, with return code `:Terminated`, the outcome of
the transport model and a warning naming the cell.

The deck is a slab whose solid converts entirely to gas, under lumped
transport: its top cell burns through the matrix. Without an endpoint the
domain ends where the matrix occupancy reaches `blend[2]` (`:matrix_exhausted`);
with one, where the cell holds more skeleton than the endpoint supports
(`:unsupported_endpoint_occupancy`). A deactivation threshold above
`blend[2]` takes the cell before its matrix runs out, and a residue keeps the
deck inside the domain.
"""

using Test
using Logging
using Pyrolysis
using Pyrolysis.Materials: R_GAS
using Pyrolysis.Integration: setup, integrate, DomainCheck, domain_margin, forget!
using Pyrolysis.Outputs: ConservationRecorder, conservation_audit
using Pyrolysis.States: mass_index

const BLEND_HI = 1e-3

function _charfree(; char = 0.0, endpoint = nothing)
    solid = SolidComponent(
        :solid;
        φ = 0.3,
        λ = 1e-5,
        κ = 1e-12,
        ρ = 1000.0,
        c = 1000.0,
        k = 0.2,
        ε = 0.9,
    )
    vol = GaseousComponent(:volatiles; M = 0.070, c = 2000.0, k = 0.03)
    n2 = GaseousComponent(:N2; M = 0.028014, c = 1100.0, k = 0.04)
    pore_n2 = z -> 0.3 * 101325.0 * 0.028014 / (R_GAS * 360.0)
    components, reaction, ξ = if char > 0
        residue = SolidComponent(
            :char;
            φ = 0.3,
            λ = 1e-5,
            κ = 1e-12,
            ρ = 200.0,
            c = 1000.0,
            k = 0.1,
        )
        (
            (solid, residue, vol, n2),
            Reaction(
                :decomp,
                1 => (2, 3);
                A = 1e5,
                E = 5e4,
                h = 1e5,
                yields = (char, 1 - char),
            ),
            [z -> 1000.0, z -> 0.0, z -> 0.0, pore_n2],
        )
    else
        (
            (solid, vol, n2),
            Reaction(:decomp, 1 => 2; A = 1e5, E = 5e4, h = 1e5),
            [z -> 1000.0, z -> 0.0, pore_n2],
        )
    end
    material = Material(
        name = :charfree,
        components = components,
        reactions = (reaction,),
        gas_transport = LumpedTransfer(; endpoint, blend = (1e-4, BLEND_HI)),
    )
    purge = GasReservoir(P = 101325.0, T = 360.0, mole_fractions = (N2 = 1.0,))
    return PyrolysisProblem(
        mesh = uniform_mesh(0.005, 8),
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = HeatFluxBC(5e4) + ConvectiveBC(h = 11.2, T_inf = 360.0),
                bottom = AdiabaticBC(),
            ),
            mass = (top = AmbientGasBC(purge), bottom = ImpermeableBC()),
        ),
        T_initial = z -> 360.0,
        ξ_initial = ξ,
        tspan = (0.0, 20.0),
    )
end

# The top cell's matrix occupancy at the last saved state.
_top_occupancy(sol) = sol.ξ[1, 8, end] / 1000.0

const _END_WARNING =
    r"solve\(\): the solve ends at t = .*, where cell 8 leaves its transport domain \(matrix_exhausted\)"

@testset "the matrix runs out: the solve ends at the crossing" begin
    sol =
        @test_logs (:warn, _END_WARNING) min_level = Logging.Warn Pyrolysis.Integration.solve(
            _charfree();
            saveat = [0.0, 20.0],
        )
    @test sol.retcode === :Terminated
    @test sol.outcome === :matrix_exhausted
    @test !succeeded(sol)
    @test !completed(sol)
    @test 10.0 < sol.t[end] < 20.0
    # The saved state lies at the crossing, on the side inside the domain.
    @test BLEND_HI <= _top_occupancy(sol) <= BLEND_HI * (1 + 1e-9)
    # The end is an event, not a run of rejected proposals.
    @test sol.statistics.domain_rejections == 0

    # The located end converges with the tolerances.
    fine = with_logger(NullLogger()) do
        Pyrolysis.Integration.solve(_charfree(); saveat = [0.0, 20.0], reltol = 1e-8)
    end
    @test fine.outcome === :matrix_exhausted
    @test abs(fine.t[end] - sol.t[end]) < 1e-4
end

@testset "the end is located on fixed steps" begin
    sol =
        @test_logs (:warn, _END_WARNING) min_level = Logging.Warn Pyrolysis.Integration.solve(
            _charfree();
            saveat = [0.0, 20.0],
            adaptive = false,
            dt = 1e-3,
        )
    @test sol.retcode === :Terminated
    @test sol.outcome === :matrix_exhausted
    @test BLEND_HI <= _top_occupancy(sol) <= BLEND_HI * (1 + 1e-9)
end

@testset "an endpoint ends where its supported occupancy does" begin
    # The default supports no skeleton, which a matrix that converts
    # continuously still holds where its occupancy reaches blend[2].
    empty_channel = LumpedEndpoint(λ = 1e-5, κ = 1e-12)
    sol = @test_logs (
        :warn,
        r"cell 8 leaves its transport domain \(unsupported_endpoint_occupancy\).*max_skeletal_occupancy = 0\.0",
    ) min_level = Logging.Warn Pyrolysis.Integration.solve(
        _charfree(endpoint = empty_channel);
        saveat = [0.0, 20.0],
    )
    @test sol.retcode === :Terminated
    @test sol.outcome === :unsupported_endpoint_occupancy

    covered = LumpedEndpoint(λ = 1e-5, κ = 1e-12, max_skeletal_occupancy = 1e-3)
    sol = @test_logs min_level = Logging.Warn Pyrolysis.Integration.solve(
        _charfree(endpoint = covered);
        saveat = [0.0, 20.0],
    )
    @test sol.retcode === :Success
    @test _top_occupancy(sol) < BLEND_HI
end

@testset "a threshold above blend[2] takes the cell first" begin
    sol = @test_logs min_level = Logging.Warn Pyrolysis.Integration.solve(
        _charfree();
        saveat = [0.0, 20.0],
        deactivation = true,
        threshold = 2 * BLEND_HI,
    )
    @test sol.retcode === :Success
    @test sol.outcome === :none
    @test length(sol.events) == 1
    @test only(sol.events).cell == 8
end

@testset "a residue keeps the deck inside the domain" begin
    sol = @test_logs min_level = Logging.Warn Pyrolysis.Integration.solve(
        _charfree(char = 0.2);
        saveat = [0.0, 20.0],
    )
    @test sol.retcode === :Success
    @test sol.outcome === :none
end

@testset "an initial state outside the domain is refused" begin
    cp = setup(_charfree())
    u0 = copy(cp.u0)
    u0[mass_index(cp.def.layout, 8, 1)] *= BLEND_HI / 2
    @test_throws ArgumentError integrate(cp; u0)
    sol = with_logger(NullLogger()) do
        integrate(cp; saveat = [0.0])
    end
    @test sol.outcome === :matrix_exhausted
end

@testset "the conservation ledgers close at the domain's end" begin
    # The recorders close the step the event cuts: the closure is the deck's
    # own at these tolerances (8e-7 energy, 1e-8 mass run to a time span just
    # short of the crossing).
    rec = ConservationRecorder()
    sol = with_logger(NullLogger()) do
        Pyrolysis.Integration.solve(
            _charfree();
            saveat = [0.0, 20.0],
            callback = rec,
            reltol = 1e-8,
        )
    end
    @test sol.outcome === :matrix_exhausted
    audit = conservation_audit(sol, rec)
    @test abs(audit.energy.relative) <= 2e-6
    @test abs(audit.mass.relative) <= 2e-8
end

@testset "a cell exactly at the domain's end still sees the next step leave" begin
    cp = setup(_charfree())
    check = DomainCheck(cp)
    k = mass_index(cp.def.layout, 8, 1)
    u0 = copy(cp.u0)
    # The top cell's solid at exactly blend[2] of its bulk density: a zero margin.
    m = u0[k] * BLEND_HI
    for _ = 1:1000
        u0[k] = m
        forget!(check)
        margin = first(domain_margin(check, u0))
        margin == 0 && break
        m = margin > 0 ? prevfloat(m) : nextfloat(m)
    end
    forget!(check)
    @test first(domain_margin(check, u0)) == 0
    sol = with_logger(NullLogger()) do
        integrate(cp; u0, saveat = [0.0, 20.0])
    end
    @test sol.retcode === :Terminated
    @test sol.outcome === :matrix_exhausted
    @test sol.t[end] < 1.0
    @test _top_occupancy(sol) >= BLEND_HI * (1 - 1e-12)
end
