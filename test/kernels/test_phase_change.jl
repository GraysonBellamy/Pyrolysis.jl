"""The phase-change kernel: its sources are the pure functions' at the derived
state, it skips dead cells, and a material without phase changes has none."""

using Test
using Pyrolysis
using Pyrolysis.Materials: n_components
using Pyrolysis.Physics: moisture_content, phase_change_rate, phase_change_heat
using Pyrolysis.States: derive!
using Pyrolysis.Kernels:
    PhaseChangeKernel, ReactionKernel, Contributions, contribute!, reset!, find_kernel
using Pyrolysis.Integration: setup

# A charring wood holding water and an oil, each with its vapour, in air: two
# phase changes of different laws beside a reaction.
function _wet(; phase_changes = true)
    wood = SolidComponent(
        :wood;
        ρ = 400.0,
        ρ_intrinsic = 1500.0,
        c = 1500.0,
        k = 0.12,
        λ = 5e-7,
    )
    char = SolidComponent(
        :char;
        ρ = 100.0,
        ρ_intrinsic = 1800.0,
        c = 1100.0,
        k = 0.08,
        λ = 1e-6,
    )
    water = LiquidComponent(:water; ρ = 1000.0, c = 4180.0, k = 0.6)
    oil = LiquidComponent(:oil; ρ = 900.0, c = 2000.0, k = 0.15)
    vapour = GaseousComponent(:vapour; M = 0.018015, c = 1900.0, k = 0.025)
    oil_vapour = GaseousComponent(:oil_vapour; M = 0.15, c = 1800.0, k = 0.02)
    gas = GaseousComponent(:gas; M = 0.03, c = 1100.0, k = 0.03)
    air = GaseousComponent(:air; M = 0.02897, c = 1006.0, k = 0.026)
    changes = (
        PhaseChange(
            :moisture,
            :water => :vapour;
            activity = CAST_SPRUCE,
            sorbents = (:wood, :char => 0.5),
        ),
        PhaseChange(
            :oil,
            :oil => :oil_vapour;
            activity = ActivityFunction((X, T) -> 2X * exp(400 * (1 / 300 - 1 / T))),
            sorbents = (:wood,),
            saturation = T -> 50.0 * exp(0.04 * (T - 300)),
            latent_heat = 3e5,
            exchange_rate = 5e2,
        ),
    )
    return Material(
        name = :wet,
        components = (wood, char, water, oil, vapour, oil_vapour, gas, air),
        reactions = (
            Reaction(
                :pyro,
                :wood => (:char, :gas);
                A = 5e9,
                E = 1.4e5,
                h = 8e5,
                yields = (0.3, 0.7),
            ),
        ),
        phase_changes = phase_changes ? changes : (),
    )
end

# A slab of six cells, temperature and moisture varying through it, the
# vapours away from equilibrium in both directions.
function _slab(mat)
    NC = n_components(mat)
    ξ = [
        z -> 400.0 - 2e4 * z,
        z -> 5.0 + 1e3 * z,
        z -> 10.0 + 4e3 * z,
        z -> 8.0,
        z -> 0.02 + 30 * z,
        z -> 0.01,
        z -> 0.05,
        z -> 0.8,
    ]
    return PyrolysisProblem(
        mesh = uniform_mesh(0.01, 6),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(1e4), bottom = AdiabaticBC()),
            mass = (top = ImpermeableBC(), bottom = ImpermeableBC()),
        ),
        T_initial = z -> 330.0 + 4e3 * z,
        ξ_initial = ξ[1:NC],
        tspan = (0.0, 1.0),
    )
end

_kernel_contribution(c, kernel, ds, topo, mat) =
    (reset!(c); contribute!(c, kernel, ds, topo, mat, 0.0))
_allocated(c, kernel, ds, topo, mat) = @allocated contribute!(c, kernel, ds, topo, mat, 0.0)

@testset "phase-change kernel: sources equal the pure functions" begin
    mat = _wet()
    cp = setup(_slab(mat))
    kernels = cp.def.kernels
    kernel = find_kernel(PhaseChangeKernel, kernels)
    @test kernel isa PhaseChangeKernel
    # after the reaction kernel, first among the others
    @test kernels[1] isa ReactionKernel && kernels[2] === kernel
    @test count(k -> k isa PhaseChangeKernel, kernels) == 1
    ds = cp.ws.derived
    topo = cp.def.setup.topology
    derive!(ds, cp.u0, cp.def.layout, cp.def.setup)
    NC = n_components(mat)
    c = Contributions{Float64}(NC, 1, 6, 7, 7)
    _kernel_contribution(c, kernel, ds, topo, mat)
    for i = 1:6
        ξ = ntuple(j -> ds.ξ[j, i], NC)
        ξ⁺ = ntuple(j -> ds.ξ⁺[j, i], NC)
        V = ds.V_eval[i]
        S_mass = zeros(NC)
        S_energy = 0.0
        for pc in mat.phase_changes
            X = moisture_content(pc, ξ, ξ⁺)
            ω = phase_change_rate(pc, X, ds.T_trial[i], ds.φ[i], ξ[pc.vapour])
            h = phase_change_heat(
                pc,
                X,
                ds.T_trial[i],
                ds.h[pc.vapour, i],
                ds.h[pc.liquid, i],
            )
            @test ω != 0
            S_mass[pc.liquid] -= ω * V
            S_mass[pc.vapour] += ω * V
            S_energy -= h * ω * V
            # each liquid and its vapour equal and opposite
            @test c.S_mass[pc.liquid, i] == -c.S_mass[pc.vapour, i]
        end
        for j = 1:NC
            @test c.S_mass[j, i] ≈ S_mass[j] rtol = 1e-14
        end
        # every other row untouched
        for j in (1, 2, 7, 8)
            @test c.S_mass[j, i] == 0
        end
        @test c.S_energy[i] ≈ S_energy rtol = 1e-14
    end
    # both directions are present: evaporation and condensation
    @test any(>(0), c.S_mass[5, :]) && any(<(0), c.S_mass[5, :])
    # the moisture content reads the sorbents' positive parts with their weights
    pc = first(mat.phase_changes)
    @test moisture_content(pc, (400.0, 10.0, 40.0), (400.0, 10.0, 40.0)) ==
          40.0 / (405.0 + 1e-3)
    # nothing else of the buffer is written
    @test all(iszero, c.F_mass) && all(iszero, c.F_energy) && all(iszero, c.θ)

    @testset "a dead cell contributes nothing" begin
        ds.active[3] = false
        _kernel_contribution(c, kernel, ds, topo, mat)
        @test all(iszero, c.S_mass[:, 3]) && c.S_energy[3] == 0
        @test all(i -> c.S_energy[i] != 0, (1, 2, 4, 5, 6))
        ds.active[3] = true
    end

    @testset "zero allocations" begin
        _allocated(c, kernel, ds, topo, mat)
        @test _allocated(c, kernel, ds, topo, mat) == 0
    end

    @testset "the reaction kernel's sources are its own" begin
        # the full assembly's sources are the two kernels' sums
        cr = Contributions{Float64}(NC, 1, 6, 7, 7)
        reset!(cr)
        contribute!(cr, kernels[1], ds, topo, mat, 0.0)
        dry_cp = setup(_slab(_wet(phase_changes = false)))
        @test find_kernel(PhaseChangeKernel, dry_cp.def.kernels) === nothing
        @test count(k -> k isa ReactionKernel, dry_cp.def.kernels) == 1
        cd = Contributions{Float64}(NC, 1, 6, 7, 7)
        dds = dry_cp.ws.derived
        derive!(dds, dry_cp.u0, dry_cp.def.layout, dry_cp.def.setup)
        reset!(cd)
        contribute!(cd, dry_cp.def.kernels[1], dds, topo, dry_cp.def.setup.material, 0.0)
        @test cr.S_mass == cd.S_mass
        @test cr.S_energy == cd.S_energy
    end
end

@testset "a material without phase changes has no phase-change kernel" begin
    for mat in (
        _wet(phase_changes = false),
        Material(
            name = :inert,
            components = (SolidComponent(:s; ρ = 1000.0, c = 1000.0, k = 1.0),),
        ),
    )
        cp = setup(
            PyrolysisProblem(
                mesh = uniform_mesh(0.01, 4),
                material = mat,
                bc_set = BoundaryConditionSet(
                    thermal = (top = HeatFluxBC(1e4), bottom = AdiabaticBC()),
                ),
                T_initial = z -> 330.0,
                ξ_initial = n_components(mat) == 1 ? [z -> 1000.0] : _slab(mat).ξ_initial,
                tspan = (0.0, 1.0),
            ),
        )
        @test find_kernel(PhaseChangeKernel, cp.def.kernels) === nothing
    end
end
