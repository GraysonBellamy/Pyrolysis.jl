"""
Precompile workload — warm up the most common Material × `SimulationMode`
combinations so first-call latency drops from minutes to seconds.

Covers the four canonical configurations:

1. PMMA cone-25, Eulerian, surface absorption
2. PMMA cone-25, ALE, surface absorption
3. Charring material (white-pine-like), ALE, Beer-Lambert + Darcy-Fick
4. Generic 2-component thermoset, ALE

Each workload runs a tiny mesh / short tspan solve to trigger compilation
of the residual, Jacobian assembly, integrator step, and linear-solve
paths. The work is run inside `@compile_workload`, so the resulting native
code is cached in the package image.

Errors during the workload are non-fatal: PrecompileTools rolls back its
recording on exception, the package still loads, and only the missed
specializations are compiled lazily on first user call.
"""

using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    # Lightweight fixtures shared across the four workloads. Allocations
    # inside `@setup_workload` are not recorded, so this is free.
    _pmma_cp = ConstantProperty(1500.0)
    _pmma_k = ConstantProperty(0.21)

    _pmma_components = (
        SolidComponent(
            :virgin,
            ρ = 1190.0,
            c = _pmma_cp,
            k = _pmma_k,
            ε = 0.96,
            α = 0.0,
            λ = 2e-5,
        ),
        SolidComponent(
            :char,
            ρ = 1190.0,
            c = _pmma_cp,
            k = _pmma_k,
            ε = 0.96,
            α = 0.0,
            λ = 2e-5,
        ),
        GaseousComponent(:MMA, M = 0.10012, c = 1500.0, k = 0.02, λ = 1e-5),
    )
    _pmma_reactions = (
        Reaction(
            :decomposition,
            1 => (2, 3),
            A = 1.5e14,
            E = 2.03e5,
            h = 820e3,
            n = 1.0,
            yields = (0.002, 0.998),
        ),
    )
    _pmma_material = Material(
        name = :PMMA_pre,
        components = _pmma_components,
        reactions = _pmma_reactions,
        mixing = MaterialMixing(k = MixingSpec(WEIGHTED, 0.5), λ = MixingSpec(PARALLEL)),
    )

    # Charring material (white-pine-like) with Darcy permeability. The
    # skeletal densities (ρ_intrinsic) are required for a permeable material:
    # without them porosity is exactly 0, which problem validation rejects.
    _wood_components = (
        SolidComponent(
            :virgin,
            ρ = 500.0,
            c = 2000.0,
            k = 0.15,
            ε = 0.9,
            α = 5.0,
            κ = 1e-15,
            ρ_intrinsic = 1500.0,
        ),
        SolidComponent(
            :char,
            ρ = 150.0,
            c = 1500.0,
            k = 0.10,
            ε = 0.95,
            α = 50.0,
            κ = 1e-12,
            ρ_intrinsic = 1950.0,
        ),
        GaseousComponent(:gas, M = 0.030, c = 1500.0, k = 0.04, λ = 1e-5),
    )
    _wood_reactions = (
        Reaction(
            :pyrolysis,
            1 => (2, 3),
            A = 1e10,
            E = 1.5e5,
            h = 200e3,
            n = 1.0,
            yields = (0.25, 0.75),
        ),
    )
    _wood_material = Material(
        name = :wood_pre,
        components = _wood_components,
        reactions = _wood_reactions,
        mixing = MaterialMixing(
            k = MixingSpec(WEIGHTED, 0.5),
            λ = MixingSpec(PARALLEL),
            κ = MixingSpec(PARALLEL),
        ),
    )

    # Generic 2-component thermoset.
    _ts_components = (
        SolidComponent(:resin, ρ = 1300.0, c = 1700.0, k = 0.20, ε = 0.9, λ = 2e-5),
        GaseousComponent(:vapor, M = 0.080, c = 1400.0, k = 0.03, λ = 1e-5),
    )
    _ts_reactions = (Reaction(:cure, 1 => 2, A = 1e12, E = 1.6e5, h = 300e3, n = 1.0),)
    _ts_material = Material(
        name = :thermoset_pre,
        components = _ts_components,
        reactions = _ts_reactions,
    )

    # Tiny BCs reused across modes.
    _top_thermal = HeatFluxBC(25e3) + RadiativeBC(ε = 0.9, T_inf = 300.0)
    _bot_thermal = AdiabaticBC()
    _bcs_eulerian =
        BoundaryConditionSet(thermal = (top = _top_thermal, bottom = _bot_thermal))
    _bcs_ale = BoundaryConditionSet(
        thermal = (top = _top_thermal, bottom = _bot_thermal),
        mass = (top = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0), bottom = ImpermeableBC()),
    )

    @compile_workload begin
        try
            # (1) PMMA Eulerian — heat-flux BC stands in for surface absorption
            # (the SURFACE_ABSORPTION trait deposits radiation in the surface
            # cell; HeatFluxBC drives the same code path here).
            let mesh = create_uniform_mesh(
                    5e-3,
                    5,
                    Val(3);
                    T_initial = 300.0,
                    ξ_initial = (1190.0, 0.0, 0.0),
                )
                bcs = _bcs_eulerian
                prob = PyrolysisProblem(
                    mesh = mesh,
                    material = _pmma_material,
                    bc_set = bcs,
                    tspan = (0.0, 0.5),
                    dt_initial = 1e-3,
                    output_times = [0.0, 0.5],
                )
                solve(
                    prob;
                    progress = false,
                    use_ale = false,
                    radiation_model = SURFACE_ABSORPTION,
                )
            end

            # (2) PMMA ALE, surface absorption.
            let mesh = create_uniform_mesh(
                    5e-3,
                    5,
                    Val(3);
                    T_initial = 300.0,
                    ξ_initial = (1190.0, 0.0, 0.0),
                )
                prob = PyrolysisProblem(
                    mesh = mesh,
                    material = _pmma_material,
                    bc_set = _bcs_ale,
                    tspan = (0.0, 0.5),
                    dt_initial = 1e-3,
                    output_times = [0.0, 0.5],
                )
                solve(
                    prob;
                    progress = false,
                    use_ale = true,
                    radiation_model = SURFACE_ABSORPTION,
                )
            end

            # (3) Charring material, ALE, Beer-Lambert + Darcy-Fick.
            let mesh = create_uniform_mesh(
                    5e-3,
                    5,
                    Val(3);
                    T_initial = 300.0,
                    ξ_initial = (500.0, 0.0, 0.0),
                )
                prob = PyrolysisProblem(
                    mesh = mesh,
                    material = _wood_material,
                    bc_set = _bcs_ale,
                    tspan = (0.0, 0.5),
                    dt_initial = 1e-3,
                    output_times = [0.0, 0.5],
                )
                solve(
                    prob;
                    progress = false,
                    use_ale = true,
                    radiation_model = BEER_LAMBERT,
                )
            end

            # (4) Generic 2-component thermoset, ALE.
            let mesh = create_uniform_mesh(
                    5e-3,
                    5,
                    Val(2);
                    T_initial = 300.0,
                    ξ_initial = (1300.0, 0.0),
                )
                prob = PyrolysisProblem(
                    mesh = mesh,
                    material = _ts_material,
                    bc_set = _bcs_ale,
                    tspan = (0.0, 0.5),
                    dt_initial = 1e-3,
                    output_times = [0.0, 0.5],
                )
                solve(prob; progress = false, use_ale = true)
            end
        catch err
            # Swallow precompile-workload failures — the package still loads
            # and missed specializations compile lazily on first call.
            @debug "Pyrolysis precompile workload failed; lazy compile fallback active" exception=(
                err,
                catch_backtrace(),
            )
        end
    end
end
