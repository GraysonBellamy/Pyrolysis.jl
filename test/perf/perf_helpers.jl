"""
Helpers shared across the performance-invariant tests.

Builds canonical configurations for `SimulationMode{...}` and returns
`(def, ws, u, du)` tuples ready to feed `residual!` into.
"""

if !isdefined(@__MODULE__, :_PYROLYSIS_PERF_HELPERS_LOADED)
    const _PYROLYSIS_PERF_HELPERS_LOADED = true

    using Pyrolysis
    using Pyrolysis.Materials: SURFACE_ABSORPTION, BEER_LAMBERT

    const PERF_NC = 2
    const PERF_N_CELLS = 100

    # ---------------------------------------------------------------------------
    # Internal: assemble the (def, ws, u, du) tuple from any built problem.
    # ---------------------------------------------------------------------------
    function _finalize(
        problem;
        use_ale,
        radiation_model = Pyrolysis.Materials.NO_RADIATION,
        radiation_intensity::Real = 0.0,
    )
        def, ws, u = initialize_state(
            problem;
            use_ale = use_ale,
            radiation_model = radiation_model,
            radiation_intensity = radiation_intensity,
        )
        du = similar(u)
        return def, ws, u, du
    end

    # ---------------------------------------------------------------------------
    # Canonical materials
    # ---------------------------------------------------------------------------

    # 2-component PMMA-like (Fickian gas transport).
    function _pmma_material()
        virgin = SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86)
        gas = GaseousComponent(:MMA, M = 0.100, c = 1500, k = 0.02, λ = 1e-5)
        rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
        return Material(name = :PMMA, components = (virgin, gas), reactions = (rxn,))
    end

    # 2-component charring solid with permeability + absorption (Beer-Lambert + Darcy).
    # The default α = 1000 m²/kg makes the 1 cm column optically thick (τ ≈ 5000);
    # κ = 0 selects Fickian gas transport instead of Darcy.
    function _wood_like_material(; α::Real = 1000.0, κ::Real = 1e-14)
        virgin = SolidComponent(
            :wood,
            ρ = 500.0,
            c = 1500.0,
            k = 0.15,
            ε = 0.9,
            α = α,             # m²/kg — absorption coefficient for Beer-Lambert
            κ = κ,             # m² — permeability for Darcy flow
            λ = 1e-5,          # mass diffusivity (matrix transport coefficient)
        )
        gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5)
        rxn = Reaction(:pyrolysis, 1 => 2, A = 1e10, E = 150e3, h = -300e3)
        return Material(name = :wood, components = (virgin, gas), reactions = (rxn,))
    end

    # ---------------------------------------------------------------------------
    # Canonical configurations — one per `SimulationMode{...}` axis combination
    # ---------------------------------------------------------------------------

    """
    	build_baseline_eulerian()

    `SimulationMode{Eulerian, NoRadiation, Fickian, NoChi}` — PMMA-like, 100 cells.
    Used by the zero-allocation and JET tests.
    """
    function build_baseline_eulerian()
        mat = _pmma_material()
        mesh = create_uniform_mesh(0.01, PERF_N_CELLS, Val(PERF_NC))
        bc =
            BoundaryConditionSet(thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()))
        problem =
            PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
        return _finalize(problem; use_ale = false)
    end

    """
    	build_ale_no_radiation()

    `SimulationMode{ALE, NoRadiation, Fickian, NoChi}` — PMMA-like with shrinking mesh.
    """
    function build_ale_no_radiation()
        mat = _pmma_material()
        mesh = create_uniform_mesh(0.01, PERF_N_CELLS, Val(PERF_NC))
        bc =
            BoundaryConditionSet(thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()))
        problem =
            PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
        return _finalize(problem; use_ale = true)
    end

    """
    	build_eulerian_surface_absorption()

    `SimulationMode{Eulerian, SurfaceAbsorption, Fickian, NoChi}` — radiation
    deposited at the surface cell only.
    """
    function build_eulerian_surface_absorption()
        mat = _pmma_material()
        mesh = create_uniform_mesh(0.01, PERF_N_CELLS, Val(PERF_NC))
        bc = BoundaryConditionSet(thermal = (top = HeatFluxBC(0.0), bottom = AdiabaticBC()))
        problem =
            PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
        return _finalize(
            problem;
            use_ale = false,
            radiation_model = SURFACE_ABSORPTION,
            radiation_intensity = 50e3,
        )
    end

    """
    	build_baseline_eulerian_darcy()

    `SimulationMode{Eulerian, NoRadiation, DarcyFick, NoChi}` — wood-like material
    (component κ > 0 triggers DarcyFick) on a non-shrinking mesh. Used by the
    `DarcyTransportOperator` correctness tests.
    """
    function build_baseline_eulerian_darcy()
        mat = _wood_like_material()
        mesh = create_uniform_mesh(0.01, PERF_N_CELLS, Val(PERF_NC))
        bc =
            BoundaryConditionSet(thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()))
        problem =
            PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
        return _finalize(problem; use_ale = false)
    end

    """
    	build_ale_beer_lambert_darcy()

    `SimulationMode{ALE, BeerLambert, DarcyFick, NoChi}` — wood-like material
    with permeability (→ Darcy flow) and a non-zero absorption coefficient
    (→ Beer-Lambert volumetric absorption).
    """
    function build_ale_beer_lambert_darcy()
        mat = _wood_like_material()
        mesh = create_uniform_mesh(0.01, PERF_N_CELLS, Val(PERF_NC))
        bc = BoundaryConditionSet(thermal = (top = HeatFluxBC(0.0), bottom = AdiabaticBC()))
        problem =
            PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
        return _finalize(
            problem;
            use_ale = true,
            radiation_model = BEER_LAMBERT,
            radiation_intensity = 50e3,
        )
    end

    # Optically moderate ALE + Beer-Lambert deck: α = 0.5 m²/kg gives a total
    # optical depth of about 2.5 over the 1 cm column, so the intensity reaching
    # every cell depends measurably on the node positions of all cells above it.
    function _build_ale_beer_lambert_moderate(κ::Real)
        mat = _wood_like_material(α = 0.5, κ = κ)
        mesh = create_uniform_mesh(0.01, PERF_N_CELLS, Val(PERF_NC))
        bc = BoundaryConditionSet(thermal = (top = HeatFluxBC(0.0), bottom = AdiabaticBC()))
        problem =
            PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
        return _finalize(
            problem;
            use_ale = true,
            radiation_model = BEER_LAMBERT,
            radiation_intensity = 50e3,
        )
    end

    """
    	build_ale_beer_lambert_darcy_moderate()

    `SimulationMode{ALE, BeerLambert, DarcyFick, NoChi}` on an optically moderate
    medium (total optical depth ≈ 2.5). Unlike the optically thick
    [`build_ale_beer_lambert_darcy`](@ref), the energy rows' dependence on the node
    positions far above each cell does not underflow, so a sparsity pattern that
    omits those z columns fails the strict comparison.
    """
    build_ale_beer_lambert_darcy_moderate() = _build_ale_beer_lambert_moderate(1e-14)

    """
    	build_ale_beer_lambert_fickian_moderate()

    `SimulationMode{ALE, BeerLambert, Fickian, NoChi}` counterpart of
    [`build_ale_beer_lambert_darcy_moderate`](@ref) (zero permeability).
    """
    build_ale_beer_lambert_fickian_moderate() = _build_ale_beer_lambert_moderate(0.0)

    """
    	mixed_reacting_state!(u, n_cells)

    Overwrite the T and first two ξ blocks of `u` with a reacting profile:
    T from 450 to 650 K, the solid density reduced by 5–35%, and a gas
    concentration from 0.02 to 0.1 kg/m³. Node positions are left unchanged.
    """
    function mixed_reacting_state!(u, n_cells)
        for i = 1:n_cells
            u[i] = 450.0 + 200.0 * i / n_cells
            u[n_cells+i] *= 0.65 + 0.3 * i / n_cells
            u[2n_cells+i] = 0.02 + 0.08 * i / n_cells
        end
        return u
    end

end # _PYROLYSIS_PERF_HELPERS_LOADED
