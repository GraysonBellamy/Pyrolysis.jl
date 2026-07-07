using Pyrolysis

if !isdefined(@__MODULE__, :fixture_problem)
    function fixture_material(; lateral::Bool = false, darcy::Bool = false)
        solid_kwargs = darcy ? (; κ = 1e-14) : (;)
        gas_kwargs = darcy ? (; κ = 0.0) : (;)
        solid = SolidComponent(
            :solid;
            ρ = 1000.0,
            c = 1000.0,
            k = 0.2,
            ε = 0.9,
            γ = 1.0,
            solid_kwargs...,
        )
        gas = GaseousComponent(
            :gas;
            M = 0.050,
            c = 1000.0,
            k = 0.02,
            λ = 1e-5,
            γ = 1.0,
            gas_kwargs...,
        )
        rxn = Reaction(:decomp, 1 => 2, A = 1e5, E = 5e4, h = 1e5)
        law = lateral ? (χ -> 1.0 - 0.2 * χ) : nothing
        return Material(
            name = :fixture_test,
            components = (solid, gas),
            reactions = (rxn,),
            lateral_shrinkage_law = law,
        )
    end

    function fixture_problem(;
        n_cells::Int = 6,
        lateral::Bool = false,
        darcy::Bool = false,
        tspan = (0.0, 0.02),
        T_initial::Real = 300.0,
    )
        material = fixture_material(lateral = lateral, darcy = darcy)
        mesh = create_uniform_mesh(
            0.01,
            n_cells,
            Val(2);
            T_initial = Float64(T_initial),
            ξ_initial = (1000.0, 0.0),
        )
        bcs = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(1e4), bottom = AdiabaticBC()),
            mass = (top = DiffusiveMassBC(0.01), bottom = ImpermeableBC()),
        )
        problem = PyrolysisProblem(mesh, material, bcs, tspan)
        Pyrolysis.Problem.apply_initial_conditions!(problem)
        return problem
    end

    function fixture_state(
        problem;
        use_ale::Bool = false,
        radiation_model = NO_RADIATION,
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
end
