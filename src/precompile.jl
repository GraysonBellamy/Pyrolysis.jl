"""
Precompile workload — solve the configurations most runs start from while the
package image is built, so a session's first solve of such a problem reuses
cached native code instead of compiling the solver.

A solve compiles two layers. The integrator (OrdinaryDiffEq's stepper, its
linear solve and the callbacks) compiles for the solver configuration only:
whether the mesh moves, whether in-depth radiation adds algebraic rows on a
fixed mesh, and whether cells are deactivated. `solve` hands it
problem-independent function and callback types (`ErasedRHS`,
`ErasedJacobian`), so the problem's types do not enter. The kernels, the state
mapping, the colored Jacobian and the solution's reconstruction compile for the
types of the problem: the component, property, reaction and boundary-condition
types, the component count, the radiation model, the gas transport and the
energy state. A workload therefore caches the integrator for every problem
solved in its configuration, and the kernels for the problems whose types it
shares.

The workloads cover fixed and moving meshes × surface absorption and
Beer–Lambert × Fickian and Darcy–Fick transport × with and without
deactivation for a two-component material, each on its own name
(`fixed_surface_fickian_2`, …, `moving_beer_lambert_darcy_deactivation_2`); a
three-component material on a fixed mesh (`fixed_surface_fickian_3`) and on a
moving mesh (`moving_beer_lambert_darcy_3`);
and the enthalpy energy state once (`moving_surface_fickian_deactivation_2_enthalpy`).
They use the constructors' plain forms (constant properties, one
decomposition reaction) with the cone boundary conditions (incident radiation,
re-radiation and convection on the exposed face, an insulated back face) on
eight cells. A deck with deactivation is a 1 mm slab at 100 kW/m² over 80 s,
whose cells deactivate within the horizon. A failing workload warns and the others still
run; the configuration it covers then compiles on first use.
"""

using PrecompileTools: @setup_workload, @compile_workload

const _WORKLOAD_CELLS = 8

# Constant-property materials in the constructors' plain forms. A permeable
# solid (`darcy`) selects Darcy–Fick transport. The two-component solid absorbs
# strongly, so that under Beer–Lambert absorption its residue keeps absorbing
# the beam and the 1 mm decks burn out within their horizon.
function _workload_material_2(; darcy::Bool)
    return Material(
        name = :workload_2,
        components = (
            SolidComponent(
                :solid,
                ρ = 1190.0,
                c = 1420.0,
                k = 0.21,
                ε = 0.86,
                α = 200.0,
                κ = darcy ? 1e-14 : 0.0,
                ρ_intrinsic = 1500.0,
                λ = 1e-6,
            ),
            GaseousComponent(:gas, M = 0.1, c = 1500.0, k = 0.02, λ = 1e-5),
        ),
        reactions = (Reaction(:pyrolysis, 1 => 2, A = 8.5e12, E = 188e3, h = 870e3),),
    )
end

function _workload_material_3(; darcy::Bool)
    return Material(
        name = :workload_3,
        components = (
            SolidComponent(
                :virgin,
                ρ = 500.0,
                c = 1500.0,
                k = 0.15,
                ε = 0.9,
                α = 5.0,
                κ = darcy ? 1e-14 : 0.0,
                ρ_intrinsic = 1500.0,
                λ = 1e-5,
            ),
            SolidComponent(
                :char,
                ρ = 150.0,
                c = 1100.0,
                k = 0.08,
                ε = 0.95,
                α = 5.0,
                κ = darcy ? 1e-12 : 0.0,
                ρ_intrinsic = 1500.0,
                λ = 1e-5,
            ),
            GaseousComponent(:gas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5),
        ),
        reactions = (
            Reaction(
                :pyrolysis,
                1 => (2, 3);
                A = 1e10,
                E = 150e3,
                h = 300e3,
                yields = (0.25, 0.75),
            ),
        ),
    )
end

# The cone boundary conditions; `darcy` adds the ambient pressure at the
# exposed face.
function _workload_bcs(q; darcy::Bool)
    thermal = (
        top = RadiativeFluxBC(flux = q) +
              RadiativeBC(T_inf = 300.0) +
              ConvectiveBC(h = 10.0, T_inf = 300.0),
        bottom = AdiabaticBC(),
    )
    top_mass =
        darcy ? DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0) :
        DiffusiveMassBC(0.05)
    return BoundaryConditionSet(;
        thermal,
        mass = (top = top_mass, bottom = ImpermeableBC()),
    )
end

function _workload_problem(material; L, q, darcy, t_end)
    return PyrolysisProblem(
        mesh = create_uniform_mesh(L, _WORKLOAD_CELLS),
        material = material,
        bc_set = _workload_bcs(q; darcy),
        tspan = (0.0, t_end),
        dt_initial = 1e-3,
        output_times = collect(range(0.0, t_end, length = 5)),
    )
end

# A 1 cm slab at 50 kW/m² over 20 s, or, with deactivation, a 1 mm slab at
# 100 kW/m² over 80 s, whose cells deactivate within the horizon.
function _workload_deck(material; darcy::Bool, deactivation::Bool)
    deactivation &&
        return _workload_problem(material; L = 1e-3, q = 100e3, darcy, t_end = 80.0)
    return _workload_problem(material; L = 0.01, q = 50e3, darcy, t_end = 20.0)
end

"""
	_precompile_workloads() -> Vector{Pair{String,Function}}

The precompile workloads, by name. Each function runs one solve and returns its
`PyrolysisSolution`; the package's precompile step runs them all, and the
tests check that each one succeeds.
"""
function _precompile_workloads()
    workloads = Pair{String,Function}[]
    function add!(
        material,
        nc;
        moving,
        radiation,
        darcy,
        deactivation,
        energy = :temperature,
    )
        name = join(
            (
                moving ? "moving" : "fixed",
                radiation === BEER_LAMBERT ? "beer_lambert" : "surface",
                darcy ? "darcy" : "fickian",
                (deactivation ? ("deactivation",) : ())...,
                string(nc),
                (energy === :enthalpy ? ("enthalpy",) : ())...,
            ),
            "_",
        )
        run() = solve(
            _workload_deck(material(; darcy); darcy, deactivation);
            use_ale = moving,
            radiation_model = radiation,
            deactivation,
            energy_state = energy,
        )
        push!(workloads, name => run)
        return nothing
    end
    for moving in (false, true),
        radiation in (SURFACE_ABSORPTION, BEER_LAMBERT),
        darcy in (false, true),
        deactivation in (false, true)

        add!(_workload_material_2, 2; moving, radiation, darcy, deactivation)
    end
    add!(
        _workload_material_3,
        3;
        moving = false,
        radiation = SURFACE_ABSORPTION,
        darcy = false,
        deactivation = false,
    )
    add!(
        _workload_material_3,
        3;
        moving = true,
        radiation = BEER_LAMBERT,
        darcy = true,
        deactivation = false,
    )
    add!(
        _workload_material_2,
        2;
        moving = true,
        radiation = SURFACE_ABSORPTION,
        darcy = false,
        deactivation = true,
        energy = :enthalpy,
    )
    return workloads
end

@setup_workload begin
    workloads = _precompile_workloads()
    @compile_workload begin
        for (name, run) in workloads
            try
                Base.CoreLogging.with_logger(run, Base.CoreLogging.NullLogger())
            catch err
                @warn "Pyrolysis precompile workload `$name` failed; its configuration " *
                      "compiles on first use instead" exception = (err, catch_backtrace())
            end
        end
    end
end
