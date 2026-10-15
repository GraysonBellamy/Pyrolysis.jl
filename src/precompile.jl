"""
Precompile workload — solve the configurations most runs start from while the
package image is built, so a session's first solve of such a problem reuses
cached native code instead of compiling the solver.

A solve compiles two layers. The integrator (OrdinaryDiffEq's stepper, its
nonlinear solver and the factorization) compiles for the solver configuration
only: the mesh mode, depletion handling, the integrator and the factorization
`solve` picks. `solve` hands it problem-independent function types
(`ErasedRHS`, `ErasedJacobian`), so the problem's types do not enter. The
Pyrolysis kernels (initialization, residual, Jacobian, solution extraction)
compile for the types of the problem: the component, property, reaction and
boundary-condition types, the component count, and the radiation and transport
modes. A workload therefore caches the integrator for every problem solved in
its configuration, and the kernels for the problems whose types it shares. The
workloads use the constructors' plain forms (constant properties, one
decomposition reaction) with the cone boundary conditions (incident radiation,
re-radiation and convection on the exposed face, an insulated back face):

| Workload | Components | Mesh | Radiation | Transport | Jacobian, factorization |
|---|---|---|---|---|---|
| `fixed_surface_2` | 2 | fixed | surface absorption | Fickian | default (`DirectSolve`, KLU) |
| `ale_surface_2` | 2 | ALE | surface absorption | Fickian | default (`ApproxSolve`, KLU) |
| `ale_surface_2_exact` | 2 | ALE | surface absorption | Fickian | `DirectSolve`, UMFPACK |
| `ale_surface_2_depletion` | 2 | ALE, parking | surface absorption | Fickian | default (`ApproxSolve`, KLU) |
| `fixed_surface_2_depletion` | 2 | fixed, parking | surface absorption | Fickian | default (`DirectSolve`, KLU) |
| `fixed_beer_lambert_darcy_2` | 2 | fixed | Beer–Lambert | Darcy | default (`DirectSolve`, KLU) |
| `ale_beer_lambert_darcy_2` | 2 | ALE | Beer–Lambert | Darcy | default (`ApproxSolve`, KLU) |
| `fixed_surface_3` | 3 | fixed | surface absorption | Fickian | default (`DirectSolve`, KLU) |
| `ale_beer_lambert_darcy_3` | 3 | ALE | Beer–Lambert | Darcy | default (`ApproxSolve`, KLU) |

Every workload uses the default integrator (`KenCarp4`). Its meshes have eight
cells: the exact ALE Jacobian takes the chunk width of its geometry
derivatives from the node count, up to 8, and a smaller mesh would compile a
width that no production mesh uses. A failing workload warns and the others
still run; the configuration it covers then compiles on first use.
"""

using PrecompileTools: @setup_workload, @compile_workload

const _WORKLOAD_CELLS = 8

# Constant-property materials in the constructors' plain forms. A permeable
# solid (`darcy`) selects Darcy transport.
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
                α = 2.4,
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

function _workload_problem(material, ::Val{NC}; L, q, darcy, t_end) where {NC}
    return PyrolysisProblem(
        mesh = create_uniform_mesh(L, _WORKLOAD_CELLS, Val(NC)),
        material = material,
        bc_set = _workload_bcs(q; darcy),
        tspan = (0.0, t_end),
        dt_initial = 1e-3,
        output_times = collect(range(0.0, t_end, length = 5)),
    )
end

"""
	_precompile_workloads() -> Vector{Pair{String,Function}}

The precompile workloads, by name. Each function runs one solve and returns its
`PyrolysisSolution`; the package's precompile step runs them all, and the
tests check that each one succeeds.
"""
function _precompile_workloads()
    two(; darcy = false, q = 50e3, L = 0.01, t_end = 20.0) =
        _workload_problem(_workload_material_2(; darcy), Val(2); L, q, darcy, t_end)
    three(; darcy = false) = _workload_problem(
        _workload_material_3(; darcy),
        Val(3);
        L = 0.01,
        q = 50e3,
        darcy,
        t_end = 20.0,
    )
    quiet = (progress = false,)
    return [
        "fixed_surface_2" =>
            () -> solve(two(); quiet..., radiation_model = SURFACE_ABSORPTION),
        "ale_surface_2" =>
            () -> solve(
                two();
                quiet...,
                use_ale = true,
                radiation_model = SURFACE_ABSORPTION,
            ),
        "ale_surface_2_exact" =>
            () -> solve(
                two();
                quiet...,
                use_ale = true,
                radiation_model = SURFACE_ABSORPTION,
                jacobian = Structured(strategy = DirectSolve()),
            ),
        # A 1 mm slab at 100 kW/m²: its surface cells burn down and park
        # within the horizon.
        "ale_surface_2_depletion" =>
            () -> solve(
                two(; q = 100e3, L = 1e-3, t_end = 40.0);
                quiet...,
                use_ale = true,
                radiation_model = SURFACE_ABSORPTION,
                handle_depletion = true,
                min_thickness = 1e-4,
            ),
        # The same slab on the fixed mesh: its top cells gasify and park.
        "fixed_surface_2_depletion" =>
            () -> solve(
                two(; q = 100e3, L = 1e-3, t_end = 40.0);
                quiet...,
                radiation_model = SURFACE_ABSORPTION,
                handle_depletion = true,
            ),
        "fixed_beer_lambert_darcy_2" =>
            () -> solve(two(; darcy = true); quiet..., radiation_model = BEER_LAMBERT),
        "ale_beer_lambert_darcy_2" =>
            () -> solve(
                two(; darcy = true);
                quiet...,
                use_ale = true,
                radiation_model = BEER_LAMBERT,
            ),
        "fixed_surface_3" =>
            () -> solve(three(); quiet..., radiation_model = SURFACE_ABSORPTION),
        "ale_beer_lambert_darcy_3" =>
            () -> solve(
                three(; darcy = true);
                quiet...,
                use_ale = true,
                radiation_model = BEER_LAMBERT,
            ),
    ]
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
