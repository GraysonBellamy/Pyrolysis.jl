"""
Save the results of a Pyrolysis.jl run to CSV files, which open directly in Excel.

    julia --project=examples examples/export_results_csv.jl

`solve` keeps its results in memory and writes no files. `write_results_csv` below
writes two files from a solution. Copy the function into your own script and call it
after `solve`:

    write_results_csv(solution, "my_run")

my_run_time.csv: one row per saved time
    t_s               time [s]
    T_surface_K       exposed-face temperature [K]
    T_back_K          back-face temperature [K]
    mass_kg_m2        sample mass per initial area [kg/m²]
    MLR_kg_m2_s       mass-loss rate per initial area [kg/(m²·s)], as a cone calorimeter
                      reports it (`solution.extras.MLR_gauge`)
    thickness_m       sample thickness [m]

my_run_profiles.csv: one row per cell per saved time, from the back face to the exposed face
    t_s               time [s]
    z_m               cell-centre height above the back face [m]
    T_K               cell temperature [K]
    xi_<name>_kg_m3   concentration of each component, named as in the Material [kg/m³]

Rows are written at the saved times: the problem's `output_times`, or `saveat` passed
to `solve`. Set those to the times you want in the files. Cells removed by
surface deactivation are left out of the profile file.

The case below is a 5.8 mm PMMA slab under 50 kW/m² for 300 s, with the boundary
conditions of `06_cone_calorimeter_surface_absorption.jl` and constant properties.
The files go to `output/` at the repository root.
"""

using Pyrolysis
using DelimitedFiles   # Julia standard library, nothing to install

const CSV_SMOKE = get(ENV, "PYROLYSIS_SMOKE", "0")=="1"
const PMMA_POROSITY = parse(Float64, get(ENV, "PYROLYSIS_PMMA_POROSITY", "0.01"))

function write_results_csv(solution, prefix)
    t = solution.t
    time_table = hcat(
        t,
        solution.surface_T,
        solution.extras.back_temperature,
        solution.extras.total_mass,
        solution.extras.MLR_gauge,
        solution.extras.thickness_history,
    )
    open("$(prefix)_time.csv", "w") do io
        println(io, "t_s,T_surface_K,T_back_K,mass_kg_m2,MLR_kg_m2_s,thickness_m")
        writedlm(io, time_table, ',')
    end

    # Saved cell centres include the current mesh motion and deactivation.
    names = [c.name for c in solution.problem.material.components]
    rows = Vector{Float64}[]
    for i in eachindex(t)
        z = get_cell_positions(solution; at = i)
        for j in axes(solution.T, 1)
            isnan(solution.T[j, i]) && continue   # deactivated cell
            push!(rows, [t[i]; z[j]; solution.T[j, i]; solution.ξ[:, j, i]])
        end
    end
    open("$(prefix)_profiles.csv", "w") do io
        println(io, join(["t_s", "z_m", "T_K", ["xi_$(n)_kg_m3" for n in names]...], ','))
        writedlm(io, permutedims(reduce(hcat, rows)), ',')
    end
    return nothing
end

pmma = Material(
    name = :PMMA,
    components = (
        SolidComponent(
            :virgin,
            φ = PMMA_POROSITY,
            λ = 1e-5,
            ρ = 1210.0,
            c = 1500.0,
            k = 0.2,
            ε = 0.96,
        ),
        SolidComponent(
            :char,
            φ = PMMA_POROSITY,
            λ = 1e-5,
            ρ = 1210.0,
            c = 1500.0,
            k = 0.2,
            ε = 0.96,
        ),
        GaseousComponent(:MMA, M = 0.10012, c = 1500.0, k = 0.02),
        GaseousComponent(:AIR, M = 0.02897, c = 1006.0, k = 0.026),
    ),
    reactions = (
        Reaction(
            :decomposition,
            1 => (2, 3);
            A = 1.5e14,
            E = 2.03e5,
            h = 820e3,
            yields = (0.002, 0.998),
        ),
    ),
)

mesh = create_uniform_mesh(5.8e-3, CSV_SMOKE ? 12 : 60)
air_initial =
    background_gas_concentration(pmma, 4, (1210.0, 0.0, 0.0, 0.0), 300.0; P = 101325.0)
t_final = CSV_SMOKE ? 20.0 : 300.0

bc_set = BoundaryConditionSet(
    thermal = (
        top = RadiativeFluxBC(flux = 50e3, absorptivity = 0.96) +
              RadiativeBC(ε = 0.96, T_inf = 300.0) +
              ConvectiveBC(h = 10.0, T_inf = 300.0),
        bottom = RadiativeBC(ε = 0.95, T_inf = 300.0) +
                 ConvectiveBC(h = 4.0, T_inf = 300.0),
    ),
    mass = (
        top = AmbientGasBC(
            GasReservoir(P = 101325.0, T = 300.0, mole_fractions = (AIR = 1.0,));
            film = MassTransferFilm(0.05*101325*0.02897/(8.31446261815324*300)),
        ),
        bottom = ImpermeableBC(),
    ),
)

problem = PyrolysisProblem(
    mesh = mesh,
    material = pmma,
    bc_set = bc_set,
    T_initial = z -> 300.0,
    ξ_initial = [z -> 1210.0, z -> 0.0, z -> 0.0, z -> air_initial],
    tspan = (0.0, t_final),
    output_times = collect(0.0:2.0:t_final),   # the rows of the CSV files
)

solution = solve(
    problem;
    progress = false,
    use_ale = true,
    reltol = 1e-9,
    abstol = 1e-11,
    pressure_tolerance = 1e-4,
    radiation_model = SURFACE_ABSORPTION,
    deactivation = true,
    min_thickness = 0.01e-3,
)
println("solve: ", solution.retcode, ", ", length(solution.t), " saved times")

output_dir = joinpath(@__DIR__, "..", "output")
mkpath(output_dir)
prefix = joinpath(output_dir, "pmma_50kW")
write_results_csv(solution, prefix)
println("wrote $(prefix)_time.csv and $(prefix)_profiles.csv")
