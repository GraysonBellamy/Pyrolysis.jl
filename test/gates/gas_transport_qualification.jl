# Joint nominal-pore-volume, mesh and pressure-error qualification of the
# shipped dense-polymer decks. Run with --project=test. Each CSV row records
# the complete solve work and conservation, not just an initial EOS check.
using Pyrolysis, Test, Printf
const ROOT=normpath(joinpath(@__DIR__, "../.."))
const DEST=get(ENV, "PYROLYSIS_QUALIFICATION_DIR", "/tmp/pyrolysis-gas-qualification")
mkpath(DEST)

# Evaluate the shipped example's definitions through problem construction.
# Plotting and file output follow its solve and are not part of this gate.
function example_problem(file)
    text=first(split(read(joinpath(ROOT, "examples", file), String), "solution = solve"))
    text=replace(text, r"^using Plots\s*$"m=>"", r"^using DelimitedFiles\s*$"m=>"")
    mod=Module(gensym(:DenseExample))
    Core.eval(mod, :(using Pyrolysis))
    redirect_stdout(devnull) do
        Base.include_string(mod, text, file)
    end
    return Base.invokelatest(getproperty, mod, :problem)
end
function pore_component(c, φ)
    is_gas(c) && return c
    # The example's condensed components all have constant bulk density.
    SolidComponent(
        c.name;
        ρ = c.density,
        ρ_intrinsic = c.density(300.0)/(1-φ),
        c = c.heat_capacity,
        k = c.thermal_conductivity,
        λ = c.gas_transfer_coeff,
        κ = c.permeability,
        ε = c.emissivity,
        absorptivity = c.surface_absorptivity,
        α = c.absorption_coeff,
        γ = c.swelling_factor,
    )
end
function dense_variant(base, φ, n; t_end = base.tspan[2])
    m=base.material
    mat=Material(
        name = m.name,
        components = map(c->pore_component(c, φ), m.components),
        reactions = m.reactions,
        mixing = MaterialMixing(; m.mixing...),
        gas_transport = m.gas_transport,
        gas_viscosity = m.gas_viscosity,
        depletion_limiter = m.depletion_limiter,
        lateral_shrinkage_law = m.lateral_shrinkage_law,
    )
    air=findfirst(c->is_gas(c) && c.name===:AIR, mat.components)
    fields=map(eachindex(mat.components)) do j
        j==air ? (z->begin
            ξ=Tuple(f(z) for f in base.ξ_initial)
            background_gas_concentration(mat, air, ξ, base.T_initial(z))
        end) : base.ξ_initial[j]
    end
    L=base.mesh.positions[end][1]
    return PyrolysisProblem(
        mesh = create_uniform_mesh(L, n),
        material = mat,
        bc_set = base.bc_set,
        T_initial = base.T_initial,
        ξ_initial = fields,
        tspan = (0.0, t_end),
        dt_initial = base.dt_initial,
        output_times = collect(range(0.0, t_end; length = 201)),
    )
end
function qualify_dense(
    file;
    meshes = (24, 48, 72),
    pores = (0.03, 0.02, 0.01),
    ptols = (0.0001, 0.00001),
    energy = :temperature,
    reltol = 1e-9,
    abstol = 1e-11,
    label = "",
)
    base=example_problem(file)
    return Base.invokelatest(
        _qualify_dense,
        base,
        file;
        meshes,
        pores,
        ptols,
        energy,
        reltol,
        abstol,
        label,
    )
end
function _qualify_dense(base, file; meshes, pores, ptols, energy, reltol, abstol, label)
    path=joinpath(DEST, replace(file, ".jl"=>"_$(energy)$(label).csv"))
    measurements=Dict{Tuple{Int,Float64,Float64},NTuple{3,Float64}}()
    open(path, "w") do io
        println(
            io,
            "porosity,cells,pressure_tolerance,reltol,abstol,retcode,outcome,time_end,peak_temperature,peak_pressure,product_escape,air_ingress,maximum_species_relative_defect,energy_relative_defect,accepted,rejected,rhs,jacobians,linear_solves,nonlinear_failures,domain_rejections,max_surface_iterations,wall_seconds",
        )
        @testset "dense polymer $file ($energy)" begin
            for n in meshes, φ in pores, ptol in ptols
                p=dense_variant(base, φ, n)
                rec=ConservationRecorder()
                elapsed=@elapsed sol=solve(
                    p;
                    energy_state = energy,
                    use_ale = true,
                    deactivation = true,
                    radiation_model = SURFACE_ABSORPTION,
                    min_thickness = 1e-5,
                    reltol,
                    abstol,
                    pressure_tolerance = ptol,
                    dtmin = 1e-14,
                    callback = rec,
                )
                @test completed(sol)
                species=gas_species_audit(sol, rec)
                defect=maximum(abs(s.relative) for s in values(species.species))
                energy_defect=conservation_audit(sol, rec).energy.relative
                println((ledger = species, energy = energy_defect))
                flush(stdout)
                finite_max(x) = maximum(filter(isfinite, vec(x)))
                product=species.species.MMA.escape
                ingress=species.species.AIR.ingress
                stats=sol.statistics
                row=(
                    φ,
                    n,
                    ptol,
                    reltol,
                    abstol,
                    sol.retcode,
                    sol.outcome,
                    sol.t[end],
                    finite_max(sol.T),
                    finite_max(get_pressure(sol)),
                    product,
                    ingress,
                    defect,
                    energy_defect,
                    stats.accepted,
                    stats.rejected,
                    stats.rhs_evaluations,
                    stats.jacobians,
                    stats.linear_solves,
                    stats.nonlinear_failures,
                    stats.domain_rejections,
                    maximum(get_boundary_exchange(sol, :top).iterations),
                    elapsed,
                )
                println(io, join(row, ','))
                flush(io)
                println((file, energy, row))
                flush(stdout)
                measurements[(n, φ, ptol)]=(
                    finite_max(sol.T),
                    finite_max(get_pressure(sol)),
                    product,
                )
                @test defect<1e-4
                @test abs(energy_defect)<1e-4
            end
            # Tolerances are stated in physical units, separately from the
            # ledger checks. Nominal storage must not change the thermal or
            # escaped-product solution at the resolution of this qualification.
            for n in meshes, ptol in ptols
                q=[measurements[(n, φ, ptol)] for φ in pores]
                @test maximum(first, q)-minimum(first, q)<0.01 # K
                @test maximum(x->x[2], q)-minimum(x->x[2], q)<100 # Pa
                @test maximum(last, q)-minimum(last, q)<1e-4*maximum(last, q)
            end
            for n in meshes, φ in pores
                a, b=measurements[(n, φ, first(ptols))], measurements[(n, φ, last(ptols))]
                @test abs(a[1]-b[1])<0.001
                @test abs(a[3]-b[3])<1e-5*b[3]
            end
            if length(meshes)>=3
                for φ in pores
                    a, b, c=(measurements[(n, φ, last(ptols))] for n in meshes[(end-2):end])
                    @test abs(c[1]-b[1])<0.8*abs(b[1]-a[1])
                    @test abs(c[3]-b[3])<0.8*abs(b[3]-a[3])
                end
            end
        end
    end
    return path
end
if abspath(PROGRAM_FILE)==@__FILE__
    files=isempty(ARGS) ?
          [
        "06_cone_calorimeter.jl",
        "06_cone_calorimeter_surface_absorption.jl",
        "export_results_csv.jl",
    ] : ARGS
    for file in files, energy in (:temperature, :enthalpy)
        qualify_dense(file; energy)
    end

end
