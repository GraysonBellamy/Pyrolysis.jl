# Run with --project=benchmark. Cost is measured only after warm-up. Each model
# is compared with its own tight trajectory against the same physical targets;
# the two models need not predict the same trajectory.
using Pyrolysis, BenchmarkTools, Printf
using Pyrolysis.Materials: R_GAS
using Pyrolysis.Assembly: residual!, buffers
using Pyrolysis.Kernels:
    GasTransportKernel, ThermalBoundaryKernel, find_kernel, contribute!, reset!

function transport_benchmark_problem(::Val{N}, dusty) where {N}
    names=ntuple(j->Symbol(:gas, j), Val(N))
    gases=ntuple(
        j->GaseousComponent(names[j]; M = 0.02+0.003j, c = 1000.0 + 20j, k = 0.03),
        Val(N),
    )
    model=dusty ?
          DustyGas(
        tortuosity = 2.0,
        D_ref = 2e-5,
        T_ref = 300.0,
        P_ref = 1e5,
        pore_diameter = 1e-6,
    ) : LumpedTransfer()
    mat=Material(
        name = :benchmark_transport,
        gas_transport = model,
        components = (
            SolidComponent(
                :matrix;
                ρ = 500.0,
                φ = 0.5,
                c = 1000.0,
                k = 0.2,
                λ = 1e-5,
                κ = 1e-14,
            ),
            gases...,
        ),
    )
    x=NamedTuple{names}(ntuple(j->N==1 ? 1.0 : j==1 ? 0.7 : 0.3/(N-1), Val(N)))
    reservoir=GasReservoir(P = 1e5, T = 300.0, mole_fractions = x)
    fields=(
        z->500.0,
        ntuple(j->z->0.5*(1.2e5+1e5z)*(0.02+0.003j)/(N*R_GAS*(500.0 + 1e4z)), Val(N))...,
    )
    PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, 8),
        material = mat,
        T_initial = z->500.0 + 1e4z,
        ξ_initial = collect(fields),
        bc_set = BoundaryConditionSet(
            thermal = (
                top = ConvectiveBC(h = 20.0, T_inf = 300.0, blowing = CouetteFilm()),
                bottom = AdiabaticBC(),
            ),
            mass = (
                top = AmbientGasBC(reservoir; film = MassTransferFilm(0.001)),
                bottom = ImpermeableBC(),
            ),
        ),
        tspan = (0.0, 0.2),
        output_times = collect(0.0:0.02:0.2),
    )
end
function kernel_cost(cp, K)
    d, w=cp.def, cp.ws
    du=zero(cp.u0)
    residual!(du, cp.u0, d, 0.0, w)
    ds, c=buffers(w, d, Float64)
    kernel=find_kernel(K, d.kernels)
    step=()->begin
        reset!(c)
        contribute!(c, kernel, ds, d.setup.topology, d.setup.material, 0.0)
        nothing
    end
    step()
    trial=@benchmark $step() seconds=0.3 samples=100 evals=1
    (;
        nanoseconds = median(trial).time,
        bytes = median(trial).memory,
        allocations = median(trial).allocs,
    )
end
function benchmark_transport(N, dusty, io)
    p=transport_benchmark_problem(Val(N), dusty)
    cp=setup(p; pressure_tolerance = 0.001)
    interior=kernel_cost(cp, GasTransportKernel)
    boundary=kernel_cost(cp, ThermalBoundaryKernel)
    tight=solve(p; reltol = 1e-10, abstol = 1e-12, pressure_tolerance = 0.0001)
    completed(tight) || error("tight reference did not complete")
    # An actual warm solve precedes timing; no compilation enters the work figure.
    options=(reltol = 1e-8, abstol = 1e-10, pressure_tolerance = 0.001)
    run=()->solve(p; options...)
    run()
    timing=@benchmark $run() samples=5 evals=1 seconds=10
    seconds=median(timing).time/1e9
    sol=run()
    completed(sol) || error("timed solve did not complete")
    p_error=maximum(abs.(get_pressure(sol) .- get_pressure(tight)))
    T_error=maximum(abs.(sol.T .- tight.T))
    x_error=maximum(abs.(get_gas_composition(sol) .- get_gas_composition(tight)))
    p_error<1.0 && T_error<1e-4 && x_error<1e-6 ||
        error("unmatched solution accuracy: $((p_error,T_error,x_error))")
    interior.bytes==0 && boundary.bytes==0 || error("warmed kernel allocations")
    st=sol.statistics
    row=(
        dusty ? "dusty" : "lumped",
        N,
        interior.nanoseconds,
        interior.bytes,
        boundary.nanoseconds,
        boundary.bytes,
        seconds,
        p_error,
        T_error,
        x_error,
        st.accepted,
        st.rejected,
        st.rhs_evaluations,
        st.jacobians,
        st.linear_solves,
        st.nonlinear_failures,
        maximum(get_boundary_exchange(sol, :top).iterations),
    )
    println(io, join(row, ','))
    flush(io)
    println(row)
    flush(stdout)
end
if abspath(PROGRAM_FILE)==@__FILE__
    output=isempty(ARGS) ? "/tmp/pyrolysis-gas-performance.csv" : ARGS[1]
    open(output, "w") do io
        println(
            io,
            "model,species,interior_ns,interior_bytes,film_boundary_ns,film_boundary_bytes,solve_seconds,max_pressure_error_Pa,max_temperature_error_K,max_mole_fraction_error,accepted,rejected,rhs,jacobians,linear_solves,nonlinear_failures,max_saved_film_iterations",
        )
        for n in (1, 2, 4, 8), dusty in (false, true)
            benchmark_transport(n, dusty, io)
        end
    end
end
