module SolverProbe

using LinearAlgebra, SparseArrays
using OrdinaryDiffEq, OrdinaryDiffEqSDIRK
import LinearSolve, SciMLBase

export Ladder,
    exact_state,
    recurrence,
    ladder_rhs!,
    ladder_jac!,
    make_problem,
    run_ladder,
    constraint_error,
    solution_error,
    reduced_jacobian,
    augmented_matrix,
    BorderCache,
    border_solve!,
    solve_mask_counterexample

# A manufactured stiff diffusion problem with a nonlinear prefix coupling:
# y' = L*y + beta*z + f(t), z_i - z_(i-1) = y_i^2/n, z_0=0.
# Lifted rows are local with interleaved ordering; eliminating z makes J dense
# below its first upper diagonal. Both forms retain precisely the same coupling.
mutable struct Ladder
    n::Int
    diffusion::Float64
    beta::Float64
    after_jump::Bool
    jumps::Int
    jac_derivative::Vector{Float64}
    jac_revision::Int
end
Ladder(n) = Ladder(n, 0.25 * (n + 1)^2, 0.4, false, 0, zeros(n), 0)
recurrence(y) = cumsum(y .^ 2) / length(y)

function exact_state(p, t; after = p.after_jump)
    x = (1:p.n) ./ (p.n + 1)
    slow = 0.2 .* sinpi.(x) .* exp(-t)
    fast = 0.1 .* sinpi.(p.n .* x) .* exp(-200t)
    y, dy = 1 .+ slow .+ fast, -slow .- 200 .* fast
    if after
        restart = 0.05 .* sinpi.(p.n .* x) .* exp(-200(t - 0.5))
        y = 0.6 .* y .+ restart
        dy = 0.6 .* dy .- 200 .* restart
    end
    return y, dy
end

function ladder_rhs!(du, u, p, t, lifted)
    n = p.n
    y = lifted ? view(u, 1:2:2n) : u
    z = lifted ? view(u, 2:2:2n) : recurrence(y)
    reference, dy = exact_state(p, t)
    zr = recurrence(reference)
    for i = 1:n
        v = y[i] - reference[i]
        left = i > 1 ? y[i-1] - reference[i-1] : zero(v)
        right = i < n ? y[i+1] - reference[i+1] : zero(v)
        value = p.diffusion * (left - 2v + right) + p.beta * (z[i] - zr[i]) + dy[i]
        if lifted
            du[2i-1] = value
            du[2i] = z[i] - (i > 1 ? z[i-1] : zero(v)) - y[i]^2 / n
        else
            du[i] = value
        end
    end
    return nothing
end

function reduced_jacobian(y, p)
    n = p.n
    J = Matrix(
        Tridiagonal(fill(p.diffusion, n-1), fill(-2p.diffusion, n), fill(p.diffusion, n-1)),
    )
    for j = 1:n, i = j:n
        J[i, j] += p.beta * 2y[j] / n
    end
    return J
end

function ladder_jac!(J, u, p, t, lifted)
    fill!(J, 0)
    n = p.n
    y = lifted ? view(u, 1:2:2n) : u
    p.jac_derivative .= 2 .* y ./ n
    p.jac_revision += 1
    if lifted
        for i = 1:n
            J[2i-1, 2i-1] = -2p.diffusion
            J[2i-1, 2i] = p.beta
            i > 1 && (J[2i-1, 2i-3] = p.diffusion)
            i < n && (J[2i-1, 2i+1] = p.diffusion)
            J[2i, 2i-1] = -p.jac_derivative[i]
            J[2i, 2i] = 1
            i > 1 && (J[2i, 2i-2] = -1)
        end
    else
        J .= reduced_jacobian(y, p)
    end
    return nothing
end

# Exact augmented solve for (mu*I + factor*J_reduced)*delta = b:
# [mu*I + factor*L   factor*beta*I] [delta] = [b]
# [-diag(2y/n)       B           ] [eta  ]   [0]
# where B is the unit lower-bidiagonal prefix recurrence. Interleaving gives
# a sparse matrix with O(n) nonzeros. This is an eliminated ODE with an exact
# augmented Newton solve, NOT another time-integrated DAE or an approximate J.
function augmented_matrix(p, mu, factor, derivative)
    rows, cols, vals = Int[], Int[], Float64[]
    function add(i, j, v)
        push!(rows, i)
        push!(cols, j)
        push!(vals, v)
    end
    for i = 1:p.n
        add(2i-1, 2i-1, mu - 2factor * p.diffusion)
        add(2i-1, 2i, factor * p.beta)
        i > 1 && add(2i-1, 2i-3, factor * p.diffusion)
        i < p.n && add(2i-1, 2i+1, factor * p.diffusion)
        add(2i, 2i-1, -derivative[i])
        add(2i, 2i, 1.0)
        i > 1 && add(2i, 2i-2, -1.0)
    end
    return sparse(rows, cols, vals, 2p.n, 2p.n)
end

mutable struct BorderCache
    factorization::Any
    key::Tuple{Float64,Float64,Int}
    rhs::Vector{Float64}
    calls::Int
    factorizations::Int
    verify::Bool
    max_residual::Float64
end
BorderCache(n; verify = false) =
    BorderCache(nothing, (NaN, NaN, -1), zeros(2n), 0, 0, verify, 0.0)

function border_solve!(cache, A, b, u, p)
    # W is either J-I/gamma or I-gamma*J depending on the integrator. Recover
    # the two coefficients from entries without assuming a sign convention.
    factor = A[1, 2] / p.diffusion
    mu = A[1, 1] - factor * (-2p.diffusion + p.beta * p.jac_derivative[1])
    key = (Float64(mu), Float64(factor), p.jac_revision)
    if cache.factorization === nothing || key != cache.key
        cache.factorization = lu(augmented_matrix(p, mu, factor, p.jac_derivative))
        cache.key = key
        cache.factorizations += 1
    end
    fill!(cache.rhs, 0)
    cache.rhs[1:2:end] .= b
    delta = cache.factorization \ cache.rhs
    u .= view(delta, 1:2:length(delta))
    cache.calls += 1
    if cache.verify
        cache.max_residual =
            max(cache.max_residual, norm(A * u - b, Inf) / max(norm(b, Inf), eps()))
    end
    return u
end

function make_problem(n, form)
    p = Ladder(n)
    lifted = form == :lifted
    y, _ = exact_state(p, 0.0)
    initial = lifted ? vec(permutedims(hcat(y, recurrence(y)))) : y
    rhs!(du, u, p, t) = ladder_rhs!(du, u, p, t, lifted)
    jac!(J, u, p, t) = ladder_jac!(J, u, p, t, lifted)
    prototype = lifted ? spzeros(2n, 2n) : zeros(n, n)
    jac!(prototype, initial, p, 0.0)
    mass = lifted ? Diagonal(repeat([1.0, 0.0], n)) : I
    f = ODEFunction(rhs!; jac = jac!, jac_prototype = prototype, mass_matrix = mass)
    return ODEProblem(f, initial, (0.0, 1.0), p)
end

function run_ladder(n, method, form; tolerance = 1e-8, masked = false, verify = false)
    prob = make_problem(n, form)
    p, lifted = prob.p, form == :lifted
    cache = BorderCache(n; verify)
    custom(A, b, u, linear_parameters, fresh, Pl, Pr, val; kwargs...) =
        border_solve!(cache, A, b, u, p)
    linear =
        form == :bordered ? LinearSolve.LinearSolveFunction(custom) :
        lifted ? LinearSolve.UMFPACKFactorization() : LinearSolve.LHLFactorization()
    method in (:KenCarp4, :KenCarp4_raw, :FBDF, :Rodas5P) || error("unknown method")
    alg =
        method in (:KenCarp4, :KenCarp4_raw) ?
        KenCarp4(;
            linsolve = linear,
            concrete_jac = true,
            smooth_est = method == :KenCarp4,
        ) :
        method == :FBDF ? FBDF(; linsolve = linear, concrete_jac = true) :
        Rodas5P(; linsolve = linear, concrete_jac = true)
    function jump!(integrator)
        state = integrator.u
        y = lifted ? view(state, 1:2:2n) : state
        y .*= 0.6
        y .+= 0.05 .* sinpi.(n .* ((1:n) ./ (n + 1)))
        lifted && (state[2:2:2n] .= recurrence(y))
        p.after_jump = true
        p.jumps += 1
        derivative_discontinuity!(integrator, true)
    end
    callback = DiscreteCallback((u, t, i) -> t == 0.5, jump!; save_positions = (true, true))
    full_norm(u, t) = norm(u) / sqrt(length(u))
    masked_norm(u, t) = u isa Number ? abs(u) : norm(view(u, 1:2:length(u))) / sqrt(n)
    sol = solve(
        prob,
        alg;
        callback,
        tstops = [0.5],
        reltol = tolerance,
        abstol = tolerance * 0.01,
        initializealg = SciMLBase.NoInit(),
        internalnorm = masked ? masked_norm : full_norm,
        maxiters = 100000,
    )
    return (; sol, p, cache)
end

function constraint_error(sol, p)
    worst = 0.0
    for u in sol.u, i = 1:p.n
        z, previous, source = u[2i], i > 1 ? u[2i-2] : 0.0, u[2i-1]^2 / p.n
        scale = max(1 / p.n, abs(z) + abs(previous) + abs(source))
        worst = max(worst, abs(z - previous - source) / scale)
    end
    return worst
end

function solution_error(sol, p, lifted)
    worst = 0.0
    for (t, u) in zip(sol.t, sol.u)
        y = lifted ? view(u, 1:2:length(u)) : u
        before = exact_state(p, t; after = t > 0.5)[1]
        err = norm(y - before, Inf)
        # Callback intentionally saves both one-sided states at the jump.
        t == 0.5 && (err = min(err, norm(y - exact_state(p, t; after = true)[1], Inf)))
        worst = max(worst, err)
    end
    return worst
end

function solve_mask_counterexample(masked)
    function rhs!(du, u, p, t)
        du[1] = 0
        du[2] = u[2]^2 - (1 + t)
    end
    f = ODEFunction(rhs!; mass_matrix = Diagonal([1.0, 0.0]))
    full(u, t) = norm(u) / sqrt(length(u))
    mask(u, t) = abs(u[1])
    return solve(
        ODEProblem(f, [1.0, 1.0], (0.0, 1.0)),
        KenCarp4();
        adaptive = false,
        dt = 0.1,
        reltol = 1e-10,
        abstol = 1e-10,
        initializealg = SciMLBase.NoInit(),
        internalnorm = masked ? mask : full,
    )
end

end
