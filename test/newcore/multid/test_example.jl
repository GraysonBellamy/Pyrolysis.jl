using Test
using LinearAlgebra
using SparseArrays
using StaticArrays
using OrdinaryDiffEq: ODEFunction, ODEProblem
import OrdinaryDiffEq
using Pyrolysis
using Pyrolysis.Materials: BEER_LAMBERT
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.States: state_length
using Pyrolysis.Assembly: residual!
using Pyrolysis.Closures: PotentialClosure
using Pyrolysis.ColoredJacobian
using Pyrolysis.Integration: setup, rodas5p_policy, default_linear_solver

# The integrator the package runs: Rodas5P with the colored Jacobian and the
# time derivative by finite differences.
const _rodas5p = rodas5p_policy(default_linear_solver()).alg

# A char-forming solid whose gas does not move: the 2D path carries no gas
# transport.
function _wood(; λ = 0.0)
    virgin = SolidComponent(:wood, ρ = 500.0, c = 1500.0, k = 0.15, ε = 0.9)
    char = SolidComponent(:char, ρ = 150.0, c = 1100.0, k = 0.08, ε = 0.95)
    gas = GaseousComponent(:gas, M = 0.06, c = 1800.0, k = 0.03, λ = λ)
    rxn = Reaction(
        :pyrolysis,
        1 => (2, 3);
        A = 1e10,
        E = 1.5e5,
        h = -3e5,
        yields = (0.25, 0.75),
    )
    return Material(name = :wood, components = (virgin, char, gas), reactions = (rxn,))
end

_rot(θ) = @SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)]
_anisotropic(θ) = _rot(θ) * (@SMatrix [4.0 0.0; 0.0 1.0]) * _rot(θ)'

# A 2 cm × 1 cm slab, distorted smoothly inside a fixed outline, heated on top.
_inside(x::SVector{2}) = SVector(
    x[1] + 0.0015sin(π * x[1] / 0.02) * sin(2π * x[2] / 0.01),
    x[2] + 0.0006sin(π * x[2] / 0.01) * sin(2π * x[1] / 0.02),
)

function _example(;
    nx = 4,
    ny = 3,
    distorted = true,
    heated = true,
    T0 = x -> 600.0 + 2e4 * x[2],
    ξ0 = (x -> 450.0, x -> 15.0, x -> 0.5),
    material = _wood(),
)
    mesh = box_mesh(
        range(0.0, 0.02; length = nx + 1),
        range(0.0, 0.01; length = ny + 1);
        transform = distorted ? _inside : identity,
    )
    top = heated ? RadiativeFluxBC(flux = 5e4) + RadiativeBC(T_inf = 300.0) : AdiabaticBC()
    return PyrolysisProblem(
        mesh = mesh,
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = top,
                bottom = AdiabaticBC(),
                left = AdiabaticBC(),
                right = AdiabaticBC(),
            ),
        ),
        T_initial = T0,
        ξ_initial = collect(ξ0),
        tspan = (0.0, 1.0),
    )
end

_slip = (:left => 1, :right => 1)

# The pattern contains the dense Jacobian's nonzeros and the colored Jacobian
# equals it per entry, at a state where the residual is finite.
function _contained(cp, u; t = 0.5)
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    du = zeros(length(u))
    residual!(du, u, def, t, ws)
    @test all(isfinite, du)
    Jd = dense_jacobian(def, ws, u, t)
    @test all(isfinite, Jd)
    @test isempty(missing_entries(Jd, cj.pattern))
    J = jacobian!(jac_prototype(cj), cj, u, t)
    @test same_structure(J, cj.pattern)
    @test compare_entries(Matrix(J), Jd; rtol = 1e-12).passed
    return Jd
end

# Physical and trial states of a set-up example: a mixed reacting profile, a
# perturbation of every unknown (the node positions by a twentieth of a cell),
# a hot state, and negative masses and temperatures on the reference mesh.
function _states(cp)
    layout = cp.def.layout
    u0 = cp.u0
    n = n_cells(layout)
    D = dimension(layout)
    setup_ = cp.def.setup
    V = setup_.reference.volume
    cell(T, ξ) = begin
        u = copy(u0)
        for i = 1:n
            m = Tuple(ξ(i) .* V[i])
            for j in eachindex(m)
                u[mass_index(layout, i, j)] = m[j]
            end
            u[energy_index(layout, i)] =
                energy_state_value(energy_state(layout), setup_, T(i), m, i)
        end
        u
    end
    perturbed = u0 .* (1 .+ 0.2 .* sin.(1:length(u0)))
    if is_moving(layout)
        h = minimum(setup_.reference.thickness)
        for nd = 1:n_nodes(layout), k = 1:D
            r = position_index(layout, nd, k)
            perturbed[r] = u0[r] + 0.05h * sin(r)
        end
    end
    negative = cell(i -> -300.0 - i, i -> (-450.0, -15.0 * i, -0.5))
    return (;
        u0,
        mixed = cell(
            i -> 450.0 + 300.0 * i / n,
            i -> (450.0 * (0.4 + 0.5 * i / n), 20.0 + 60.0 * i / n, 0.3),
        ),
        hot = cell(i -> 900.0, i -> (300.0, 40.0, 1.0)),
        perturbed,
        negative,
    )
end

@testset "the 2D example sets up with its multipliers consistent: $(energy) form" for energy in
                                                                                      (
    :temperature,
    :enthalpy,
)
    cp = setup(
        _example();
        energy_state = energy,
        use_ale = true,
        slip = _slip,
        tensor = _anisotropic(π / 6),
    )
    layout = cp.def.layout
    @test dimension(layout) == 2
    @test closure_unknowns(layout) === :multipliers
    @test cp.def.closure isa PotentialClosure
    du = zeros(state_length(layout))
    residual!(du, cp.u0, cp.def, 0.0, cp.ws)
    θV = sum(abs, cp.ws.contributions.θ .* cp.ws.derived.V_eval)
    @test θV > 0
    @test maximum(abs, du[multiplier_index(layout, i)] for i = 1:n_cells(layout)) <=
          1e-10 * θV
    @test all(isfinite, du)
    @test all(>(0), cp.abstol)
end

@testset "the 2D example's refusals" begin
    p = _example()
    @test_throws ArgumentError setup(p; radiation_model = BEER_LAMBERT)
    @test_throws ArgumentError setup(_example(material = _wood(λ = 1e-5)))
    @test_throws ArgumentError setup(p; use_ale = true, closure = :balance)
    @test_throws ArgumentError setup(p; use_ale = true, w_override = (x, t) -> 0 * x)
    @test_throws ArgumentError setup(p; use_ale = true, closure = :lagrangian)
    @test_throws ArgumentError setup(
        p;
        use_ale = true,
        slip = _slip,
        tensor = -_anisotropic(0.3),
    )
end

@testset "containment of the 2D example: $(name), $(energy) form" for (name, options) in (
        ("potential closure", (use_ale = true, slip = _slip, tensor = _anisotropic(π / 6))),
        ("fixed mesh", (tensor = _anisotropic(π / 6),)),
    ),
    energy in (:temperature, :enthalpy)

    problem = _example(; T0 = x -> 600.0 + 2e4 * x[2] + 5e3 * x[1])
    if !get(options, :use_ale, false)
        m = problem.mesh
        cells = [m.topology.cells[i].nodes for i = 1:n_cells(m.topology)]
        problem = PyrolysisProblem(
            mesh = cell_mesh(
                m.positions,
                cells;
                tag = (i, k, s) -> ((:left, :right), (:bottom, :top))[k][s+1],
                fixed = true,
            ),
            material = problem.material,
            bc_set = problem.bc_set,
            T_initial = problem.T_initial,
            ξ_initial = problem.ξ_initial,
            tspan = problem.tspan,
        )
    end
    cp = setup(problem; energy_state = energy, options...)
    @info "2D example, $(name), $(energy) form: $(cp.jacobian.n_colors) colors, " *
          "$(nnz(cp.jacobian.pattern)) stored entries for $(state_length(cp.def.layout)) unknowns"
    for (state, u) in pairs(_states(cp))
        @testset "$state" begin
            Jd = _contained(cp, u)
            layout = cp.def.layout
            if state === :hot && is_moving(layout)
                # A multiplier row reads its cell's state through the strain
                # rate, and a position row the multipliers around its node.
                r = multiplier_index(layout, 5)
                @test Jd[r, energy_index(layout, 5)] != 0
                @test Jd[r, mass_index(layout, 5, 1)] != 0
                @test Jd[position_index(layout, 7, 2), multiplier_index(layout, 5)] != 0
            end
        end
    end
end

# The rows of a set-up problem's algebraic unknowns at the accepted states of
# a fixed-step Rodas5P integration, relative to the largest `θ V`.
function _algebraic_residuals(cp, dt; tend = 0.4)
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    layout = def.layout
    f = ODEFunction{true}(
        (du, u, p, t) -> residual!(du, u, def, t, ws);
        mass_matrix = mass_matrix(layout),
        jac = (J, u, p, t) -> (jacobian!(J, cj, u, t); nothing),
        jac_prototype = jac_prototype(cj),
    )
    sol = OrdinaryDiffEq.solve(
        ODEProblem(f, copy(cp.u0), (0.0, tend)),
        _rodas5p;
        dt,
        adaptive = false,
    )
    @test sol.retcode == OrdinaryDiffEq.ReturnCode.Success
    algebraic = findall(!, is_differential(layout))
    du = zeros(state_length(layout))
    worst = 0.0
    for (k, u) in enumerate(sol.u)
        residual!(du, u, def, sol.t[k], ws)
        c, ds = ws.contributions, ws.derived
        @test maximum(c.θ) - minimum(c.θ) <= 1e-12 * maximum(abs, c.θ)
        worst = max(worst, maximum(abs, du[algebraic]) / maximum(abs, c.θ .* ds.V_eval))
    end
    return worst, sol
end

# A uniform slab at a uniform temperature, insulated on every face: every cell
# reacts alike, so the strain rate is uniform, and the potential closure
# between slip walls lifts the slab uniformly. Its balances hold to round-off
# at the consistent state; at the accepted states of a Rosenbrock integration,
# which is linearly implicit, every algebraic row holds to the method's local
# error, as the line's balance rows do under the same steps.
@testset "the per-cell volume balance: round-off at the consistent state, the integrator's order at the accepted states: $(energy) form" for energy in
                                                                                                                                             (
    :temperature,
    :enthalpy,
)
    problem = _example(; distorted = false, heated = false, T0 = x -> 750.0)
    cp = setup(
        problem;
        energy_state = energy,
        use_ale = true,
        slip = _slip,
        abstol = 1e-10,
        reltol = 1e-10,
    )
    layout = cp.def.layout
    du = zeros(state_length(layout))
    residual!(du, cp.u0, cp.def, 0.0, cp.ws)
    c, ds = cp.ws.contributions, cp.ws.derived
    @test maximum(abs, c.r_multiplier) <= 1e-12 * maximum(abs, c.θ .* ds.V_eval)

    line = PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, 3, Val(3)),
        material = _wood(),
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
        ),
        T_initial = z -> 750.0,
        ξ_initial = [z -> 450.0, z -> 15.0, z -> 0.5],
        tspan = (0.0, 1.0),
    )
    cp1 = setup(line; energy_state = energy, use_ale = true, abstol = 1e-10, reltol = 1e-10)
    r2, r1 = Float64[], Float64[]
    for dt in (0.1, 0.05, 0.025)
        worst, sol = _algebraic_residuals(cp, dt)
        push!(r2, worst)
        push!(r1, first(_algebraic_residuals(cp1, dt)))
        @test sol.u[end][position_index(layout, n_nodes(layout), 2)] <
              cp.u0[position_index(layout, n_nodes(layout), 2)]
    end
    @info "largest per-cell balance residual at the accepted states, $(energy) form" r2 r1
    @test all(>(3.5), log2.(r2[1:(end-1)] ./ r2[2:end]))
    @test all(r2 .<= 5 .* r1)
end

# The largest number of stored entries in a row: the lower bound of any
# column coloring.
_largest_row(S) = maximum(diff(sparse(transpose(S)).colptr))

# The colors of the 2D example's Jacobian on meshes too large for every pair of
# unknowns to interact: the stencil's lower bound, the largest row, is constant
# in the cell count, and the greedy coloring stays within a fixed factor of it.
@testset "the 2D example's colors are bounded independently of the cell count" begin
    rows = Dict(:moving => Int[], :fixed => Int[])
    colors = Dict(:moving => Int[], :fixed => Int[])
    for (nx, ny) in ((6, 4), (10, 8), (16, 12), (24, 18))
        cp = setup(
            _example(; nx, ny);
            use_ale = true,
            slip = _slip,
            tensor = _anisotropic(π / 6),
        )
        push!(rows[:moving], _largest_row(cp.jacobian.colored))
        push!(colors[:moving], cp.jacobian.n_colors)
        m = cp.problem.mesh
        cells = [m.topology.cells[i].nodes for i = 1:n_cells(m.topology)]
        fixed = PyrolysisProblem(
            mesh = cell_mesh(
                m.positions,
                cells;
                tag = (i, k, s) -> ((:left, :right), (:bottom, :top))[k][s+1],
                fixed = true,
            ),
            material = cp.problem.material,
            bc_set = cp.problem.bc_set,
            T_initial = cp.problem.T_initial,
            ξ_initial = cp.problem.ξ_initial,
            tspan = cp.problem.tspan,
        )
        cpf = setup(fixed; tensor = _anisotropic(π / 6))
        push!(rows[:fixed], _largest_row(cpf.jacobian.colored))
        push!(colors[:fixed], cpf.jacobian.n_colors)
    end
    @info "2D example at 6×4, 10×8, 16×12 and 24×18 cells: largest rows and colors" rows colors
    for kind in (:moving, :fixed)
        @test all(==(rows[kind][end]), rows[kind][2:end])
        @test all(colors[kind] .<= 1.3 .* rows[kind])
    end
end
