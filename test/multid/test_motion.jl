using Test
using LinearAlgebra
using StaticArrays
using OrdinaryDiffEq: ODEFunction, ODEProblem
import OrdinaryDiffEq
using Pyrolysis
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.States: state_length
using Pyrolysis.Assembly: residual!
using Pyrolysis.ColoredJacobian
using Pyrolysis.Integration: setup, rodas5p_policy, default_linear_solver

# The integrator the package runs: Rodas5P with the colored Jacobian and the
# time derivative by finite differences.
const _rodas5p = rodas5p_policy(default_linear_solver()).alg

const _L = 0.01

# One inert solid at rest on a smoothly distorted square of `n × n` cells.
function _uniform(n; ξ = 600.0)
    s = SolidComponent(:s, ρ = 800.0, c = 1000.0, k = 0.3)
    inside(x) = SVector(
        x[1] + 0.0006sin(π * x[1] / _L) * sin(2π * x[2] / _L),
        x[2] + 0.0006sin(π * x[2] / _L) * sin(2π * x[1] / _L),
    )
    mesh = box_mesh(
        range(0.0, _L; length = n + 1),
        range(0.0, _L; length = n + 1);
        transform = inside,
    )
    adiabatic = AdiabaticBC()
    return PyrolysisProblem(
        mesh = mesh,
        material = Material(name = :uniform, components = (s,), reactions = ()),
        bc_set = BoundaryConditionSet(
            thermal = (
                top = adiabatic,
                bottom = adiabatic,
                left = adiabatic,
                right = adiabatic,
            ),
        ),
        T_initial = x -> 500.0,
        ξ_initial = [x -> ξ],
        tspan = (0.0, 10.0),
    )
end

# A fixed-step Rodas5P integration of a set-up problem from `u0`.
function _integrate(cp, dt; u0 = cp.u0, tspan = (0.0, 10.0))
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    f = ODEFunction{true}(
        (du, u, p, t) -> residual!(du, u, def, t, ws);
        mass_matrix = mass_matrix(def.layout),
        jac = (J, u, p, t) -> (jacobian!(J, cj, u, t); nothing),
        jac_prototype = jac_prototype(cj),
    )
    sol =
        OrdinaryDiffEq.solve(ODEProblem(f, copy(u0), tspan), _rodas5p; dt, adaptive = false)
    @test sol.retcode == OrdinaryDiffEq.ReturnCode.Success
    return sol
end

# The derived state of `cp` at `u`.
_derived(cp, u) = (derive!(cp.ws.derived, u, cp.def.layout, cp.def.setup); cp.ws.derived)

# A uniform stretch of the mesh over the solid at rest: the 2D volumes are
# quadratic in the positions, so the free stream holds to the integrator's
# order, as on a cylinder.
_stretch(t) = 1 - 0.15 * (1 - cos(2π * t / 20))
const _compress = (x, t) -> (-0.15 * (2π / 20) * sin(2π * t / 20) / _stretch(t)) * x

@testset "free stream under a uniform stretch: the integrator's order" begin
    errors = Float64[]
    for dt in (0.8, 0.4, 0.2)
        cp = setup(
            _uniform(6);
            use_ale = true,
            w_override = _compress,
            abstol = 1e-12,
            reltol = 1e-12,
        )
        sol = _integrate(cp, dt)
        layout = cp.def.layout
        top = position_index(layout, n_nodes(layout), 2)
        @test sol.u[end][top] ≈ cp.u0[top] * _stretch(10.0) rtol = 1e-6
        ξ_err, T_err = 0.0, 0.0
        for u in sol.u
            ds = _derived(cp, u)
            ξ_err = max(ξ_err, maximum(abs, ds.ξ[1, :] ./ 600 .- 1))
            T_err = max(T_err, maximum(abs, ds.T .- 500))
        end
        @test T_err <= 1e-10
        push!(errors, ξ_err)
    end
    orders = log2.(errors[1:2] ./ errors[2:3])
    @info "2D free-stream error under a uniform stretch" errors orders
    @test errors[3] > 1e-11
    @test all(>(3.5), orders)
end

# An interior motion that leaves the outline in place, over a history pair
# that varies across the square: the inventories cross the faces, and the
# column sums, so the mass-weighted mean, stay put.
const _swirl =
    (x, t) -> 2e-4 * sin(π * x[1] / _L) * sin(π * x[2] / _L) * SVector(cos(0.3t), sin(0.3t))

@testset "the history pair's mass-weighted mean under a prescribed motion" begin
    cp = setup(
        _uniform(6);
        use_ale = true,
        w_override = _swirl,
        history = true,
        abstol = 1e-12,
        reltol = 1e-10,
    )
    layout = cp.def.layout
    n = n_cells(layout)
    u0 = copy(cp.u0)
    g = cp.ws.derived.geometry
    for i = 1:n
        M = u0[history_index(layout, i, 2)]
        @test M ≈ u0[mass_index(layout, i, 1)]
        χ = 0.5 + 0.4 * sin(2π * g.cell_center[i][1] / _L) * g.cell_center[i][2] / _L
        u0[history_index(layout, i, 1)] = χ * M
    end
    sol = _integrate(cp, 0.5; u0)
    G(u) = sum(u[history_index(layout, i, 1)] for i = 1:n)
    M(u) = sum(u[history_index(layout, i, 2)] for i = 1:n)
    m(u) = sum(u[mass_index(layout, i, 1)] for i = 1:n)
    moved = maximum(
        abs(sol.u[end][history_index(layout, i, 1)] - u0[history_index(layout, i, 1)])
        for i = 1:n
    )
    @test moved > 1e-6 * G(u0) / n
    for u in sol.u
        @test G(u) ≈ G(u0) rtol = 1e-13
        @test M(u) ≈ M(u0) rtol = 1e-13
        @test m(u) ≈ m(u0) rtol = 1e-13
        @test G(u) / M(u) ≈ G(u0) / M(u0) rtol = 1e-13
    end
end

# A single skewed cell whose node velocities follow a linear field `v = L x`:
# the cell's mean velocity gradient is `L` exactly, and its deformation
# gradient's row is `L F`.
@testset "the deformation gradient's row carries L F: $(D)D" for D in (2, 3)
    s = SolidComponent(:s, ρ = 800.0, c = 1000.0, k = 0.3)
    shear(x) =
        D == 2 ? SVector(x[1] + 0.4x[2], 0.9x[2] + 0.1x[1]) :
        SVector(x[1] + 0.4x[3], x[2] - 0.2x[3] + 0.1x[1], 0.9x[3])
    mesh =
        D == 2 ? box_mesh([0.0, 1.0], [0.0, 1.0]; transform = shear) :
        box_mesh([0.0, 1.0], [0.0, 1.0], [0.0, 1.0]; transform = shear)
    tags = boundary_tags(mesh.topology)
    problem = PyrolysisProblem(
        mesh = mesh,
        material = Material(name = :inert, components = (s,), reactions = ()),
        bc_set = BoundaryConditionSet(
            thermal = NamedTuple{Tuple(tags)}(Tuple(AdiabaticBC() for _ in tags)),
        ),
        T_initial = x -> 400.0,
        ξ_initial = [x -> 800.0],
        tspan = (0.0, 1.0),
    )
    cp = setup(problem; use_ale = true, closure = :balance, deformation = true)
    layout = cp.def.layout
    u = copy(cp.u0)
    @test [u[deformation_index(layout, 1, a, b)] for a = 1:D, b = 1:D] == I(D)
    L =
        D == 2 ? (@SMatrix [0.3 -0.7; 0.2 0.5]) :
        (@SMatrix [0.3 -0.7 0.1; 0.2 0.5 -0.4; 0.6 0.0 -0.2])
    F =
        D == 2 ? (@SMatrix [1.1 0.2; -0.1 0.9]) :
        (@SMatrix [1.1 0.2 0.0; -0.1 0.9 0.3; 0.05 0.0 1.2])
    for (nd, x) in enumerate(mesh.positions), k = 1:D
        u[velocity_index(layout, nd, k)] = (L*x)[k]
    end
    for a = 1:D, b = 1:D
        u[deformation_index(layout, 1, a, b)] = F[a, b]
    end
    du = zeros(state_length(layout))
    residual!(du, u, cp.def, 0.0, cp.ws)
    Ḟ = [du[deformation_index(layout, 1, a, b)] for a = 1:D, b = 1:D]
    @test Ḟ ≈ L * F atol = 1e-14
end
