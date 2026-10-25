using Test
using LinearAlgebra
using StaticArrays
using Pyrolysis
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.States: state_length
using Pyrolysis.Assembly: residual!
using Pyrolysis.Closures: BalanceRowClosure
using Pyrolysis.Integration: setup

function _wood()
    virgin = SolidComponent(
        :wood,
        ρ_intrinsic = 1500.0,
        λ = 0.0,
        ρ = 500.0,
        c = 1500.0,
        k = 0.15,
        ε = 0.9,
    )
    char = SolidComponent(
        :char,
        ρ_intrinsic = 1800.0,
        λ = 0.0,
        ρ = 150.0,
        c = 1100.0,
        k = 0.08,
        ε = 0.95,
    )
    gas = GaseousComponent(:gas, M = 0.06, c = 1800.0, k = 0.03)
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

_top() =
    RadiativeFluxBC(flux = 5e4) +
    RadiativeBC(T_inf = 300.0) +
    ConvectiveBC(h = 10.0, T_inf = 300.0)
_T(z) = 600.0 + 2e4 * z
_ξ = (z -> 450.0 - 1e4 * z, z -> 15.0 + 1e3 * z, z -> 0.5)

# A 1 cm column of `n` stretched cells on a line, and the same column as `n`
# rectangles one cell wide, its width the line's cross section, its lateral
# faces adiabatic (and, the gas not moving, impermeable).
function _columns(n; energy)
    line = PyrolysisProblem(
        mesh = create_stretched_mesh(0.01, n, 1.1),
        material = _wood(),
        bc_set = BoundaryConditionSet(thermal = (top = _top(), bottom = AdiabaticBC())),
        T_initial = _T,
        ξ_initial = collect(_ξ),
        tspan = (0.0, 1.0),
    )
    cp1 = setup(line; energy_state = energy, use_ale = true)
    z = cp1.def.setup.reference.positions
    width = cp1.def.setup.reference.cross_section
    column = PyrolysisProblem(
        mesh = box_mesh([0.0, width], z),
        material = _wood(),
        bc_set = BoundaryConditionSet(
            thermal = (
                top = _top(),
                bottom = AdiabaticBC(),
                left = AdiabaticBC(),
                right = AdiabaticBC(),
            ),
        ),
        T_initial = x -> _T(x[2]),
        ξ_initial = [x -> f(x[2]) for f in _ξ],
        tspan = (0.0, 1.0),
    )
    cp2 = setup(column; energy_state = energy, use_ale = true, closure = :balance)
    return cp1, cp2, width
end

# A state of the line as a state of the column: each level's two nodes at the
# line's node, moving with its velocity, laterally at rest.
function _lift(u1, L1, L2, width)
    u2 = zeros(state_length(L2))
    for i = 1:n_cells(L1)
        for j = 1:n_components(L1)
            u2[mass_index(L2, i, j)] = u1[mass_index(L1, i, j)]
        end
        u2[energy_index(L2, i)] = u1[energy_index(L1, i)]
    end
    for nd = 1:n_nodes(L1), (side, x) in ((1, 0.0), (2, width))
        node = 2(nd - 1) + side
        u2[position_index(L2, node, 1)] = x
        u2[position_index(L2, node, 2)] = u1[position_index(L1, nd)]
        u2[velocity_index(L2, node, 2)] = u1[velocity_index(L1, nd)]
    end
    return u2
end

# Physical and trial states of the line: its start, a mixed reacting profile,
# a hot state, a perturbation of every unknown (the nodes by a twentieth of a
# cell, so they stay ordered), negative masses and temperatures.
function _states(cp)
    L = cp.def.layout
    n = n_cells(L)
    s = cp.def.setup
    V = s.reference.volume
    cell(T, ξ) = begin
        u = copy(cp.u0)
        for i = 1:n
            m = Tuple(ξ(i) .* V[i])
            for j in eachindex(m)
                u[mass_index(L, i, j)] = m[j]
            end
            u[energy_index(L, i)] = energy_state_value(energy_state(L), s, T(i), m, i)
        end
        u
    end
    perturbed = cp.u0 .* (1 .+ 0.2 .* sin.(1:length(cp.u0)))
    h = minimum(s.reference.thickness)
    for nd = 1:n_nodes(L)
        r = position_index(L, nd)
        perturbed[r] = cp.u0[r] + 0.05h * sin(r)
    end
    return (;
        start = cp.u0,
        mixed = cell(
            i -> 450.0 + 300.0 * i / n,
            i -> (450.0 * (0.4 + 0.5i / n), 20.0 + 60.0i / n, 0.3),
        ),
        hot = cell(i -> 950.0, i -> (300.0, 40.0, 1.0)),
        perturbed,
        negative = cell(i -> -300.0 - i, i -> (-450.0, -15.0i, -0.5)),
    )
end

@testset "an N×1 column reproduces the line's rows: $(energy) form" for energy in (
    :temperature,
    :enthalpy,
)
    cp1, cp2, width = _columns(8; energy)
    L1, L2 = cp1.def.layout, cp2.def.layout
    @test cp2.def.closure isa BalanceRowClosure
    @test length(cp2.def.closure.ties) == n_cells(L1)
    # The set-up states agree: the line's start lifted is the column's start.
    @test _lift(cp1.u0, L1, L2, width) ≈ cp2.u0 rtol = 1e-13
    for (name, u1) in pairs(_states(cp1))
        @testset "$name" begin
            u2 = _lift(u1, L1, L2, width)
            d1 = zeros(state_length(L1))
            d2 = zeros(state_length(L2))
            residual!(d1, u1, cp1.def, 0.5, cp1.ws)
            # The scale of a velocity row: the largest volume change rate.
            θV = maximum(abs, cp1.ws.contributions.θ .* cp1.ws.derived.V_eval)
            residual!(d2, u2, cp2.def, 0.5, cp2.ws)
            @test all(isfinite, d1) && all(isfinite, d2)
            rows = Dict{Symbol,Vector{Tuple{Float64,Float64}}}()
            put!(kind, a, b) = push!(get!(rows, kind, Tuple{Float64,Float64}[]), (a, b))
            for i = 1:n_cells(L1)
                for j = 1:n_components(L1)
                    put!(:mass, d1[mass_index(L1, i, j)], d2[mass_index(L2, i, j)])
                end
                put!(:energy, d1[energy_index(L1, i)], d2[energy_index(L2, i)])
            end
            for nd = 1:n_nodes(L1)
                left, right = 2nd - 1, 2nd
                put!(:position, d1[position_index(L1, nd)], d2[position_index(L2, left, 2)])
                put!(
                    :position,
                    d1[position_index(L1, nd)],
                    d2[position_index(L2, right, 2)],
                )
                put!(:velocity, d1[velocity_index(L1, nd)], d2[velocity_index(L2, left, 2)])
                # The lateral rows and the ties vanish on a lifted state.
                @test d2[position_index(L2, left, 1)] == 0
                @test d2[position_index(L2, right, 1)] == 0
                @test d2[velocity_index(L2, left, 1)] == 0
                @test d2[velocity_index(L2, right, 1)] == 0
                nd > 1 && @test d2[velocity_index(L2, right, 2)] == 0
            end
            for (kind, p) in rows
                scale = max(maximum(abs(a) for (a, _) in p), kind === :velocity ? θV : 0.0)
                worst = maximum(abs(a - b) for (a, b) in p)
                @test worst <= 1e-12 * scale
            end
        end
    end
end
