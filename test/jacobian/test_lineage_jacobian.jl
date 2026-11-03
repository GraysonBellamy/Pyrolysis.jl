# The colored Jacobian of a material with a reaction order other than one
# against central finite differences and the dense reference: at the initial
# state and near the conversion floor.

using Test
using SparseArrays
using Pyrolysis
using Pyrolysis.States
using Pyrolysis.Meshes
using Pyrolysis.Assembly: residual!
using Pyrolysis.ColoredJacobian
using Pyrolysis.Integration: setup, consistent_state!

# A uniform slab burning to gas at order `n` and a rate constant 0.02 s⁻¹.
function _uniform(cells; n = 2.0, inert = 0.0)
    L, T0, ρ = 10e-3, 600.0, 500.0
    a = SolidComponent(:a, ρ = ρ, φ = 0.3, c = 1000.0, k = 0.2, λ = 1e-6)
    b = SolidComponent(:b, ρ = ρ, φ = 0.3, c = 1000.0, k = 0.2, λ = 1e-6)
    gas = GaseousComponent(:gas, M = 0.03, c = 1000.0, k = 0.02)
    material = Material(
        name = :uniform,
        components = (a, b, gas),
        reactions = (Reaction(:burn, :a => :gas, A = 0.02, E = 0.0, h = 0.0, n = n),),
    )
    return PyrolysisProblem(
        mesh = uniform_mesh(L, cells),
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (
                top = AmbientGasBC(GasReservoir(P = 0.0, T = T0)),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> T0,
        ξ_initial = [z -> ρ * (1 - inert), z -> ρ * inert, z -> 0.0],
        tspan = (0.0, 300.0),
    )
end

# The Jacobian by central differences, the step scaled by each column's value
# and tolerance, by the cell size on a node position, whose own value is no
# scale (the back node sits at zero).
function _differences(cp, u, t, columns)
    def, ws = cp.def, cp.ws
    layout = def.layout
    L = length(u)
    J = zeros(L, L)
    plus, minus = zeros(L), zeros(L)
    v = copy(u)
    positions =
        is_moving(layout) ? [position_index(layout, nd) for nd = 1:n_nodes(layout)] : Int[]
    size = minimum(def.setup.reference.thickness)
    for c in columns
        δ = c in positions ? 1e-5 * size : 1e-6 * max(abs(u[c]), 1e2 * cp.abstol[c])
        v[c] = u[c] + δ
        residual!(plus, v, def, t, ws)
        v[c] = u[c] - δ
        residual!(minus, v, def, t, ws)
        v[c] = u[c]
        J[:, c] .= (plus .- minus) ./ 2δ
    end
    return J
end

function _checks(cp, u; t = 1.0, rtol = 1e-4, columns = eachindex(u))
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    Jd = dense_jacobian(def, ws, u, t)
    @test isempty(missing_entries(Jd, cj.pattern))
    J = jacobian!(jac_prototype(cj), cj, u, t)
    @test same_structure(J, cj.pattern)
    @test compare_entries(Matrix(J), Jd; rtol = 1e-12).passed
    # Against each row's own scale: the other columns keep the colored values.
    # At order zero a mass row's derivative on a position is exactly zero, which
    # differences of the row leave at 7e-11 of round-off.
    differences = Matrix(J)
    differences[:, columns] .= _differences(cp, u, t, columns)[:, columns]
    @test compare_entries(Matrix(J), differences; rtol, atol = 1e-9).passed
    return Matrix(J)
end

# The state of `cp` with the reactant at the conversion `x` of its reference in
# every cell, the rest of its mass in the gas, the mesh shrunk to the condensed
# volume that is left and the algebraic rows consistent.
function _at_conversion(cp, x)
    layout = cp.def.layout
    u = copy(cp.u0)
    for i = 1:n_cells(layout)
        m = u[mass_index(layout, i, 1)]
        u[mass_index(layout, i, 1)] = x * m
        u[mass_index(layout, i, 3)] += (1 - x) * m
    end
    return consistent_state!(u, cp.def, cp.ws, 1.0)
end

@testset "the colored Jacobian at order $n, $(energy) form" for n in (0.0, 0.5, 2.0),
    energy in (:temperature, :enthalpy)

    @testset "at the start" begin
        cp = setup(_uniform(8; n); energy_state = energy, use_ale = true)
        _checks(cp, copy(cp.u0))
    end
    @testset "on a fixed mesh" begin
        cp = setup(_uniform(8; n); energy_state = energy)
        _checks(cp, copy(cp.u0))
    end
    @testset "near the conversion floor" begin
        cp = setup(_uniform(8; n); energy_state = energy, use_ale = true)
        for x in (5e-4, 1e-6)
            _checks(cp, _at_conversion(cp, x))
        end
    end
end

@testset "the reference adds nothing to the pattern" begin
    second = setup(_uniform(6; n = 2.0); use_ale = true)
    first_order = setup(_uniform(6; n = 1.0); use_ale = true)
    @test second.jacobian.pattern == first_order.jacobian.pattern
end
