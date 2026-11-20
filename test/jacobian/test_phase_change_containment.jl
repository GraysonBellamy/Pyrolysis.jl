"""The colored Jacobian of a material with a phase change: its pattern
contains the dense Jacobian's nonzeros, it equals the dense Jacobian per entry,
and the dense Jacobian equals central differences, on fixed and moving meshes
in both energy forms. (The containment suite of `test_containment.jl` is
recorded Jacobian by Jacobian by the deletion identity gate, so these cases
live apart from it.)"""

using Test
using ForwardDiff
using SparseArrays
using Pyrolysis
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.Assembly
using Pyrolysis.Assembly: residual!
using Pyrolysis.ColoredJacobian
using Pyrolysis.Integration: setup

# The state of `cp` with cell temperatures `T(i)` and concentrations `ξ(i)`
# [kg/m³]; every other unknown as at the start.
function _cell_state(cp, T, ξ)
    u = copy(cp.u0)
    layout = cp.def.layout
    V = cp.def.setup.reference.volume
    for i = 1:n_cells(layout)
        m = Tuple(ξ(i) .* V[i])
        for j in eachindex(m)
            u[mass_index(layout, i, j)] = m[j]
        end
        u[energy_index(layout, i)] =
            energy_state_value(energy_state(layout), cp.def.setup, T(i), m, i)
    end
    return u
end

# The pattern contains the dense Jacobian's nonzeros and the colored Jacobian
# equals it per entry, at a state where the residual is finite. Returns the
# dense Jacobian.
function _contained(cp, u; t = 1.0)
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

# Central differences of the residual over every column, the step following
# each unknown's magnitude: temperatures, enthalpies, masses and node
# positions differ by orders of magnitude.
function _central_difference_jacobian(cp, u; t = 1.0)
    J = zeros(length(u), length(u))
    du_p, du_m, du0, v = similar(u), similar(u), similar(u), copy(u)
    layout = cp.def.layout
    zero_gas = Dict(
        mass_index(layout, i, g) => i for i = 1:n_cells(layout) for g in
        Pyrolysis.Materials.gas_component_indices(cp.def.setup.material.components) if
        iszero(u[mass_index(layout, i, g)])
    )
    residual!(du0, u, cp.def, t, cp.ws)
    for j in eachindex(u)
        if haskey(zero_gas, j)
            # The C¹ positive part has the same derivative from either side
            # at zero. Stay on its negative branch: a central step crossing
            # the pressure-scale layer measures a secant, not that derivative.
            h = 1e-7 * cp.def.setup.reference.solid_mass[zero_gas[j]]
            v[j] = -h
            residual!(du_p, v, cp.def, t, cp.ws)
            v[j] = -2h
            residual!(du_m, v, cp.def, t, cp.ws)
            J[:, j] .= (3 .* (du0 .- du_p) .- (du_p .- du_m)) ./ (2h)
        else
            h = max(1e-7 * abs(u[j]), 1e-9)
            v[j] = u[j] + h
            residual!(du_p, v, cp.def, t, cp.ws)
            v[j] = u[j] - h
            residual!(du_m, v, cp.def, t, cp.ws)
            J[:, j] .= (du_p .- du_m) ./ (2h)
        end
        v[j] = u[j]
    end
    return J
end

# Wood holding water, its vapour and the pyrolysis gas: the phase change beside
# the pyrolysis, under the CAST spruce isotherm and under an activity function
# of temperature, whose heat takes a derivative nested inside the Jacobian's
# dual numbers.
function _wet(activity)
    wood = SolidComponent(
        :wood,
        ρ = 500.0,
        ρ_intrinsic = 1500.0,
        c = (1500.0, 0.5),
        k = 0.15,
        ε = 0.9,
        κ = 1e-14,
        λ = 1e-5,
    )
    water = LiquidComponent(:water, ρ = 1000.0, c = 4180.0, k = 0.6)
    vapour = GaseousComponent(:vapour, M = 0.018015, c = 1900.0, k = 0.025)
    gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03)
    return Material(
        name = :wet,
        components = (wood, water, vapour, gas),
        reactions = (Reaction(:pyrolysis, 1 => 4, A = 1e10, E = 150e3, h = -300e3),),
        phase_changes = (PhaseChange(:moisture, 2 => 3; activity, sorbents = (1,)),),
    )
end

const _wet_laws = (
    cast = CAST_SPRUCE,
    activity_function = ActivityFunction((X, T) -> 4X * exp(300 * (1 / 298.15 - 1 / T))),
)

_wet_column(activity, n) = PyrolysisProblem(
    mesh = uniform_mesh(0.01, n),
    material = _wet(activity),
    bc_set = BoundaryConditionSet(
        thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()),
    ),
    T_initial = z -> 300.0,
    ξ_initial = [z -> 450.0, z -> 45.0, z -> 0.01, z -> 0.0],
    tspan = (0.0, 100.0),
)

# Physical states on both sides of equilibrium and above the isotherm's cap,
# cells without water, water and vapour at trace levels of either sign, and
# the trial states of the perf configurations.
function _wet_states(cp)
    n = n_cells(cp.def.layout)
    u0 = cp.u0
    uniform(T, ξ) = _cell_state(cp, i -> T, i -> ξ)
    return (;
        u0,
        mixed = _cell_state(
            cp,
            i -> 330.0 + 80.0 * i / n,
            i -> (
                450.0 * (0.7 + 0.3 * i / n),
                5.0 + 40.0 * i / n,
                0.02 + 0.3 * i / n,
                0.01 * i / n,
            ),
        ),
        evaporating = uniform(420.0, (420.0, 30.0, 0.05, 0.2)),
        condensing = uniform(320.0, (450.0, 45.0, 1.0, 0.0)),
        above_the_cap = uniform(370.0, (400.0, 150.0, 0.2, 0.0)),
        dry_cells = _cell_state(
            cp,
            i -> 380.0,
            i -> (450.0, isodd(i) ? 0.0 : 30.0, 0.1, 0.0),
        ),
        trace = _cell_state(
            cp,
            i -> 360.0,
            i -> (450.0, (isodd(i) ? -1 : 1) * 1e-9, (iseven(i) ? -1 : 1) * 1e-9, 0.0),
        ),
        perturbed = u0 .* (1 .+ 0.2 .* sin.(1:length(u0))),
        zero_state = zeros(length(u0)),
        negative = -abs.(u0),
    )
end

@testset "a phase change, $(law): $(moving ? "moving" : "fixed") mesh, $(energy) form" for law in
                                                                                           keys(
        _wet_laws,
    ),
    moving in (false, true),
    energy in (:temperature, :enthalpy)

    cp = setup(_wet_column(_wet_laws[law], 8); energy_state = energy, use_ale = moving)
    layout = cp.def.layout
    states = _wet_states(cp)
    for (state, u) in pairs(states)
        @testset "$state" begin
            Jd = _contained(cp, u)
            if state === :mixed
                # the step couples each cell's water, vapour and energy
                i = 3
                w, v, e = mass_index(layout, i, 2),
                mass_index(layout, i, 3),
                energy_index(layout, i)
                @test Jd[v, w] != 0 && Jd[w, v] != 0
                @test Jd[v, e] != 0 && Jd[e, v] != 0 && Jd[e, w] != 0
            end
        end
    end
    # the dense Jacobian against central differences
    for state in (:mixed, :condensing, :above_the_cap)
        @testset "central differences, $state" begin
            u = states[state]
            Jd = dense_jacobian(cp.def, cp.ws, u, 1.0)
            @test compare_entries(Jd, _central_difference_jacobian(cp, u); rtol = 1e-6).passed
        end
    end
end
