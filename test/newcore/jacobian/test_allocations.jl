using Test
using Pyrolysis
using Pyrolysis.Materials: BEER_LAMBERT, SURFACE_ABSORPTION
using Pyrolysis.States
using Pyrolysis.Assembly
using Pyrolysis.Assembly: residual!
using Pyrolysis.ColoredJacobian
using Pyrolysis.Integration: setup

# The colored Jacobian evaluates the residual at one dual-number type in
# buffers built once for it, and allocates nothing per call: the budget is the
# measurement. A per-cell allocation (a boxed element of a heterogeneous
# tuple, an untyped buffer) costs tens to hundreds of KiB at these sizes.
const JACOBIAN_BUDGET = 0
const N_CELLS = 100

function _pmma()
    virgin = SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86)
    gas = GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)
    rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
    return Material(name = :PMMA, components = (virgin, gas), reactions = (rxn,))
end

function _wood(; α = 1000.0, κ = 1e-14, law = nothing)
    virgin = SolidComponent(
        :wood,
        ρ = 500.0,
        c = (1500.0, 0.5),
        k = 0.15,
        ε = 0.9,
        α = α,
        κ = κ,
        λ = 1e-5,
        ρ_intrinsic = 1500.0,
    )
    gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5)
    rxn = Reaction(:pyrolysis, 1 => 2, A = 1e10, E = 150e3, h = -300e3)
    return Material(
        name = :wood,
        components = (virgin, gas),
        reactions = (rxn,),
        lateral_shrinkage_law = law,
    )
end

_column(mat; top) = PyrolysisProblem(
    mesh = create_uniform_mesh(0.01, N_CELLS, Val(2)),
    material = mat,
    bc_set = BoundaryConditionSet(thermal = (top = top, bottom = AdiabaticBC())),
    tspan = (0.0, 100.0),
)

# The configurations of `test/perf/perf_helpers.jl`, a variable cross section,
# and inventory fluxes relative to a prescribed mesh velocity. The wood's
# heat capacity depends on the temperature, so the per-component tuples mix
# constants and functions of it under automatic differentiation.
const _incident = RadiativeFluxBC(flux = 50e3)
const _configurations = (
    baseline_eulerian = (
        problem = () -> _column(_pmma(); top = HeatFluxBC(50e3)),
        options = (;),
    ),
    eulerian_surface_absorption = (
        problem = () -> _column(_pmma(); top = _incident),
        options = (radiation_model = SURFACE_ABSORPTION,),
    ),
    baseline_eulerian_darcy = (
        problem = () -> _column(_wood(); top = HeatFluxBC(50e3)),
        options = (;),
    ),
    ale_no_radiation = (
        problem = () -> _column(_pmma(); top = HeatFluxBC(50e3)),
        options = (use_ale = true,),
    ),
    ale_beer_lambert_darcy = (
        problem = () -> _column(_wood(α = 0.5); top = _incident),
        options = (use_ale = true, radiation_model = BEER_LAMBERT),
    ),
    ale_cross_section = (
        problem = () -> _column(_wood(law = χ -> 1 - 0.2χ); top = HeatFluxBC(50e3)),
        options = (use_ale = true,),
    ),
    relative_fluxes = (
        problem = () -> _column(
            Material(name = :inert, components = _pmma().components, reactions = ());
            top = HeatFluxBC(5e3),
        ),
        options = (use_ale = true, w_override = (x, t) -> -1e-3 * x),
    ),
)

# Allocations of one residual and one colored Jacobian after warm-up, measured
# inside a function.
function _allocations(cp)
    u = copy(cp.u0)
    du = similar(u)
    J = jac_prototype(cp.jacobian)
    residual!(du, u, cp.def, 1.0, cp.ws)
    jacobian!(J, cp.jacobian, u, 1.0)
    jacobian!(J, cp.jacobian, u, 1.0)
    r = @allocated residual!(du, u, cp.def, 1.0, cp.ws)
    j = @allocated jacobian!(J, cp.jacobian, u, 1.0)
    return r, j
end

@testset "allocations: $(name), $(energy) form" for name in keys(_configurations),
    energy in (:temperature, :enthalpy)

    config = _configurations[name]
    cp = setup(config.problem(); energy_state = energy, config.options...)
    r, j = _allocations(cp)
    @test r == 0
    @test j <= JACOBIAN_BUDGET
end

@testset "one residual method" begin
    # Every physics option is a type parameter of the definition, not a
    # method of the residual.
    @test length(methods(residual!)) == 1
end
