"""
Type analysis of the two entry points, run in a fresh
process (`test/runtests.jl`), so no other test file's definitions enter the
analysis.

JET's `@test_call` analyses `residual!` and the colored `jacobian!` at the
concrete types `setup` builds, on a fixed mesh, a moving mesh with in-depth
radiation and lumped wall and Darcy transport, and a moving mesh under the
equidistribution correction, in both energy forms. Reports are limited to the
package's own modules.
"""

using Test
using Pyrolysis
using Pyrolysis.Materials: BEER_LAMBERT
using Pyrolysis.Assembly: residual!
using Pyrolysis.ColoredJacobian: jac_prototype, jacobian!
using Pyrolysis.Closures: EquidistributionCorrection
using Pyrolysis.Integration: setup

function _load_jet()
    try
        @eval using JET
    catch _
        @info "JET is not available; the type analysis is skipped"
        return false
    end
    # JET loads as an inert stub on Julia versions it does not support.
    if !(@eval JET.JET_AVAILABLE)
        @info "JET does not support Julia $(VERSION); the type analysis is skipped"
        return false
    end
    return true
end

function _pmma()
    virgin = SolidComponent(
        :virgin,
        φ = 0.05,
        λ = 1e-5,
        ρ = 1190.0,
        c = 1420.0,
        k = 0.21,
        ε = 0.86,
    )
    gas = GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02)
    rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
    return Material(name = :PMMA, components = (virgin, gas), reactions = (rxn,))
end

function _wood()
    virgin = SolidComponent(
        :wood,
        ρ = 500.0,
        c = (1500.0, 0.5),
        k = 0.15,
        ε = 0.9,
        α = 0.5,
        κ = 1e-14,
        λ = 1e-5,
        ρ_intrinsic = 1500.0,
    )
    gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03)
    rxn = Reaction(:pyrolysis, 1 => 2, A = 1e10, E = 150e3, h = -300e3)
    return Material(name = :wood, components = (virgin, gas), reactions = (rxn,))
end

_column(mat; top, stretch = nothing) = PyrolysisProblem(
    mesh = stretch === nothing ? create_uniform_mesh(0.01, 8) :
           create_stretched_mesh(0.01, 8, stretch),
    material = mat,
    bc_set = BoundaryConditionSet(thermal = (top = top, bottom = AdiabaticBC())),
    tspan = (0.0, 100.0),
)

const _configurations = (
    fixed = (problem = () -> _column(_pmma(); top = HeatFluxBC(50e3)), options = (;)),
    moving_beer_lambert_darcy = (
        problem = () -> _column(_wood(); top = RadiativeFluxBC(flux = 50e3)),
        options = (use_ale = true, radiation_model = BEER_LAMBERT),
    ),
    equidistribution = (
        problem = () -> _column(_pmma(); top = HeatFluxBC(50e3), stretch = 1.01),
        options = (
            use_ale = true,
            correction = EquidistributionCorrection(
                activation = (2.0, 0.5),
                monitor = (1.0, 1.0),
            ),
        ),
    ),
)

if _load_jet()
    include("type_analysis.jl")
else
    @testset "type analysis (skipped)" begin
        @test true
    end
end
