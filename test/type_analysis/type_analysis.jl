"""
The type analysis itself; `runtests.jl` includes it once JET has loaded, since
`@test_call` cannot be expanded without it.
"""

using JET

@testset "type analysis: $(name), $(energy) form" for name in keys(_configurations),
    energy in (:temperature, :enthalpy)

    config = _configurations[name]
    cp = setup(config.problem(); energy_state = energy, config.options...)
    u = copy(cp.u0)
    du = similar(u)
    J = jac_prototype(cp.jacobian)
    residual!(du, u, cp.def, 1.0, cp.ws)
    jacobian!(J, cp.jacobian, u, 1.0)
    @test_call target_modules = (Pyrolysis,) residual!(du, u, cp.def, 1.0, cp.ws)
    @test_call target_modules = (Pyrolysis,) jacobian!(J, cp.jacobian, u, 1.0)
end
