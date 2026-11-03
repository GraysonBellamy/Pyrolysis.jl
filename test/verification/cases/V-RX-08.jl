# =============================================================================
# V-RX-08 — a graded conversion against the closed form
# =============================================================================
#
# A 10 mm slab at a uniform 600 K, adiabatic, of a reactant A at the fraction
# f(z₀) of the column's mass, falling to the top, and a filler of the same density
# stripped at 2 s⁻¹. With `basis = 1` the reference is the cell's whole initial
# mass, so A starts at the conversion x₀ = f(z₀), and every cell, keeping its
# material, follows
#
#   D_t x = −k xⁿ,   x(t) = [x₀^(1−n) − (1−n) k t]^(1/(1−n)).
#
# V-RX-06 holds one conversion over the column; here the column holds a range of
# them, the state a heated deck reaches. The window keeps x at 20 x_th or more,
# where the onset is exactly one, and the rate is the kernel's own.

function _vrx08_problem(n_cells)
    L0, T0, ρ, f_top = 10e-3, 600.0, 500.0, 0.2
    A = SolidComponent(:A, ρ = ρ, φ = 0.3, c = 1000.0, k = 0.2, λ = 1e-6)
    filler = SolidComponent(:filler, ρ = ρ, φ = 0.3, c = 1000.0, k = 0.2, λ = 1e-6)
    gas = GaseousComponent(:gas, M = 0.03, c = 1000.0, k = 0.02)
    material = Material(
        name = :VRX08,
        components = (A, filler, gas),
        reactions = (
            Reaction(:burn, 1 => 3; A = 0.02, E = 0.0, h = 0.0, n = 0.5, basis = 1.0),
            Reaction(:strip, 2 => 3; A = 2.0, E = 0.0, h = 0.0),
        ),
    )
    f(z) = 0.99 - (0.99 - f_top) * z / L0
    return PyrolysisProblem(
        mesh = uniform_mesh(L0, n_cells),
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (
                top = AmbientGasBC(GasReservoir(P = 0.0, T = T0)),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> T0,
        ξ_initial = [z -> ρ * f(z), z -> ρ * (1 - f(z)), z -> 0.0],
        tspan = (0.0, 30.0),
        output_times = collect(range(0.0, 30.0, length = 31)),
    )
end

# The conversion of a material point that starts at `x0`.
function _vrx08_x(x0, t)
    p = sqrt(x0) - 0.5 * 0.02 * t
    return p > 0 ? p^2 : 0.0
end

# The cells' references: the condensed mass each started with.
_vrx08_reference(sol) = [sum(view(sol.m, 1:2, i, 1)) for i in axes(sol.m, 2)]

# The largest relative errors over the saved times of the conversions and of
# the kernel's rate of A.
function _vrx08_errors(sol, material)
    R = _vrx08_reference(sol)
    x0 = view(sol.m, 1, :, 1) ./ R
    rxn = first(material.reactions)
    conversion, rate = 0.0, 0.0
    for j in eachindex(sol.t)
        x = _vrx08_x.(x0, sol.t[j])
        conversion = max(conversion, maximum(abs, view(sol.m, 1, :, j) ./ R ./ x .- 1))
        kernel = sum(
            Pyrolysis.Physics.reaction_rate(
                rxn,
                (sol.m[1, i, j], sol.m[2, i, j], sol.m[3, i, j]),
                600.0,
                R[i];
                limiter = material.depletion_limiter,
            ) for i in eachindex(R)
        )
        rate = max(rate, abs(kernel / (0.02 * sum(R .* sqrt.(x))) - 1))
    end
    return (; conversion, rate)
end

function _verification_case_vrx08()
    function qois(sol)
        m = _vrx08_errors(sol, _vrx08_problem(20).material)
        drift = maximum(abs(T - 600.0) for T in sol.T)
        return [
            QoIResult(
                name = "x(t) of every cell",
                value = m.conversion,
                exact = 0.0,
                error = m.conversion,
                tol = 7e-13,  # 3× observed 2.3e-13
                comment = "max relative error against the closed form",
            ),
            QoIResult(
                name = "the kernel's rate",
                value = m.rate,
                exact = 0.0,
                error = m.rate,
                tol = 3.2e-14,  # 3× observed 1.1e-14
                comment = "max relative error of the rate of A against the sum",
            ),
            QoIResult(
                name = "isothermal: T stays at 600 K",
                value = drift,
                exact = 0.0,
                error = drift,
                tol = 1.9e-11,  # 3× observed 6.0e-12
                comment = "K",
            ),
        ]
    end
    return VerificationCase(
        id = "V-RX-08",
        title = "a graded conversion on the material",
        tags = [:ale, :kinetics, :reaction_order, :exact],
        references = ["closed form x(t) at each material point"],
        build = _vrx08_problem,
        solve_kwargs = (use_ale = true, reltol = 1e-8, abstol = 1e-10),
        qois = qois,
        ci = CIConfig(n_cells = 20, budget_s = 30.0),
        report = (
            curves = sol -> begin
                R = _vrx08_reference(sol)
                [
                    ReportCurve(
                        name = :conversion_history,
                        title = "V-RX-08: conversion of the cells",
                        xlabel = "t (s)",
                        ylabel = "x",
                        series = [
                            CurveSeries(
                                label = "cell $i",
                                x = sol.t,
                                y = sol.m[1, i, :] ./ R[i],
                                exact = t -> _vrx08_x(sol.m[1, i, 1] / R[i], t),
                            ) for i in (1, 10, 20)
                        ],
                    ),
                ]
            end,
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vrx08)
