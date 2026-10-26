"""
PyrolysisPlotsExt — Plots.jl recipes for the solution of `Pyrolysis`.

Loaded automatically when both `Pyrolysis` and `Plots` (or `RecipesBase`)
are available in the user's environment. Consumes the solution's `t`, `T`,
`surface_T` and `mass_loss_rate`.

# Usage

```julia
using Pyrolysis, Plots
sol = solve(problem)

plot(sol)                   # surface temperature trace (default)
plot(sol; kind = :heatmap)  # T(z, t) heatmap
plot(sol; kind = :mlr)      # mass-loss-rate curve
plot(sol; kind = :surface_T)# surface-temperature trace (explicit)
```
"""
module PyrolysisPlotsExt

using Pyrolysis: PyrolysisSolution
using RecipesBase: @recipe, @series

# Default `plot(sol)` is the surface-temperature trace — inexpensive and
# gives a quick at-a-glance view of the solve.
@recipe function f(sol::PyrolysisSolution; kind = :surface_T)
    if kind === :surface_T
        seriestype --> :line
        xguide --> "time [s]"
        yguide --> "surface temperature [K]"
        title --> "Surface temperature"
        legend --> false
        return sol.t, sol.surface_T
    elseif kind === :mlr
        seriestype --> :line
        xguide --> "time [s]"
        yguide --> "mass-loss rate [kg/(m²·s)]"
        title --> "Mass-loss rate"
        legend --> false
        return sol.t, sol.mass_loss_rate
    elseif kind === :heatmap
        seriestype --> :heatmap
        xguide --> "time [s]"
        yguide --> "cell index (depth)"
        title --> "Temperature T(z, t)"
        colorbar_title --> "T [K]"
        # `sol.T` is (n_cells, n_steps); pass time on x, cell index on y.
        return sol.t, axes(sol.T, 1), sol.T
    else
        throw(
            ArgumentError(
                "Unknown plot kind: $(repr(kind)). " *
                "Expected :surface_T, :mlr, or :heatmap.",
            ),
        )
    end
end

end # module
