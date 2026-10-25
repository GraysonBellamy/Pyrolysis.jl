# 9. The Solution Object and Plotting

Every call to `solve` returns a single `PyrolysisSolution`. This chapter is a complete
reference for that object: the fields it always carries, the `extras` series, the
deactivation events and the solver statistics, the accessor functions for pulling out
profiles and reconstructing the fields at arbitrary times, the conservation recorder and
audit, and the `Plots.jl` recipes shipped with the package. It closes with practical
post-processing recipes: surface temperature and mass-loss rate, moving meshes, and
curves for comparison against cone-calorimeter or TGA data. The numerics that produce
these quantities are derived in the Technical Reference (energy assembly in Technical
Reference §7, the surface energy balance in §12, deactivation in §14, the conservation
ledgers and the audit in §16).

A note on coordinates that governs every array below: `z = 0` is the **bottom /
substrate** (the face tagged `:bottom`) and `z = L` is the **top / exposed surface** (the
face tagged `:top`). Heat enters from the top and the material recedes downward.
"Surface" always means the top, exposed face; "back face" means the bottom.

---

## 9.1 What `solve` returns

`solve(problem; ...) -> PyrolysisSolution`. The return type is an immutable struct defined
in `src/Outputs/solution.jl`. You never construct it yourself: after the integration,
`solve` reconstructs the fields at every saved time from the saved states (one residual
evaluation per saved state, under the deactivation mask of that time) and hands the
result back.

The examples in this chapter use a 2 mm PMMA layer on 20 cells under 50 kW/m², on a
moving mesh with deactivation and a conservation recorder:

```julia
using Pyrolysis

pmma = Material(
    name = :PMMA,
    components = (
        SolidComponent(:virgin, ρ = 1190.0, φ=.3, c = 1420.0, k = 0.19, ε = 0.95, λ = 1e-6),
        GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02,),
    ),
    reactions = (Reaction(:decomposition, 1 => 2, A = 8.5e12, E = 188e3, h = 870e3),),
)
problem = PyrolysisProblem(
    mesh = create_uniform_mesh(2e-3, 20),
    material = pmma,
    bc_set = BoundaryConditionSet(
        thermal = (
            top = RadiativeFluxBC(flux = 50e3) + RadiativeBC(ε = 0.95, T_inf = 300.0) +
                  ConvectiveBC(h = 10.0, T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
        mass = (top = AmbientGasBC(GasReservoir(P=0., T=300.)), bottom = ImpermeableBC()),
    ),
    T_initial = z -> 300.0,
    ξ_initial = [z -> 1190.0, z -> 0.0],
    tspan = (0.0, 300.0),
    output_times = collect(0.0:2.0:300.0),
)

recorder = ConservationRecorder()
solution = solve(problem; use_ale = true, deactivation = true,
                 radiation_model = SURFACE_ABSORPTION, callback = recorder)
```

Its summary, printed at the REPL:

```julia
solution
# PyrolysisSolution(:Terminated, outcome :min_thickness, 55 saved times, 12.774 s)

display(solution)
# PyrolysisSolution
#   Return code: Terminated
#   Outcome: min_thickness
#   Saved times: 55
#   Time range: [0.0, 107.92160864026599] s
#   Live cells: 1 of 20
#   Invocations: 1
#   Accepted steps: 325
#   Solve time: 12.774 s
```

The layer burned through in about 108 s: one deactivation invocation removed 19 of the 20
cells together, and the thickness stop ended the solve when the last cell thinned below
the default `min_thickness`. The solve time depends on the machine.

`PyrolysisProblem`, `PyrolysisSolution`, `solve`, every accessor of §9.4 and the
conservation recorder and audit of §9.8 are exported by `using Pyrolysis`.

---

## 9.2 Core fields

Every field is present after every solve. The time-indexed arrays are filled at every
saved time (`solution.t`), so an index `i` into the time vector selects a consistent
snapshot across all of them. `n_t` is the number of saved times, `NC` the number of
components.

| Field | Type | Shape | Meaning / units |
|---|---|---|---|
| `problem` | `PyrolysisProblem` | — | The problem that was solved; its mesh is the mesh you built, unchanged |
| `t` | `Vector{Float64}` | `(n_t,)` | Saved times [s] |
| `T` | `Matrix{Float64}` | `(n_cells, n_t)` | Cell temperatures [K]; `NaN` in a deactivated cell |
| `ξ` | `Array{Float64,3}` | `(NC, n_cells, n_t)` | Concentrations `ξ_j` [kg·m⁻³]; `NaN` in a deactivated cell |
| `m` | `Array{Float64,3}` | `(NC, n_cells, n_t)` | Component masses per cell [kg]; zero in a deactivated cell |
| `z` | `Matrix{Float64}` | `(n_nodes, n_t)`; `0×0` on a fixed mesh | Node positions [m]; `NaN` beyond the outermost live node |
| `correction_velocity` | `Matrix{Float64}` | `(n_nodes, n_t)`; `0×0` on a fixed mesh | Node corrections `c = w − v_s`, the mesh velocity relative to the solid's [m/s]; zero without a correction; `NaN` beyond the outermost live node |
| `cell_centers` | `Matrix{Float64}` | `(n_cells, n_t)` | Cell-centre positions [m]; `NaN` in a deactivated cell |
| `boundary_T` | `Dict{Symbol,Vector{Float64}}` | per tag, `(n_t,)` | Area-weighted surface temperature of each boundary tag's faces [K]; after cells beyond a tag's faces are deactivated, of the face bounding the live cells on that side |
| `surface_T` | `Vector{Float64}` | `(n_t,)` | Surface temperature of the `:top` faces, `T_s` [K] |
| `mass_loss_rate` | `Vector{Float64}` | `(n_t,)` | Mass leaving through every boundary face, per unit current area of the `:top` faces [kg·m⁻²·s⁻¹] |
| `boundary_mass_flux` | `Dict{Symbol,Vector{Float64}}` | per tag, `(n_t,)` | Outward mass flux of every component summed, area-weighted over the tag's faces [kg·m⁻²·s⁻¹] |
| `boundary_blowing_ratio` | `Dict{Symbol,Vector{Float64}}` | per tag, `(n_t,)` | Area-weighted ratio `h/h₀` the tag's convective condition applied (1 without outflow or without a blowing model) |
| `extras` | `SolutionExtras` | — | Secondary series (§9.3) |
| `events` | `Vector{DeactivationEvent}` | — | The deactivation events, in order (§9.3) |
| `outcome` | `Symbol` | — | `:running` unless the deactivation or the thickness stop ended the solve: `:burned_out`, `:interior` or `:min_thickness` |
| `retcode` | `Symbol` | — | The integrator's return code |
| `solve_time` | `Float64` | — | Wall time of the solve [s] |
| `statistics` | `SolverStatistics` | — | The solver's counters (§9.3) |
| `raw_solution` | ODE solution | — | The integrator's solution: the raw states and the interpolant |
| `definition` | `CoreDefinition` | — | The definition of the solve, through which `solution(t)` reconstructs the fields |

A few points worth internalizing:

- **`T`, `ξ` and `m` are cell-centred.** Row `j` is cell `j`, counting from the bottom
  (`z = 0`) toward the top (`z = L`), so `solution.T[1, i]` is the cell at the substrate.
  The cell at the exposed surface is the outermost live one,
  `solution.T[solution.extras.live_cells[i], i]`; while every cell is live that is
  `solution.T[end, i]`. The *surface* temperature on the exposed face, as opposed to the
  adjacent cell centre, is the separately solved `solution.surface_T[i]` (§9.4).

- **`ξ` and `m` are component-major in their first index.** `solution.ξ[k, j, i]` is the
  concentration of component `k` in cell `j` at saved time `i`, in the order you listed the
  components in the `Material` constructor. These are *mass concentrations*
  `ξ_j = Y_j ρ` in kg·m⁻³, not mass fractions (Technical Reference §3, §4). The masses
  `m = ξ V` are the solver's states, per cell, for the mesh's cross section (1 m² by
  default for a planar mesh, so then they are also per unit area).

- **`z` is empty on a fixed mesh.** Without `use_ale = true` no node moves and
  `solution.z` is a `0×0` matrix; `solution.cell_centers` (or `get_cell_positions`,
  §9.4) gives the cell centres on either mesh.

- **Deactivated cells read `NaN`.** With `deactivation = true`, a deactivated cell reads
  `NaN` in `T`, `ξ` and `cell_centers` from the time it goes, and the nodes beyond the
  outermost live node read `NaN` in `z` and `correction_velocity`;
  `extras.live_cells` gives the live count (§9.7, Chapter 10).

- **`retcode` and `outcome` together say how the solve ended** (§8.14); `succeeded` and
  `completed` (§9.4) combine them.

### Direct field access and accessors

The fields are part of the solution's fixed schema, so reading them directly is always
safe:

```julia
times      = solution.t                 # Vector{Float64}
T_field    = solution.T                 # (n_cells, n_t)
T_at_t5    = solution.T[:, 5]           # profile at saved time index 5
T_surface  = solution.surface_T         # Vector{Float64}
mlr        = solution.mass_loss_rate    # Vector{Float64}
xi_virgin  = solution.ξ[1, :, end]      # component 1 profile at the last saved time
```

The accessor functions of §9.4 return the same arrays, or slices of them, and add the
reconstruction between saved times.

---

## 9.3 The `extras` series, the events and the statistics

`solution.extras` is a `SolutionExtras`, a typed NamedTuple whose keys are always all
present. Each series holds one value per saved time, indexed like the core fields, unless
noted. Quantities of the `:top` (surface) or `:bottom` (back) faces are `NaN` without that
tag, or once no live cell borders it (after a complete burnout, for example). Per-area
quantities divide by the current area of the faces they name unless noted.

| `extras` key | Units | Meaning |
|---|---|---|
| `gas_generation_rate` | kg·m⁻²·s⁻¹ | Gas the reactions of the live cells produce, per unit initial area of the `:top` faces |
| `back_temperature` | K | Temperature of the `:bottom` faces (`== boundary_T[:bottom]`) |
| `total_mass` | kg·m⁻² | The summed component masses, per unit initial area of the `:top` faces |
| `thickness_history` | m | Distance from the first node to the outermost live node |
| `min_cell_thickness` | m | Thickness of the thinnest live cell |
| `surface_mass_flux` | kg·m⁻²·s⁻¹ | Outward mass flux of the `:top` faces only (cf. `mass_loss_rate`) |
| `surface_heat_flux_net` | W·m⁻² | Heat flux the `:top` condition applies at the face, positive into the domain (under a prescribed temperature, the conducted flux); under `BEER_LAMBERT` the incident flux enters in depth and is not part of it |
| `surface_heat_flux_absorbed` | W·m⁻² | Incident radiation absorbed at the face |
| `surface_heat_flux_reradiation` | W·m⁻² | Re-radiative loss, positive leaving |
| `surface_heat_flux_convection` | W·m⁻² | Convective loss, positive leaving |
| `surface_heat_flux_conduction` | W·m⁻² | Flux conducted from the `:top` faces into their cells, `k (T_s − T_cell)/d` |
| `surface_blowing_ratio`, `surface_pooling_factor`, `surface_face_fraction`, `surface_face_buoyancy` | – | Blowing diagnostics of the `:top` surface balance (1, 1, 0 and 1 without outflow or without a blowing model; §5.2.3) |
| `volumetric_absorption` | W·m⁻² | Intensity the cells absorbed in depth, per unit exposed area; zero unless `radiation_model = BEER_LAMBERT` |
| `radiation_transmitted` | W·m⁻² | Intensity leaving through the last face of the in-depth beam, per unit exposed area; zero unless `BEER_LAMBERT` |
| `surface_emissivity` | – | Emissivity of the cell under the `:top` faces |
| `surface_conductivity` | W·m⁻¹·K⁻¹ | Conductivity of the cell under the `:top` faces |
| `back_heat_flux` | W·m⁻² | Heat flux the `:bottom` condition applies, positive into the domain |
| `MLR_gauge` | kg·m⁻²·s⁻¹ | Mass-loss rate per unit initial area of the `:top` faces (cone gauge convention) |
| `cross_section_area_history` | m² | The cross-section measure, `A(t)` (varies only under a lateral-shrinkage law) |
| `live_cells` | – | Number of live cells |
| `deactivation_times` | s | Time of every deactivation, one entry per deactivated cell, in order; not indexed by saved time |
| `live_mass` | kg | Mass of the live cells |
| `live_energy` | J | Energy of the live cells (§9.8) |
| `mass_outflow` | kg·s⁻¹ | Mass outflow of the live region through its boundary faces |
| `energy_rate` | W | Enthalpy rate of the live region (§9.8) |

For a top condition of incident radiation, re-radiation and convection, the surface
decomposition closes as `net = absorbed − reradiation − convection`, and the surface
balance makes `net` equal `surface_heat_flux_conduction` at convergence. Under
`BEER_LAMBERT`, `volumetric_absorption + radiation_transmitted` is the absorbed share
`a·q` of the incident flux.

```julia
gen_rate  = solution.extras.gas_generation_rate    # kg/(m²·s)
T_back    = solution.extras.back_temperature       # K
mass_m2   = solution.extras.total_mass             # kg/m²
thickness = solution.extras.thickness_history      # m
q_net     = solution.extras.surface_heat_flux_net  # W/m²
```

`get_gas_generation_rate(solution)` returns `solution.extras.gas_generation_rate`. Mind
its meaning: it is the *volumetric* production of the live cells, **not** the mass-loss
rate through the surface. The mass leaving the sample is `solution.mass_loss_rate` (per
current area) and `extras.MLR_gauge` (per initial area); §9.5 says which to use.

### The deactivation events

`solution.events` lists one `DeactivationEvent` per deactivated cell, in order (the type
is reachable through `Pyrolysis.Internal`). Its fields:

| Field | Meaning |
|---|---|
| `t` | Time of the event [s] |
| `invocation` | The deactivation callback's invocation it belongs to, numbered from one over the solve and shared by every cell one invocation deactivates |
| `cell` | The deactivated cell |
| `receiver` | The live cell that took its inventories (zero when none was left, the inventories then lost) |
| `fraction` | The cell's condensed mass as a fraction of its initial condensed mass when it went |
| `rows` | The state rows masked |
| `created_volume` | Volume the merge created to keep the pore space [m³] |
| `lost_mass` | Masses lost per component [kg] (the vented gas on a fixed mesh) |
| `lost_enthalpy` | Sensible enthalpy lost [J] |
| `regularization` | Caloric regularization the cells held that the receiver does not [J] |

A merge of several cells into one receiver records one event per dead cell, the outermost
first; the merge's created volume, losses and regularization are written to that first
event, so the sums over the events are the merges'. In the example of §9.1 the 19 events
share invocation 1:

```julia
length(solution.events)                    # 19
solution.statistics.invocations            # 1
unique(e.invocation for e in solution.events)
```

### The solver statistics

`solution.statistics` is a `SolverStatistics` with the fields `accepted`, `rejected`
(steps), `rhs_evaluations`, `jacobians`, `linear_solves`, `nonlinear_failures`,
`domain_rejections` (proposals the domain test rejected), `largest_constraint_residual`
(the largest scaled algebraic-row residual of the accepted states, one or below passing;
zero for a fixed-step solve) and `invocations` (deactivation invocations that deactivated
cells). §8.3 explains the domain test behind two of them.

---

## 9.4 Accessor functions and reconstruction

All accessors are exported by `using Pyrolysis`:

| Function | Returns |
|---|---|
| `get_times(sol)` | `sol.t` |
| `get_temperatures(sol)` | `sol.T` |
| `get_concentrations(sol)` | `sol.ξ` |
| `get_gas_generation_rate(sol)` | `sol.extras.gas_generation_rate` |
| `n_timesteps(sol)` | `length(sol.t)` |
| `get_cell_positions(sol; at = lastindex(sol.t))` | The cell centres at saved time `at` |
| `get_profile(sol, field, i)` | One field at saved time `i` |
| `sol(t)`, `interpolate_solution(sol, t)` | The fields reconstructed at any time `t` |
| `succeeded(sol)`, `completed(sol)` | How the solve ended |
| `conservation_audit(sol[, recorder])` | The mass and energy closures (§9.8) |

### `get_profile(solution, field, i)`

Pulls a single field at one saved time index. `field` is `:T` (the cell temperatures),
`:ξ` (the concentrations, `NC × n_cells`) or `:ξj` / `:xij` for component `j`, `:m` (the
component masses, `NC × n_cells`) or `:mj`, or `:z` (the node positions; on a fixed mesh
the mesh's). Deactivated cells are `NaN`, the masses included, and so are the nodes beyond
the outermost live one. Any other symbol raises an `ArgumentError`.

```julia
T_prof   = get_profile(solution, :T, 10)     # == solution.T[:, 10]
xi_all   = get_profile(solution, :ξ, 10)     # == solution.ξ[:, :, 10]   (NC × n_cells)
xi_gas   = get_profile(solution, :ξ2, 10)    # == solution.ξ[2, :, 10]
m_virgin = get_profile(solution, :m1, 10)    # component-1 masses [kg], NaN where dead
z_nodes  = get_profile(solution, :z, 10)     # node positions [m]
```

`i` is an index into `solution.t`, not a time in seconds. To find the index nearest a
target time, use `argmin(abs.(solution.t .- t_target))`.

### `get_cell_positions(solution; at)`

Returns the cell-centre positions at the saved time `solution.t[at]` [m], `NaN` for a
deactivated cell (and so for every cell beyond the outermost live node). The default is
the last saved time; the initial positions are `get_cell_positions(solution; at = 1)`
when the first saved time is the start of the time span. This is the depth axis for
plotting a profile on either mesh, fixed or moving.

```julia
z_end   = get_cell_positions(solution)            # at the last saved time
z_start = get_cell_positions(solution; at = 1)    # the initial cell centres
```

### `solution(t)` and `interpolate_solution(solution, t)`

Both return the fields at an arbitrary time `t` within the solution's span, as one saved
time of the solution holds them: a NamedTuple `(; t, T, ξ, m, z, cell_centers)`. They are
reconstructed from the integrator's interpolant between the saved states, through the
same state mapping the solver uses, with the deactivation mask of the last event at or
before `t` applied; `z` is empty on a fixed mesh. At a saved time they return exactly the
saved fields.

```julia
f = solution(42.7)                     # == interpolate_solution(solution, 42.7)
f.T                                    # cell temperatures at t = 42.7 s [K]
f.ξ[1, :]                              # component-1 concentrations [kg/m³]
f.z                                    # node positions [m] (moving mesh)
```

The interpolant between two saved states that a deactivation event separates is not the
trajectory, since the event changed the state inside that interval. A `t` strictly inside
such an interval is refused with an error naming the event:

```
ArgumentError: t = 0.125 lies between the saved states at 0.1 and 0.15 s, which the
deactivation of cell 12 at t = 0.14432365151720558 s (invocation 1) separates: the
interpolant there is not the trajectory
```

A `t` outside the span is refused as well. When you need values near an event, save at
the times you need (`output_times` or `saveat`, §8.12); every other interval interpolates
freely.

### Surface and back-face temperatures

The temperature of each boundary is solved at every saved time from that boundary's
surface energy balance: `solution.surface_T` for the top, and `solution.boundary_T` per
tag for every boundary (`solution.boundary_T[:bottom] == solution.extras.back_temperature`).
After surface cells are deactivated, `surface_T` is the temperature of the face bounding
the outermost live cell.

```julia
T_s_history    = solution.surface_T              # top face [K]
T_back_history = solution.boundary_T[:bottom]    # bottom face [K]
```

### `succeeded` and `completed`

`succeeded(solution)` is `true` when the solve finished: return code `:Success`, or
`:Terminated` by a complete burnout (outcome `:burned_out`) or the `min_thickness` stop
(outcome `:min_thickness`). An interior void (`:interior`), a `:Failure`, and a
termination by a callback of your own are not. `completed(solution)` is the strict test:
return code `:Success`, the end of the time span reached. The solve of §9.1 is
`succeeded` but not `completed`.

---

## 9.5 Mass-loss rate (MLR) and gas generation

Three related but distinct quantities are available. Choosing correctly matters when the
cross section changes (a lateral-shrinkage law) or for cone-gauge comparison.

1. **`solution.mass_loss_rate`**: the mass leaving through every boundary face, per unit
   *current* area `A(t)` of the `:top` faces [kg·m⁻²·s⁻¹]. The back face contributes only
   when its mass condition passes gas (zero for the default impermeable holder); the
   strictly top-face flux is `solution.extras.surface_mass_flux`. With a constant cross
   section this is the MLR.

2. **`solution.extras.MLR_gauge`**: the same mass per unit *initial* area `A₀`, the
   cone-calorimeter "gauge" convention: `MLR_gauge = mass_loss_rate · A(t)/A₀`. Under a
   lateral-shrinkage law this differs from the local MLR; with a constant area the two
   coincide. **For cone-calorimeter comparison, use `MLR_gauge`** (the lab measures mass
   per fixed sample-holder area).

3. **`solution.extras.gas_generation_rate`**: the gas the reactions of the live cells
   produce, per unit initial area [kg·m⁻²·s⁻¹]. This is the total gas being *created*, not
   the gas crossing the surface. It differs transiently from the surface flux because gas
   takes time to diffuse or flow out and can accumulate in the pore space. Use it to see
   where in the reaction the production peaks; use the surface quantities for what leaves
   the sample.

After a complete burnout no live cell borders the top, and the three read `NaN` from then
on; filter with `isfinite` before reducing.

```julia
times     = solution.t
mlr_local = solution.mass_loss_rate          # per A(t)
mlr_gauge = solution.extras.MLR_gauge        # per A₀ (cone convention)
gen_rate  = solution.extras.gas_generation_rate

# Peak MLR and the time it occurs:
i_peak   = argmax(x -> isfinite(mlr_gauge[x]) ? mlr_gauge[x] : -Inf, eachindex(mlr_gauge))
t_peak   = times[i_peak]
peak_mlr = mlr_gauge[i_peak]
```

### An HRR proxy

The package does not compute a heat-release rate: that requires a gas-phase combustion
model, which the condensed-phase solver omits (Technical Reference §1). The standard
*proxy* is the mass-loss rate scaled by an assumed effective heat of combustion `ΔH_c`:

```julia
ΔH_c      = 25e6                                 # J/kg, effective heat of combustion (assumed)
hrr_proxy = solution.extras.MLR_gauge .* ΔH_c    # W/m²
```

State this assumption plainly in any comparison; the proxy is only as good as `ΔH_c`.

---

## 9.6 Plotting with the `Plots.jl` recipes

The `PyrolysisPlotsExt` extension attaches a `Plots.jl` recipe for `PyrolysisSolution`. It
loads automatically when both `Pyrolysis` and `Plots` (or `RecipesBase`) are present in
your environment; no explicit `using PyrolysisPlotsExt` is needed. The recipe reads a
`kind` keyword to route to one of three plot types; the default is the
surface-temperature trace.

```julia
using Pyrolysis, Plots

plot(solution)                      # default: surface temperature vs time
plot(solution; kind = :surface_T)   # surface temperature vs time (explicit)
plot(solution; kind = :mlr)         # mass-loss rate vs time
plot(solution; kind = :heatmap)     # temperature field T(cell, t) as a heatmap
```

| `kind` | x | y | source field |
|---|---|---|---|
| `:surface_T` (default) | `t` [s] | surface T [K] | `sol.t`, `sol.surface_T` |
| `:mlr` | `t` [s] | MLR [kg·m⁻²·s⁻¹] | `sol.t`, `sol.mass_loss_rate` |
| `:heatmap` | `t` [s] | cell index (depth) | `sol.t`, `axes(sol.T, 1)`, `sol.T` |

Any other `kind` raises an `ArgumentError` listing the valid options.

Notes that follow from the recipe:

- The recipes read only `t`, `T`, `surface_T` and `mass_loss_rate`.
- `:mlr` plots `mass_loss_rate` (per current area). For the cone-gauge curve, plot
  `extras.MLR_gauge` yourself (§9.5).
- `:heatmap` puts the **cell index** on the y axis, not the physical depth; deactivated
  cells (`NaN`) are left blank. For a depth axis in metres, plot profiles against
  `get_cell_positions` (below).

Because these are standard recipes, the usual `Plots.jl` attributes pass through:

```julia
plot(solution; kind = :mlr, label = "50 kW/m²", lw = 2, color = :firebrick)
plot!(solution2; kind = :mlr, label = "25 kW/m²")   # overlay a second run
savefig("mlr_comparison.png")
```

### Plotting depth profiles by hand

There is no built-in recipe for "temperature versus depth at a fixed time", but it takes
a few lines. `get_cell_positions(solution; at = i)` gives the cell centres at saved time
`i` on a fixed or a moving mesh, with `NaN` for deactivated cells, which `Plots` skips:

```julia
using Pyrolysis, Plots

plt = plot(xguide = "z [mm]", yguide = "T [K]", title = "Temperature profiles")
for t_target in (10.0, 30.0, 60.0)
    i = argmin(abs.(solution.t .- t_target))
    z = get_cell_positions(solution; at = i)     # cell centres at that time [m]
    plot!(plt, z .* 1e3, solution.T[:, i]; label = "t = $(solution.t[i]) s")
end
plt
```

To plot a species field instead, swap `solution.T[:, i]` for `solution.ξ[k, :, i]`
(component `k`).

---

## 9.7 Moving meshes: post-processing ALE solutions

When you solve with `use_ale = true`, the node positions are part of the state and the
mesh recedes or expands in time (Technical Reference §11). This changes how you read
profiles:

- **`solution.z`** is `(n_nodes, n_t)`: the node positions in metres at each saved time,
  with `n_nodes = n_cells + 1`, `NaN` beyond the outermost live node.
- **`solution.cell_centers`** is `(n_cells, n_t)`, `NaN` in deactivated cells.
- **`solution.extras.thickness_history`** is the live thickness, from the first node to
  the outermost live node. On a fixed mesh it is the domain length until cells are
  deactivated, then the live span.
- **`solution.extras.min_cell_thickness`** is the thinnest live cell, and
  **`solution.correction_velocity`** the mesh's velocity relative to the solid at every
  node (zero without a correction; Chapter 10).
- With **deactivation** (`deactivation = true`) the live cell count drops over the run.
  The arrays keep the full cell count and the rows of deactivated cells read `NaN`.
  Filter `NaN` (`isfinite`) before reducing a moving-mesh profile.
- **`solution.extras.cross_section_area_history`** carries `A(t)` under a
  lateral-shrinkage law; otherwise it stays at `A₀`.

A depth profile at a fixed time uses the cell centres of that time:

```julia
using Pyrolysis, Plots

i = argmin(abs.(solution.t .- 30.0))         # nearest saved time to 30 s
z = get_cell_positions(solution; at = i)     # cell centres at this time [m]
T = solution.T[:, i]
live = isfinite.(T)                          # drop deactivated cells
plot(z[live] .* 1e3, T[live];
     xguide = "z [mm]", yguide = "T [K]", legend = false,
     title = "Temperature profile at t = $(round(solution.t[i], digits = 1)) s")
```

Surface recession over the whole run:

```julia
plot(solution.t, solution.extras.thickness_history .* 1e3;
     xguide = "time [s]", yguide = "thickness [mm]", legend = false,
     title = "Sample thinning")
```

To visualize the moving mesh itself, plot each node's trajectory (every column of
`solution.z` is the node set at one time):

```julia
plt = plot(xguide = "time [s]", yguide = "node position z [m]", legend = false)
for j in axes(solution.z, 1)
    plot!(plt, solution.t, solution.z[j, :])   # a node beyond the live mesh breaks off (NaN)
end
plt
```

On the deck of §9.1 the node trajectories show the equidistribution correction at work:
the nodes of the burning layer draw together evenly as it thins, and the layer leaves in
one deactivation invocation instead of cell by cell.

If the sample thins below `min_thickness`, or its last cell deactivates, the run ends
early with `solution.retcode == :Terminated` and `solution.outcome` set to
`:min_thickness` or `:burned_out`; `solution.t[end]` is then the burn-through time, not
`tspan[2]`, because the state a solve ends on is always saved.

---

## 9.8 The conservation recorder and audit

![Mass and sensible-enthalpy inventories with internal transfers, boundary outflow, event losses and energy sources including sensible-enthalpy formation.](../../assets/figures/14_conservation.svg)

**Figure 14. Conservation of the live inventory (schematic).** Shown for mass-balanced reactions without external mass sources. Internal face fluxes and conservative receiver transfers cancel over live cells. Boundary outflow and recorded event losses leave that inventory. Reactions conserve total mass; energy sources include reaction heat, in-depth radiation and sensible-enthalpy formation $\sum_{i,j}h_j(T_i)S_{\mathrm{mass},ij}$. Event energy losses include enthalpy and caloric regularization. The audit compares inventory change with integrated transport, sources and event losses.

`conservation_audit` closes the mass and energy ledgers of the live region: the change in
its inventories against the time integral of what entered it and what the deactivation
events took out. It reads its integrals from a `ConservationRecorder` that observed the
solve, or, without one, from the series the solution saves.

```julia
recorder = ConservationRecorder()
solution = solve(problem; use_ale = true, deactivation = true,
                 radiation_model = SURFACE_ABSORPTION, callback = recorder)

audit = conservation_audit(solution, recorder)
audit.method            # :recorder
audit.mass.relative     # relative mass defect
audit.energy.relative   # relative energy defect

rough = conservation_audit(solution)   # no recorder: the trapezoid over the saved states
rough.method            # :trapezoid
```

The result is `(; method, mass, energy)`, with `mass` and `energy` each a NamedTuple:

| Field | Meaning (mass in kg, energy in J) |
|---|---|
| `initial`, `final` | The live inventories at the start and at the end |
| `change` | `final − initial` |
| `inflow` | Time integral of the rate entering the live region: the negated mass outflow; the enthalpy rate |
| `lost` | What the deactivation events took out: their lost masses; their lost enthalpy and caloric regularization |
| `defect` | `change − inflow + lost`, zero for an exact ledger |
| `relative` | `defect` over the largest of `\|change\|`, `\|inflow\|` and `\|lost\|` |

The enthalpy rate is the energy entering through the boundary faces of the live cells,
by conduction or the boundary condition and with the enthalpy their mass fluxes carry,
plus the live cells' reaction heats. The live energy is the enthalpy inventory in the
enthalpy form, and in the temperature form the cells' energy at their temperatures, which
equals their sensible enthalpy `Σⱼ mⱼ hⱼ(T)` away from empty cells and is what a
deactivation event conserves. On a moving mesh a merge loses nothing; on a fixed mesh the
vented gas and its enthalpy are the `lost` terms, so the ledger still closes.

**With a recorder** (`method = :recorder`), the integrals come from three-point
Gauss–Legendre quadrature on every accepted step's interpolant, up to the time where an
event ends a step early, and the inventories are recorded at the start and after every
step; the defect is then the integrator's error and falls with the tolerance. On the deck
of §9.1 at the default tolerances it closes mass to 4e-7 and energy to 5e-5 relative
(measured on that run). The recorder's own fields are `mass_outflow` [kg],
`energy_inflow` [J], `initial` and `final` (each `(energy, mass)`) and `steps` (the
accepted steps it saw). A step a callback of your own ends early is one it cannot account
for; its integrals, and with them the audit, become `NaN`.

**Without a recorder** (`method = :trapezoid`), the integrals are the trapezoid rule over
the saved series `extras.mass_outflow` and `extras.energy_rate`, from the first to the last
saved state, and the inventories are `extras.live_mass` and `extras.live_energy`. Its
error is the saving interval's: on a 12-cell ablating test column saved every 0.05 s,
the trapezoid's relative mass defect was 2e-2 against the recorder's 2e-6 on the same
solve. Use it for a quick look; pass a recorder whenever the closure matters.

The series behind the trapezoid are available for your own reductions:

```julia
ts   = solution.t
M    = solution.extras.live_mass      # kg
E    = solution.extras.live_energy    # J
mdot = solution.extras.mass_outflow   # kg/s
Edot = solution.extras.energy_rate    # W
```

The ledger equations, the energy convention they book (the gas enthalpy carried across
the faces by the mass fluxes, DSC-style reaction heats, and the conditional
LTE-film conversion to the reservoir control surface) and the event records are derived in Technical Reference §16.

---

## 9.9 Post-processing for cone / TGA comparison

Putting the pieces together, here is a typical reduction of a cone-calorimeter run to the
curves you would overlay on experimental data. The setup (material, BCs, mesh) is shown
in full in the worked-examples chapter; here the focus is the *output* side.

```julia
using Pyrolysis, Plots

solution = solve(problem;
                 abstol = 1e-8, reltol = 1e-6,
                 use_ale = true, deactivation = true,
                 min_thickness = 0.01e-3)

@assert succeeded(solution)

times      = solution.t
T_surface  = solution.surface_T                 # K, exposed-face temperature
T_back     = solution.extras.back_temperature   # K, substrate-face temperature
mlr_gauge  = solution.extras.MLR_gauge          # kg/(m²·s), cone-gauge normalized
mass_m2    = solution.extras.total_mass         # kg/m², remaining mass per initial area

# Mass-loss fraction relative to the initial mass:
mass_frac  = mass_m2 ./ mass_m2[1]

# Surface and back temperature, then MLR, stacked:
p1 = plot(times, T_surface; label = "T_surface",
          xguide = "time [s]", yguide = "T [K]")
plot!(p1, times, T_back; label = "T_back")
p2 = plot(times, mlr_gauge; legend = false,
          xguide = "time [s]", yguide = "MLR [kg/(m²·s)]")
plot(p1, p2; layout = (2, 1), size = (700, 600))
```

For a **TGA-style** idealization (a thin, near-isothermal sample), the natural output is
the remaining mass fraction versus temperature. Because the package solves a spatially
resolved slab, take the column total mass and pair it with a representative temperature
(the surface, or a mean over the live cells):

```julia
mass_frac = solution.extras.total_mass ./ solution.extras.total_mass[1]
T_repr    = solution.surface_T     # or [sum(filter(isfinite, c)) / count(isfinite, c) for c in eachcol(solution.T)]
plot(T_repr, mass_frac; xguide = "T [K]", yguide = "m/m₀", legend = false,
     title = "Mass-loss curve")
```

For overlaying simulated and measured curves on a common time base, make the solver save
exactly where you want to compare: set `output_times` on the problem (or `saveat` on
`solve`) to the experiment's times, which avoids relying on the interpolant (and on
`solution(t)`'s refusal across an event's interval, §9.4):

```julia
problem = PyrolysisProblem(; mesh, material, bc_set,
                           tspan = (0.0, t_final),
                           output_times = exp_times)
solution = solve(problem)
# solution.extras.MLR_gauge is aligned 1:1 with exp_times (plus the end of a
# terminated solve, which is always saved).
```

---

## 9.10 Quick reference

```julia
# Core fields (always present; index i selects a saved time):
solution.t                      # times [s]
solution.T[:, i]                # cell temperatures at time i [K]   (bottom→top; NaN where dead)
solution.ξ[k, :, i]             # component k concentration profile [kg/m³]
solution.m[k, :, i]             # component k masses per cell [kg]
solution.z[:, i]                # node positions [m] (moving mesh; 0×0 on a fixed one)
solution.cell_centers[:, i]     # cell centres [m]
solution.surface_T[i]           # exposed-face surface temperature [K]
solution.boundary_T[:bottom][i] # back-face temperature [K]
solution.mass_loss_rate[i]      # MLR per current area (every boundary face) [kg/(m²·s)]
solution.events                 # deactivation events
solution.retcode                # :Success / :Terminated / :Failure / ...
solution.outcome                # :running / :burned_out / :interior / :min_thickness
solution.statistics             # accepted, rejected, jacobians, ..., invocations
solution.solve_time             # wall-clock [s]

# extras (always-present keys):
solution.extras.MLR_gauge                 # MLR per initial area (cone gauge)
solution.extras.back_temperature          # back-face T [K]
solution.extras.total_mass                # mass per A₀ [kg/m²]
solution.extras.thickness_history         # live thickness [m]
solution.extras.live_cells                # live cell count
solution.extras.deactivation_times        # one time per deactivated cell [s]
solution.extras.gas_generation_rate       # volumetric gas production [kg/(m²·s)]
solution.extras.surface_heat_flux_net     # net surface flux [W/m²]

# Accessors (exported by `using Pyrolysis`):
get_profile(solution, :T, i)              # field at one time index (:T, :ξ, :ξj, :m, :mj, :z)
get_cell_positions(solution; at = i)      # cell centres at time index i [m] (default: last)
solution(t)                               # fields at continuous t (== interpolate_solution(solution, t))
succeeded(solution), completed(solution)  # finished / reached tspan[2]
conservation_audit(solution, recorder)    # mass and energy closures

# Plotting recipes (need `using Plots`):
plot(solution)                  # :surface_T (default)
plot(solution; kind = :mlr)     # mass-loss rate
plot(solution; kind = :heatmap) # T(cell, t) by cell index
```

See the Technical Reference for the physics behind these outputs: energy assembly (§7),
the surface energy balance that produces `surface_T` (§12), mesh motion behind `z` and
`thickness_history` (§10, §11), deactivation and its events (§14), and the ledgers behind
the conservation audit (§16).
