# 9. The Solution Object and Plotting

Every call to `solve` returns a single `PyrolysisSolution`. This chapter is a complete
reference for that object: the concrete fields it always carries, the typed `extras`
NamedTuple of opt-in diagnostics, the accessor functions for pulling out profiles and
interpolating at arbitrary times, and the `Plots.jl` recipes shipped with the package.
It closes with practical post-processing recipes — computing surface temperature and
mass-loss rate, reading the conservation ledger, handling moving (ALE) meshes, and
preparing curves for comparison against cone-calorimeter or TGA data. The numerics that
produce these quantities are derived in the Technical Reference (energy assembly in
Technical Reference §7, the boundary-temperature solve in §12, conservation diagnostics
in §16).

A note on coordinates that governs every array below: `z = 0` is the **bottom /
substrate** (boundary id 2, tag `:bottom`) and `z = L` is the **top / exposed surface**
(boundary id 1, tag `:top`). Heat enters from the top and the material recedes downward.
"Surface" always means the top, exposed face; "back face" means the bottom.

---

## 9.1 What `solve` returns

`solve(problem; ...) -> PyrolysisSolution`. The return type is a single mutable struct.
You never construct it yourself — `solve` allocates it, populates the arrays during
solution extraction, and hands it back.

```julia
using Pyrolysis

solution = solve(problem; abstol = 1e-8, reltol = 1e-6)

solution                       # PyrolysisSolution(:Success, 601 timesteps, 12.4s)
typeof(solution)               # PyrolysisSolution{...}
```

The struct is defined in `src/Problem/pyrolysis_problem.jl`. Its layout is a **fixed core
schema** of concretely typed fields plus one typed `extras` NamedTuple. The fixed core is
what the plotting recipes and accessors rely on; `extras` holds secondary diagnostics that
are always present (the keys never change) but whose contents depend on the run.

`PyrolysisProblem`, `PyrolysisSolution`, and `solve` are all top-level exports of
`Pyrolysis`. Many of the convenience accessors discussed in §9.4 are *not* — they live in
submodules. The distinction matters for which `using` line you write, and is called out
explicitly where it applies.

---

## 9.2 Core fields

These fields are always present after a successful `solve`. They are populated for every
saved output time (`solution.t`), so an index `i` into the time vector selects a consistent
snapshot across all of them.

| Field | Type | Shape | Meaning / units |
|---|---|---|---|
| `problem` | `PyrolysisProblem` | scalar | The original problem (mesh, material, BCs) |
| `t` | `Vector{Float64}` | `(n_steps,)` | Saved solution times [s] |
| `T` | `Matrix{Float64}` | `(n_cells, n_steps)` | Cell-center temperatures [K] |
| `ξ` | `Array{Float64,3}` | `(NC, n_cells, n_steps)` | Mass concentrations `ξ_j` [kg·m⁻³] |
| `z` | `Matrix{Float64}` | `(n_nodes, n_steps)` | Node positions [m]; `0×0` in Eulerian mode |
| `surface_T` | `Vector{Float64}` | `(n_steps,)` | Top surface temperature `T_s` [K] |
| `mass_loss_rate` | `Vector{Float64}` | `(n_steps,)` | MLR per current area (top + permeable-bottom) [kg·m⁻²·s⁻¹] |
| `energy_residual` | `Vector{Float64}` | `(n_steps,)` | Energy-conservation residual; zeros if tracker off |
| `retcode` | `Symbol` | scalar | ODE-solver return code |
| `solve_time` | `Float64` | scalar | Wall-clock solve time [s] |
| `raw_solution` | `Any` | — | Raw `DifferentialEquations.jl` solution (interpolation source) |
| `extras` | `NamedTuple` | — | Opt-in diagnostics (see §9.3) |

A few points worth internalizing:

- **`T` and `ξ` are cell-centered.** Row `j` is cell `j`, counting from the bottom
  (`z = 0`) toward the top (`z = L`). So `solution.T[end, i]` is the cell nearest the
  exposed surface and `solution.T[1, i]` is the cell at the substrate. The actual *surface*
  temperature (on the exposed face, not the adjacent cell center) is the separately solved
  `solution.surface_T[i]` — see §9.5.

- **`ξ` is component-major in its first index.** `solution.ξ[k, j, i]` is the
  concentration of component `k` in cell `j` at saved time `i`. The component order matches
  the order you listed components in the `Material` constructor. These are *mass
  concentrations* `ξ_j = Y_j ρ` in kg·m⁻³, the primary species state variable — not mass
  fractions (Technical Reference §3, §4).

- **`z` is empty in Eulerian mode.** When you solve without `use_ale=true`, the mesh is
  fixed, so node positions are constant and `solution.z` is a `0×0` matrix. Use
  `get_cell_positions(solution)` (§9.4) to get the fixed cell centers instead. In ALE mode
  `z` is `(n_nodes, n_steps)`; see §9.7.

- **`retcode`** is `Symbol(raw_solution.retcode)`. Common values:
  `:Success` (completed to `t_end`), `:Terminated` (a callback stopped it — e.g. the ALE
  thickness-termination callback fired because the sample fully depleted),
  `:MaxIters`, `:DtLessThanMin`.

- **`energy_residual`** is all zeros unless a conservation tracker wrote into it. It is the
  series the `:residual` plot recipe consumes (§9.6). For the richer per-snapshot ledger,
  enable `diagnostics=true` and read `extras.diagnostics_log` (§9.8).

### Direct field access vs. accessors

The shipped examples read fields directly, which is the most robust path because the field
names are part of the fixed schema:

```julia
times      = solution.t                 # Vector{Float64}
T_field    = solution.T                 # (n_cells, n_steps)
T_at_t5    = solution.T[:, 5]           # profile at saved time index 5
T_surface  = solution.surface_T         # Vector{Float64}
mlr        = solution.mass_loss_rate    # Vector{Float64}
xi_virgin  = solution.ξ[1, :, end]      # component 1 profile at last time
```

There is also a small set of accessor functions (§9.4) — but be careful which package
namespace exports them.

### Display and quick checks

`PyrolysisSolution` has a pretty `show`, and a couple of boolean helpers:

```julia
solution                     # one-line summary
display(solution)            # multi-line summary (retcode, time range, peaks)

solution.retcode == :Success # success check that needs no import
```

Two helper functions, `succeeded(sol)` (true iff `retcode == :Success`) and
`n_timesteps(sol)` (`== length(sol.t)`), are exported by the `Problem` submodule, **not** by
top-level `Pyrolysis`. With plain `using Pyrolysis` you must qualify them
(`Pyrolysis.succeeded(solution)`) or pull them in via `using Pyrolysis.Internal`. Reading
`solution.retcode` directly is the simplest success check and needs no import.

---

## 9.3 The `extras` NamedTuple

`solution.extras` is a typed NamedTuple whose **keys are fixed at construction** — they are
always all present, even if a particular run leaves some of them at their default (zeros, or
an empty `diagnostics_log`). The numeric vectors are sized `(n_steps,)` and indexed by the
same time index as the core fields. This is where most cone-calorimeter post-processing
quantities live.

| `extras` key | Units | Meaning |
|---|---|---|
| `gas_generation_rate` | kg·m⁻²·s⁻¹ | Volumetric gas production integrated over the column, gauge-normalized |
| `back_temperature` | K | Back-face (bottom) surface temperature at `z=0` |
| `total_mass` | kg·m⁻² | Total domain mass per initial cross-section area |
| `thickness_history` | m | Sample thickness `z[end]-z[1]`; constant in Eulerian mode |
| `surface_mass_flux` | kg·m⁻²·s⁻¹ | Per-species summed **top-face** mass flux (excludes any bottom venting; cf. `mass_loss_rate`) |
| `surface_heat_flux_net` | W·m⁻² | Heat flux the boundary balance applies at the top face (`net = conduction + transpiration`); under `BEER_LAMBERT` it excludes the incident radiation, which is reported in `volumetric_absorption` |
| `surface_heat_flux_absorbed` | W·m⁻² | Absorbed incident (radiative) flux |
| `surface_heat_flux_reradiation` | W·m⁻² | Re-radiation loss `εσ(T_∞⁴−T_s⁴)` |
| `surface_heat_flux_convection` | W·m⁻² | Convective exchange `h_conv(T_∞−T_s)` |
| `surface_heat_flux_conduction` | W·m⁻² | Heat conducted into the slab at the surface |
| `surface_heat_flux_transpiration` | W·m⁻² | Enthalpy carried by escaping gas (zero with default S_conv) |
| `volumetric_absorption` | W·m⁻² | Beer–Lambert in-depth absorbed flux `α_surf·q_ext·(1−e^{−τ_tot})`; zero unless `radiation_model = BEER_LAMBERT` |
| `surface_emissivity` | – | Composition-weighted `ε_eff` at the surface |
| `surface_conductivity` | W·m⁻¹·K⁻¹ | Effective surface conductivity `k_eff` |
| `back_heat_flux` | W·m⁻² | Heat flux at the back face |
| `MLR_gauge` | kg·m⁻²·s⁻¹ | MLR normalized to initial area `A_0` (cone gauge convention) |
| `cross_section_area_history` | m² | Cross-section `A(t)` history (varies only in WithChi mode) |
| `diagnostics_log` | — | `Vector{DiagnosticsRow}`; empty unless `diagnostics=true` (§9.8) |

```julia
gen_rate  = solution.extras.gas_generation_rate    # kg/(m²·s)
T_back    = solution.extras.back_temperature       # K
mass_m2   = solution.extras.total_mass             # kg/m²
thickness = solution.extras.thickness_history      # m
q_net     = solution.extras.surface_heat_flux_net  # W/m²
```

There is one exported accessor for `extras`:

```julia
# Exported by the Problem submodule; qualify or `using Pyrolysis.Internal`.
Pyrolysis.get_gas_generation_rate(solution)   # == solution.extras.gas_generation_rate
```

Note its docstring caveat: `gas_generation_rate` is the *volumetric* production summed over
the domain, **not** the surface mass-loss rate. The canonical MLR is
`solution.mass_loss_rate` (per current area) and `extras.MLR_gauge` (per initial area). See
§9.5 for which one to use.

---

## 9.4 Accessor functions and interpolation

The package ships a handful of convenience functions. **None of them are top-level exports
of `Pyrolysis`.** `get_profile`, `get_cell_positions`, and `interpolate_solution` are
exported by the `Solver` submodule; `get_surface_temperature` and
`get_back_face_temperature` by the `Discretization` submodule. To call them unqualified,
bring them in via the umbrella namespace:

```julia
using Pyrolysis
using Pyrolysis.Internal     # re-exports every submodule's public symbol
```

or qualify each call (`Pyrolysis.get_profile(...)`). The umbrella namespace is documented as
the "escape hatch" for advanced users; symbols reached through it are not part of the stable
public API. For production scripts that only ever need the fixed fields, plain field access
(§9.2) avoids the issue entirely.

### `get_profile(solution, field, time_idx)`

Pulls a single field at one saved time index. `field` is `:T` for temperature, `:ξ` for the
full `(NC, n_cells)` concentration slab, or `:ξ1`, `:ξ2`, … (equivalently `:xi1`, `:xi2`,
…) for an individual component.

```julia
T_prof   = get_profile(solution, :T, 10)     # == solution.T[:, 10]
xi_all   = get_profile(solution, :ξ, 10)     # == solution.ξ[:, :, 10]  (NC × n_cells)
xi_char  = get_profile(solution, :ξ2, 10)    # == solution.ξ[2, :, 10]
```

`time_idx` is an index into `solution.t`, not a time in seconds. To find the index nearest a
target time, use `argmin(abs.(solution.t .- t_target))`.

### `get_cell_positions(solution)`

Returns the cell-center positions from the problem's mesh
(`solution.problem.mesh.cell_centers`). This is the right depth axis for plotting Eulerian
profiles, where the mesh is fixed and `solution.z` is empty.

```julia
z_centers = get_cell_positions(solution)     # Vector, length n_cells [m]
```

In ALE mode the mesh moves, so for moving-mesh plots you should reconstruct positions per
time step from `solution.z` instead (§9.7).

### `interpolate_solution(solution, t)` and the call syntax

Both interpolate the **raw state vector** at an arbitrary continuous time using the
underlying `DifferentialEquations.jl` dense interpolant. The state vector is packed
block-major `[T | ξ | z (ALE) | χ (WithChi)]` (Technical Reference §13, §15) — it is *not*
the convenient `(n_cells, n_steps)` view.

```julia
u_t = interpolate_solution(solution, 42.7)   # raw state vector at t = 42.7 s
u_t = solution(42.7)                          # identical: PyrolysisSolution is callable
```

Both error if `raw_solution === nothing` (e.g. a stub solution that was never solved). To
turn a raw state vector into physical quantities, use the `StateLayout` accessors. These
*are* top-level exports of `Pyrolysis` (and are also re-exported by the `Problem` submodule),
so plain `using Pyrolysis` is enough to call them unqualified:
`get_temperature`, `get_concentration`, and `get_node_positions`.

Mind their signatures — they are **per-cell / per-component scalar accessors**, not bulk
extractors:

- `get_temperature(u, layout, i)` returns the temperature of a *single* cell `i` (one scalar).
- `get_concentration(u, layout, i, j)` returns the concentration of component `j` in a
  *single* cell `i` (one scalar). Note the order: cell index first, component index second.
- `get_node_positions(u, layout)` is the exception — it returns a view of *all* node
  positions (ALE layouts only).

To rebuild a whole profile, loop over the cells. The `layout` comes from a freshly built
workspace:

```julia
def     = to_problem_def(solution.problem)
ws      = build_workspace(solution.problem, def)
layout  = ws.layout
n_cells = layout.n_cells

T_cells = [get_temperature(u_t, layout, i) for i in 1:n_cells]      # cell temperatures [K]
xi_1    = [get_concentration(u_t, layout, i, 1) for i in 1:n_cells] # component 1 [kg/m³]

# A single cell, if that is all you need:
T_cell1 = get_temperature(u_t, layout, 1)          # temperature of cell 1 [K]
```

For most workflows you do not need any of this — saving at the times you care about
(`output_times` / `saveat`) and indexing the core arrays is simpler and avoids rebuilding a
layout. Reach for interpolation, and these scalar accessors, only when you need a value
strictly between saved times.

### `get_surface_temperature(mesh)` / `get_back_face_temperature(mesh)`

These take a **mesh**, not a solution, and read the last-solved boundary state
(`mesh.boundary_states`). They reflect the mesh's state at the *end* of the solve (the
step-update callback syncs the mesh to the final integrator state), so they are a snapshot,
not a time history. For the time history of the surface temperature, use
`solution.surface_T` (top) and `solution.extras.back_temperature` (bottom), which are solved
at every saved time during extraction.

```julia
using Pyrolysis.Internal: get_surface_temperature, get_back_face_temperature
mesh = solution.problem.mesh
T_s_final    = get_surface_temperature(mesh)       # top face, final state [K]
T_back_final = get_back_face_temperature(mesh)     # bottom face, final state [K]
```

---

## 9.5 Mass-loss rate (MLR) and gas generation

Three related but distinct quantities are available. Choosing correctly matters when the
cross-section changes (variable-`A` / WithChi mode) or for cone-gauge comparison.

1. **`solution.mass_loss_rate`** — MLR per *current* area `A(t)`, [kg·m⁻²·s⁻¹]. It is the
   per-species-summed gas flux leaving the slab — the exposed face plus, when a permeable
   mass BC is configured there, the back face (zero for the default impermeable holder).
   The strictly top-face flux is `solution.extras.surface_mass_flux`. In constant-area runs
   (Eulerian, or ALE without a lateral-shrinkage law), this is the MLR.

2. **`solution.extras.MLR_gauge`** — surface MLR normalized to the *initial* area `A_0`, the
   cone-calorimeter "gauge" convention: `MLR_gauge = mass_loss_rate · A(t)/A_0`. With a
   variable cross-section this differs from the local MLR; with constant area the two
   coincide. **For cone-calorimeter comparison, use `MLR_gauge`** (the lab measures mass per
   fixed sample-holder area), unless your area is constant — in which case either works.

3. **`solution.extras.gas_generation_rate`** — *volumetric* gas production integrated over
   the whole column and divided by `A_0`, [kg·m⁻²·s⁻¹]. This is the total gas being
   *created* by reactions, not the gas crossing the surface. It can differ transiently from
   the surface flux because gas takes time to diffuse/advect out and can accumulate in the
   pore space. Use it to diagnose where in the reaction the production peaks; use the surface
   flux quantities for what actually leaves the sample.

```julia
times     = solution.t
mlr_local = solution.mass_loss_rate          # per A(t)
mlr_gauge = solution.extras.MLR_gauge        # per A_0 (cone convention)
gen_rate  = solution.extras.gas_generation_rate

# Peak MLR and the time it occurs:
i_peak    = argmax(mlr_gauge)
t_peak    = times[i_peak]
peak_mlr  = mlr_gauge[i_peak]
```

### An HRR proxy

The package does not compute a heat-release rate — that requires a gas-phase combustion
model the solid solver deliberately omits (Technical Reference §1). The standard *proxy* is
the mass-loss rate scaled by an assumed effective heat of combustion `ΔH_c`:

```julia
ΔH_c      = 25e6                                 # J/kg, effective heat of combustion (assumed)
hrr_proxy = solution.extras.MLR_gauge .* ΔH_c    # W/m²
```

State this assumption plainly in any comparison; the proxy is only as good as `ΔH_c`.

---

## 9.6 Plotting with the `Plots.jl` recipes

The `PyrolysisPlotsExt` extension attaches a `Plots.jl` recipe for `PyrolysisSolution`. It
loads automatically when both `Pyrolysis` and `Plots` (or `RecipesBase`) are present in your
environment — no explicit `using PyrolysisPlotsExt`. The recipe interprets a `kind` keyword
to route to one of four plot types; the default is the surface-temperature trace.

```julia
using Pyrolysis, Plots

solution = solve(problem)

plot(solution)                      # default: surface temperature vs time
plot(solution; kind = :surface_T)   # surface temperature vs time (explicit)
plot(solution; kind = :mlr)         # mass-loss rate vs time
plot(solution; kind = :residual)    # energy-conservation residual vs time
plot(solution; kind = :heatmap)     # temperature field T(z, t) as a heatmap
```

| `kind` | x | y | source field |
|---|---|---|---|
| `:surface_T` (default) | `t` [s] | surface T [K] | `sol.t`, `sol.surface_T` |
| `:mlr` | `t` [s] | MLR [kg·m⁻²·s⁻¹] | `sol.t`, `sol.mass_loss_rate` |
| `:residual` | `t` [s] | energy residual [W·m⁻²] | `sol.t`, `sol.energy_residual` |
| `:heatmap` | `t` [s] | cell index (depth) | `sol.t`, `axes(sol.T,1)`, `sol.T` |

Any unrecognized `kind` raises an `ArgumentError` listing the valid options.

Notes that follow from the recipe implementation:

- The recipes consume only the fixed core schema (`t`, `T`, `surface_T`, `mass_loss_rate`,
  `energy_residual`), so they work regardless of which `extras` you enabled.
- `:mlr` plots `mass_loss_rate` (per current area). For the cone-gauge curve, plot
  `extras.MLR_gauge` yourself (see §9.5).
- `:residual` plots `energy_residual`, which is zeros unless a conservation tracker filled
  it; a flat line at zero simply means the tracker was off.
- `:heatmap` uses **cell index** on the y axis, not physical depth. For a depth-resolved
  heatmap with metric units, build it yourself from `get_cell_positions` (Eulerian) or
  `solution.z` (ALE).

Because these are standard recipes, all the usual `Plots.jl` attributes pass through:

```julia
plot(solution; kind = :mlr, label = "25 kW/m²", lw = 2, color = :firebrick)
plot!(solution2; kind = :mlr, label = "50 kW/m²")   # overlay a second run
savefig("mlr_comparison.png")
```

### Plotting depth profiles by hand

There is no built-in recipe for "temperature vs depth at a fixed time," but it is a couple of
lines. In **Eulerian** mode use the fixed cell centers:

```julia
using Pyrolysis, Plots
using Pyrolysis.Internal: get_cell_positions

z = get_cell_positions(solution)            # fixed cell centers [m]
plt = plot(xguide = "depth z [mm]", yguide = "T [K]", title = "Temperature profiles")
for t_target in (10.0, 30.0, 60.0)
    i = argmin(abs.(solution.t .- t_target))
    plot!(plt, z .* 1e3, solution.T[:, i]; label = "t = $(solution.t[i]) s")
end
plt
```

To plot a species field instead, swap `solution.T[:, i]` for `solution.ξ[k, :, i]` (component
`k`).

---

## 9.7 Moving meshes: post-processing ALE solutions

When you solve with `use_ale=true`, node positions become part of the state and the mesh
recedes/expands in time (Technical Reference §11). This changes how you read profiles:

- **`solution.z`** is now `(n_nodes, n_steps)` — node positions in metres at each saved
  time. `n_nodes = n_cells + 1`.
- **`solution.extras.thickness_history`** is the sample thickness `z[end] − z[1]` at each
  time (it shrinks as the surface recedes). In Eulerian mode this is a constant equal to the
  domain length.
- If surface-cell **depletion** is enabled (`handle_depletion=true`, which requires
  `use_ale=true`), the cell count drops over the run. The output arrays are pre-allocated at
  the *maximum* cell count and rows for cells that no longer exist are filled with `NaN`.
  Always filter `NaN` before plotting or reducing a moving-mesh profile.
- **`solution.extras.cross_section_area_history`** carries `A(t)` if you used a
  lateral-shrinkage law (WithChi mode); otherwise it stays at `A_0`.

A depth profile at a fixed ALE time uses the cell midpoints reconstructed from the node
positions, and must skip the `NaN`-padded entries:

```julia
using Pyrolysis, Plots

i = argmin(abs.(solution.t .- 30.0))         # nearest saved time to 30 s
z_nodes = solution.z[:, i]                    # node positions at this time [m]
T_cells = solution.T[:, i]

valid = .!isnan.(T_cells)                     # drop depleted/padded cells
z_mid = 0.5 .* (z_nodes[1:end-1] .+ z_nodes[2:end])   # cell-center positions
plot(z_mid[valid] .* 1e3, T_cells[valid];
     xguide = "depth z [mm]", yguide = "T [K]",
     title = "ALE temperature profile at t = $(round(solution.t[i], digits=1)) s")
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
    plot!(plt, solution.t, solution.z[j, :])   # NaN gaps render as breaks after depletion
end
plt
```

If the sample fully depletes, the run terminates early (the ALE thickness-termination
callback fires when thickness drops below `min_thickness`), so `solution.retcode` will be
`:Terminated` and `solution.t[end]` will be the burn-through time, not `tspan[2]`.

---

## 9.8 Reading the conservation diagnostics

Pass `diagnostics=true` to `solve` to log snapshots of the conservation ledger. The solver
records a `DiagnosticsRow` every `diagnostics_stride` accepted steps (default 10) into
`solution.extras.diagnostics_log` (Technical Reference §16).

```julia
solution = solve(problem; diagnostics = true, diagnostics_stride = 10)

log = solution.extras.diagnostics_log    # Vector{DiagnosticsRow}; empty if diagnostics=false
```

Each row is a typed `NamedTuple` with these fields:

| Field | Units | Meaning |
|---|---|---|
| `t` | s | Simulation time of the snapshot |
| `E_matrix` | J | Matrix-only sensible energy `Σ ρc_p^eff T V` |
| `E_total` | J | Total sensible energy `Σ ρc_p^total T V` (matrix + gas) |
| `discrete_residual` | J | Discrete energy-balance residual |
| `discrete_relative` | – | `discrete_residual` relative to `|E_matrix(0)|` |
| `physical_residual` | J | Physical energy-balance residual |
| `physical_relative` | – | `physical_residual` relative to `|E_total(0)|` |
| `gas_storage_gap` | J | `ΔE_gas = E_total − E_matrix` (sensible heat in gas, not evolved) |
| `mass_max_relative_error` | – | Max relative mass-conservation error over components |

**Important caveat about snapshot mode.** With `diagnostics=true`, `solve` uses the
*snapshot* recorder, which evaluates only the instantaneous energies. In a snapshot row, the
cumulative-ledger fields — `discrete_residual`, `discrete_relative`, `physical_residual`,
`physical_relative`, `mass_max_relative_error` — are filled with `NaN`, because the
cumulative source integrals are not accumulated in this path. The fields you can trust from
`solve(...; diagnostics=true)` are `t`, `E_matrix`, `E_total`, and `gas_storage_gap`. The
full cumulative residuals require wiring `record_substep!` into the residual layer, which is
beyond the standard `solve` interface.

Reading the energies and the gas-storage gap (the quantity that justifies the matrix-only
heat-capacity convention, Technical Reference §16):

```julia
log      = solution.extras.diagnostics_log
ts       = [row.t for row in log]
E_matrix = [row.E_matrix for row in log]
E_total  = [row.E_total for row in log]
gas_gap  = [row.gas_storage_gap for row in log]

# The gas-storage gap should be a small fraction of the matrix energy
# (< ~0.2% at 1 atm; the physical-form residual should stay < 1%).
rel_gap = gas_gap ./ E_matrix
@show maximum(abs, rel_gap)
```

If you want a true running energy balance with non-`NaN` residuals, the lower-level
`Diagnostics` API (`EnergyConservationTracker`, `update_energy_tracker!`,
`check_energy_conservation`, `energy_conservation_report`) is available via
`using Pyrolysis.Internal`; that path is for advanced auditing and is described in the
Technical Reference §16.

A reminder tied to the energy convention (Technical Reference §7, §16): the default
`energy_form = :standard` carries gas sensible enthalpy through the advection source
`S_conv` and **excludes** a separate gas-generation sink. The surface transpiration BC is
mutually exclusive with that convention — `solve` rejects `use_transpiration_bc=true` to
avoid double-counting — so under default settings
`extras.surface_heat_flux_transpiration` is zero. If you switch to
`energy_form = :with_generation_sink`, the conservative-form sink is added back; the
diagnostics ledger reflects whichever form you chose.

---

## 9.9 Post-processing for cone / TGA comparison

Putting the pieces together, here is a typical reduction of a cone-calorimeter run to the
curves you would overlay on experimental data. The setup (material, BCs, mesh) is shown in
full in the worked-examples chapter; here we focus on the *output* side.

```julia
using Pyrolysis, Plots

solution = solve(problem;
                 abstol = 1e-8, reltol = 1e-6,
                 use_ale = true, handle_depletion = true,
                 min_thickness = 0.01e-3)

@assert solution.retcode in (:Success, :Terminated)

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
remaining mass fraction versus temperature. Because the package solves a spatially resolved
slab, take the column total mass and pair it with a representative temperature (the surface,
or a column average):

```julia
mass_frac = solution.extras.total_mass ./ solution.extras.total_mass[1]
T_repr    = solution.surface_T                  # or vec(sum(solution.T, dims=1)) ./ size(solution.T, 1)
plot(T_repr, mass_frac; xguide = "T [K]", yguide = "m/m₀", legend = false,
     title = "Mass-loss curve")
```

For overlaying simulated and measured curves on a common time base, interpolate onto the
experiment's sample times rather than the solver's saved times. The cleanest approach is to
make the solver save exactly where you want to compare — set `output_times` on the problem
(or `saveat` on `solve`) to the experiment's times — which avoids any reliance on the raw
interpolant:

```julia
problem = PyrolysisProblem(; mesh, material, bc_set,
                           tspan = (0.0, t_final),
                           output_times = exp_times)
solution = solve(problem)
# Now solution.extras.MLR_gauge is aligned 1:1 with exp_times.
```

---

## 9.10 Quick reference

```julia
# Core fields (always present; index i selects a saved time):
solution.t                      # times [s]
solution.T[:, i]                # cell temperatures at time i [K]   (bottom→top)
solution.ξ[k, :, i]             # component k concentration profile [kg/m³]
solution.z[:, i]                # node positions [m] (ALE; 0×0 in Eulerian)
solution.surface_T[i]           # exposed-face surface temperature [K]
solution.mass_loss_rate[i]      # MLR per current area (top + permeable-bottom) [kg/(m²·s)]
solution.retcode                # :Success / :Terminated / ...
solution.solve_time             # wall-clock [s]

# extras (always-present keys; some default to zeros):
solution.extras.MLR_gauge                 # MLR per initial area (cone gauge)
solution.extras.back_temperature          # back-face T [K]
solution.extras.total_mass                # domain mass per A_0 [kg/m²]
solution.extras.thickness_history         # sample thickness [m]
solution.extras.gas_generation_rate       # volumetric gas production [kg/(m²·s)]
solution.extras.surface_heat_flux_net     # net surface flux [W/m²]
solution.extras.diagnostics_log           # Vector{DiagnosticsRow} (needs diagnostics=true)

# Accessors (Solver/Discretization submodule exports — qualify or `using Pyrolysis.Internal`):
get_profile(solution, :T, i)              # field at one time index
get_cell_positions(solution)              # fixed cell centers (Eulerian) [m]
interpolate_solution(solution, t)         # raw state vector at continuous t  (== solution(t))
get_surface_temperature(mesh)             # final-state top face T [K]
get_back_face_temperature(mesh)           # final-state bottom face T [K]

# Plotting recipes (need `using Plots`):
plot(solution)                  # :surface_T (default)
plot(solution; kind = :mlr)     # mass-loss rate
plot(solution; kind = :residual)# energy residual
plot(solution; kind = :heatmap) # T(z,t) by cell index
```

See the Technical Reference for the physics behind these outputs: energy assembly and the
`S_conv` convention (§7), the boundary-temperature Newton solve that produces `surface_T`
(§12), ALE node motion behind `z` and `thickness_history` (§10, §11), and the dual energy
ledgers behind the diagnostics rows (§16).
