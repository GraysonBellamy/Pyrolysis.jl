# 6. Meshes and Geometry

The mesh is the spatial backbone of every Pyrolysis.jl simulation: it discretizes the
through-thickness coordinate `z` into the finite-volume cells over which temperature and
species concentrations are solved. Pyrolysis.jl is a strictly **one-dimensional** (1D)
through-thickness solver, so "geometry" here means a stack of cells along `z`, plus a single
scalar cross-sectional area `A` that turns those 1D cells into physical volumes. This chapter
shows how to create meshes (uniform and stretched), how to choose resolution, what the
cross-sectional area `A_0` means in a 1D model, how to inspect a mesh, and how the spatial
coordinate is laid out in the raw state vector — including the special case where the mesh
itself becomes part of the solved state (ALE). For the underlying finite-volume theory see
Technical Reference §13; for the moving-mesh (ALE) formulation see Technical Reference §11.

> **Coordinate convention (read this first).** Throughout Pyrolysis.jl, `z = 0` is the
> **bottom / substrate** boundary (tag `:bottom`, internal `boundary_id = 2`) and `z = L` is
> the **top / exposed surface** (tag `:top`, internal `boundary_id = 1`). Heat enters from the
> top; the material recedes toward the bottom. Cell indices increase from bottom (cell 1, near
> `z = 0`) to top (cell `n`, near `z = L`). Positive flux is transport in the `+z` direction.
> (You may see a stray comment in `examples/06_cone_calorimeter.jl` saying "z=0 is the exposed
> TOP surface" — that comment is wrong; the mesh generator and the rest of the package use the
> convention stated here.)

---

## 6.1 The mesh in one dimension

A 1D mesh in Pyrolysis.jl consists of:

- **Nodes** — the vertices that bound cells. For `n` cells there are `n + 1` nodes. Node 1
  sits at `z = 0` (bottom), node `n+1` sits at `z = L` (top).
- **Cells** — the finite-volume control volumes. Cell `i` spans nodes `i` and `i+1`, i.e. the
  interval `[z_i, z_{i+1}]`. Temperature `T_i` and the per-component mass concentrations
  `ξ_{j,i}` are stored at the cell **center** `z_i^c = (z_i + z_{i+1}) / 2`.
- **Faces** — in 1D each face is a single point that co-locates with a node, so there are also
  `n + 1` faces. Fluxes (conduction, gas transport) are evaluated at faces. Faces 1 and `n+1`
  are the boundary faces.

The mesh container type is `Mesh1D{NC}`, parameterized by the number of chemical components
`NC` (the number of `SolidComponent`/`LiquidComponent`/`GaseousComponent` entries in your
`Material`). You never construct `Mesh1D` directly; you call one of the two mesh-generation
functions described in §6.2. The relationship `n_faces = n_nodes = n_cells + 1` always holds.

Each cell carries a `CellState{NC}` (temperature, an `NTuple{NC}` of concentrations, and a
cached volume); `solve` populates them from the problem's initial conditions on its own copy
of the mesh (the mesh you built is left unchanged). The
volume of cell `i` is computed geometrically from the node positions and the cross-sectional
area:

```
V_i = |z_{i+1} − z_i| · A
```

See Technical Reference §13 for the full description of the cell/node/face layout and the
divergence operator built on it.

---

## 6.2 Creating a mesh

Two mesh creators are exported at the top level (`using Pyrolysis` makes both available):

| Function | Purpose |
|---|---|
| `create_uniform_mesh(L, n_cells, Val(NC); ...)` | Equal cell spacing `Δz = L / n_cells`. |
| `create_stretched_mesh(L, n_cells, Val(NC), stretch; ...)` | Geometrically graded cells (fine near one surface, coarse near the other). |

Both return a `Mesh1D{NC}`.

### 6.2.1 Uniform meshes

```julia
create_uniform_mesh(L::Real, n_cells::Int, ::Val{NC};
                    T_initial::Real = 300.0,
                    ξ_initial = nothing,
                    cross_section_area::Real = 1.0) -> Mesh1D{NC}
```

Positional arguments:

- `L` — domain thickness in metres (`m`). This is the total material length `z = L`.
- `n_cells` — number of finite-volume cells.
- `Val(NC)` — the number of chemical components, passed as a **type parameter** via
  `Val(...)`. `NC` must equal the number of components in the `Material` you intend to use.
  Passing a plain integer (e.g. `3` instead of `Val(3)`) will not compile.

Keyword arguments and their defaults:

| Keyword | Default | Meaning |
|---|---|---|
| `T_initial` | `300.0` | Uniform initial cell temperature `[K]`. Note: `solve` re-applies the problem's `T_initial` closure, so this matters mainly when you inspect the bare mesh. |
| `ξ_initial` | `nothing` | Initial mass concentrations `[kg·m⁻³]`. `nothing` → all zero; a scalar → that value for every component; an `NTuple`/vector of length `NC` → per-component values. |
| `cross_section_area` | `1.0` | Cross-sectional area `A_0` `[m²]` (see §6.4). |

A 1 cm slab with 100 cells and 2 components:

```julia
using Pyrolysis

mesh = create_uniform_mesh(0.01, 100, Val(2))   # L = 1 cm, 100 cells, NC = 2
```

A typical cone-calorimeter sample (5.8 mm, 72 cells, 3 components — virgin solid, char,
pyrolysis gas — initialized with the virgin solid at its bulk density):

```julia
sample_thickness = 5.8e-3      # 5.8 mm  [m]
n_cells          = 72          # 0.08 mm resolution

mesh = create_uniform_mesh(
    sample_thickness,
    n_cells,
    Val(3);                              # NC = 3: virgin, char, gas
    T_initial = 300.0,
    ξ_initial = (1210.0, 0.0, 0.0),      # full virgin density, no char, no gas
)
```

The node positions are `z_i = (i − 1)·Δz` with `Δz = L / n_cells`, so cell `i` occupies
`[(i−1)Δz, iΔz]` and every cell has the same volume `V_i = Δz · A`.

### 6.2.2 Stretched (graded) meshes

In a cone-calorimeter or fire-exposure problem the steepest temperature and reaction-rate
gradients are at the **top** (exposed) surface, where the pyrolysis front lives. A graded mesh
puts small cells where the action is and large cells deep in the cold substrate, giving the
same accuracy for fewer total cells.

```julia
create_stretched_mesh(L::Real, n_cells::Int, ::Val{NC}, stretch::Real;
                      T_initial::Real = 300.0,
                      ξ_initial = nothing,
                      cross_section_area::Real = 1.0,
                      stretch_from::Symbol = :top) -> Mesh1D{NC}
```

The extra arguments beyond the uniform case:

| Argument | Default | Meaning |
|---|---|---|
| `stretch` (positional) | — | Geometric stretching factor. `stretch = 1` ⇒ uniform (any value with `\|stretch − 1\| < 1e-10` is treated exactly as uniform). `stretch > 1` makes successive cells larger by that factor along the stretching direction. Must be `≥ 1` — a sub-unity value throws an `ArgumentError`; choose which end is refined with `stretch_from`, never with `stretch < 1`. |
| `stretch_from` (keyword) | `:top` | Which surface the **fine** end is anchored to: `:top` puts the smallest cells at `z = L` (the exposed surface); `:bottom` puts the smallest cells at `z = 0`. |

The cell thicknesses follow a geometric series. With `stretch_from = :top` the smallest cell is
adjacent to the exposed surface and the spacing grows geometrically away from it; the first cell
thickness `Δz_0` is chosen so the cells sum exactly to `L`:

```
Δz_0 = L · (stretch − 1) / (stretch^n − 1)
```

Fine cells at the exposed surface, coarse toward the cold back:

```julia
# 50 cells over 1 cm, smallest cells at the top (exposed) surface
mesh = create_stretched_mesh(0.01, 50, Val(2), 1.1; stretch_from = :top)
```

A modest factor (`1.05`–`1.15`) is usually enough. Very aggressive stretching produces large
cell-to-cell size jumps, which can hurt accuracy and conditioning; the conduction and transport
schemes use distance-weighted harmonic-mean face properties precisely so that graded meshes
stay second-order accurate (Technical Reference §7 and §13), but extreme ratios still degrade
the gradient reconstruction. If you stretch hard, sanity-check the cell-size jump with
`Pyrolysis.Adaptivity.mesh_quality_summary` (§6.5.4).

> **Tip.** When in doubt, start uniform, run a convergence study (§6.3), and only switch to a
> stretched mesh if a uniform mesh with the same near-surface resolution is too expensive.

---

## 6.3 Choosing `n_cells`: resolution versus cost

The number of cells controls both accuracy and run time. There is no universal "correct"
`n_cells`; it depends on how thin the reacting/thermal layer is relative to the sample. Some
practical guidance:

- **Resolve the thermal/reaction front.** For high-flux cone tests the front can be a fraction
  of a millimetre thick. As a rule of thumb you want several cells across that layer. The cone
  example above uses `72` cells over `5.8 mm` (≈ `0.08 mm`/cell). The Douglas-fir convergence
  example (`examples/convergence_douglas_fir.jl`) sweeps `n_cells` and reports cell size in mm,
  peak mass-loss rate, and time-to-ignition so you can see when the answer stops moving.
- **Cost scales with the state size.** The solved state has length
  `n · (1 + NC)` for a fixed (Eulerian) mesh, plus `n_nodes` more for ALE and `n` more for a
  variable cross-section (see §6.6). The stiff implicit solver factorizes a Jacobian whose
  size grows with that state, so doubling `n_cells` more than doubles the per-step cost.
- **Run a convergence study.** Solve the same problem at increasing `n_cells` and watch a
  scalar output of interest (peak MLR, ignition time, back-face temperature) converge. Stop
  refining when successive results agree to your tolerance. This is the recommended way to pick
  a production resolution. See Technical Reference §18 for the verification/convergence
  methodology and the `convergence_douglas_fir` example for a worked spatial-convergence sweep.

A minimal convergence-loop sketch (build the same problem at each resolution, extract one
scalar, compare):

```julia
using Pyrolysis

results = Dict{Int,Float64}()
for n in (24, 48, 96, 192)
    mesh = create_uniform_mesh(5.8e-3, n, Val(3); ξ_initial = (1210.0, 0.0, 0.0))
    problem = PyrolysisProblem(
        mesh      = mesh,
        material  = material,                       # your Material{3,...}
        bc_set    = bcs,                            # your BoundaryConditionSet
        T_initial = z -> 300.0,
        ξ_initial = [z -> 1210.0, z -> 0.0, z -> 0.0],
        tspan     = (0.0, 1600.0),
    )
    sol = solve(problem; radiation_model = SURFACE_ABSORPTION, progress = false)
    results[n] = maximum(sol.mass_loss_rate)        # peak MLR  [kg·m⁻²·s⁻¹]
end
```

When the entries of `results` agree to your target tolerance, you have enough cells.

---

## 6.4 Cross-sectional area `A_0` in a 1D model

A 1D through-thickness model has no real lateral extent, so the cross-sectional area is a
**gauge / normalization constant** rather than a physical width. It enters the model only
through cell volumes and face areas:

```
V_i = |z_{i+1} − z_i| · A          A_face = A   (same A on every face in 1D)
```

Because the governing equations are written per unit volume and per unit area, the area cancels
out of intensive results (temperature profiles, MLR per unit area, surface temperature) **as
long as `A` is constant**. The default `cross_section_area = 1.0 m²` is therefore the natural
choice for almost every simulation, and it makes per-area diagnostics read directly. Set it to
something else only if you have a specific reason to track absolute (extensive) quantities such
as total mass or total energy in physical units.

The mesh stores two area fields:

- `cross_section_area` — the **current** area `A(t)` `[m²]`. For a fixed-cross-section problem
  this never changes.
- `cross_section_area_initial` — the fixed **calibration reference** `A_0` `[m²]`, set once at
  mesh creation.

These differ only when you enable a **lateral-shrinkage law** on the material, which makes the
area evolve as `A(t) = A_0 · law(χ̄)` with `χ̄` the column-averaged pyrolysis progress. That is
the "variable cross-section" (`WithChi`) mode; it is triggered by attaching a
`lateral_shrinkage_law` to the `Material` (covered in Chapter 7 and Technical Reference §10).
When the area varies, mass-loss-rate diagnostics are normalized to the original area `A_0`
(the `MLR_gauge` extra), so they remain comparable to fixed-area runs.

You pass the area at mesh-creation time:

```julia
# 100 cm² sample face, otherwise identical to the default-area mesh
mesh = create_uniform_mesh(0.01, 100, Val(2); cross_section_area = 1e-2)
```

---

## 6.5 Inspecting a mesh

Once you have a `Mesh1D`, you can query its size and geometry. **Note on imports:** the basic
size accessors and geometry queries are exported by the `Geometry` submodule but are *not* part
of the curated top-level `Pyrolysis` surface. Reach them either fully qualified
(`Pyrolysis.Geometry.cell_size`, etc.) or by bringing in the umbrella namespace with
`using Pyrolysis.Internal`, which re-exports every submodule's public symbols. The two mesh
*creators* (`create_uniform_mesh`, `create_stretched_mesh`) are top-level exported; the mesh
*type* (`Mesh1D`) is exported by the `Geometry` submodule — reach it as
`Pyrolysis.Geometry.Mesh1D` or through the `Pyrolysis.Internal` umbrella — not from the
curated top-level surface.

```julia
using Pyrolysis
using Pyrolysis.Internal   # makes the accessors below available unqualified
```

### 6.5.1 Size accessors

| Call | Returns |
|---|---|
| `n_cells(mesh)` | Number of cells `n`. |
| `n_nodes(mesh)` | Number of nodes (`= n + 1`). |
| `n_faces(mesh)` | Number of faces (`= n + 1`). |

```julia
mesh = create_uniform_mesh(0.01, 100, Val(2))
n_cells(mesh)        # 100
n_nodes(mesh)        # 101
n_faces(mesh)        # 101
```

### 6.5.2 Geometry queries

| Call | Returns |
|---|---|
| `cell_size(mesh, i)` | Thickness `Δz_i = \|z_{i+1} − z_i\|` of cell `i` `[m]`. |
| `domain_bounds(mesh)` | Tuple `(z_min, z_max)` from node positions `[m]`. |

The cached geometry is held in plain `Mesh1D` fields, which you read directly:

| Field | Contents |
|---|---|
| `mesh.cell_centers` | Center positions `z_i^c` of the cells `[m]`. |
| `mesh.cell_volumes` | Cell volumes `V_i = A·Δz_i` `[m³]`. |
| `mesh.face_areas` | Face areas `A` `[m²]` (constant in 1D). |
| `mesh.boundary_faces` | `Dict` from boundary id (`1` = top, `2` = bottom) to face indices. |
| `mesh.cross_section_area` | Current area `A(t)` `[m²]` (`mesh.cross_section_area_initial` holds `A_0`). |

```julia
domain_bounds(mesh)              # (0.0, 0.01)
mesh.cell_centers[1]             # ≈ 5.0e-5   (first cell, near z = 0, the bottom)
mesh.cell_centers[end]           # ≈ 9.95e-3 (last cell, near z = L, the top/exposed surface)
cell_size(mesh, 1)               # 1.0e-4    (uniform Δz)
sum(mesh.cell_volumes)           # 0.01 * 1.0  (length × default area = 0.01 m³)

# Face indices for the two boundaries (id 1 = top/exposed, id 2 = bottom/substrate)
mesh.boundary_faces[1]           # [101]  — top face at z = L
mesh.boundary_faces[2]           # [1]    — bottom face at z = 0
```

For the user-facing boundary tags (`:top`, `:bottom`) that you use when specifying boundary
conditions, see Chapter 5; internally `:top` maps to id `1` and `:bottom` to id `2`.

### 6.5.3 Cell-by-cell positions (for plotting)

To get all cell-center positions as a vector — for plotting a temperature profile against
depth — read `mesh.cell_centers` directly, or, after a solve, use the solution accessor
`get_cell_positions` (exported by the `Solver` submodule; reach it through
`Pyrolysis.Internal`):

```julia
z = mesh.cell_centers                 # Vector{Float64} of length n_cells  [m]

# after solving:
z = get_cell_positions(solution)      # same values (Eulerian), read from solution.problem.mesh
```

### 6.5.4 Mesh-quality metrics (graded meshes)

For stretched meshes you can check grading quality with the metrics in the `Adaptivity`
submodule (exported from it, so also available through `Pyrolysis.Internal`):

| Call | Meaning |
|---|---|
| `Pyrolysis.Adaptivity.compute_aspect_ratios(mesh)` | Per-cell largest size ratio to a neighbour (`1.0` = uniform locally). |
| `Pyrolysis.Adaptivity.mesh_quality_summary(mesh)` | NamedTuple with `n_active_cells`, `min_cell_size`, `max_cell_size`, `max_aspect_ratio`, `max_smoothness`, `all_positive_volume`. |

```julia
using Pyrolysis
mesh = create_stretched_mesh(0.01, 50, Val(2), 1.1; stretch_from = :top)
q = Pyrolysis.Adaptivity.mesh_quality_summary(mesh)
q.max_aspect_ratio   # ≈ 1.1 for stretch = 1.1
q.min_cell_size      # the smallest cell, at the exposed surface
```

A `max_aspect_ratio` close to your `stretch` factor is expected and fine. Values far above ~2
indicate aggressive grading; consider more cells or a gentler stretch.

---

## 6.6 The state-vector layout (for users who index raw state)

Most users never touch the raw ODE state vector — `solve` returns a `PyrolysisSolution` with
named fields (`T`, `ξ`, `z`, `surface_T`, `mass_loss_rate`, …; see Chapter 9). But if you write
a callback, a custom diagnostic, or a sensitivity output function, you need to know how the flat
state vector `u` is packed. The packing is **block-major** and is described by a `StateLayout`,
which is itself derived from the problem's `SimulationMode` (Chapter 7). With `n = n_cells` and
`NC` components, the blocks appear in this order:

| Block | Present when | Indices | Contents |
|---|---|---|---|
| `T` | always | `1 .. n` | Cell temperatures `[K]` |
| `ξ_j` | always | `n + (j−1)·n + (1..n)` for `j = 1..NC` | Mass concentrations `ξ_{j,i}` `[kg·m⁻³]` |
| `z` | ALE only | `n·(1+NC) + (1..n_nodes)` | Node positions `[m]` |
| `χ` | variable cross-section only | trailing `n` entries | Pyrolysis progress `χ_i` (dimensionless) |

The total length is `state_length(layout) = n·(1+NC) [+ n_nodes if ALE] [+ n if WithChi]`.

### 6.6.1 Reading state safely with accessors

Rather than computing offsets by hand, use the layout-aware accessors (top-level exported):

```julia
# `layout` is a StateLayout; `u` is a raw state vector (e.g. from a callback)
T_i  = get_temperature(u, layout, i)        # temperature of cell i  [K]
ξ_ji = get_concentration(u, layout, i, j)   # concentration of comp j in cell i  [kg·m⁻³]

# block views (no copies): Val(:T), Val(:ξ) (with component index), Val(:z), Val(:χ)
T_block  = block_view(u, layout, Val(:T))       # view of all temperatures
ξ1_block = block_view(u, layout, Val(:ξ), 1)    # view of component 1's concentrations
```

ALE-only accessors return zero / empty for non-ALE layouts, so they are safe to call
unconditionally:

```julia
z_i        = get_node_position(u, layout, i)         # node position  [m]  (ALE only)
z_nodes    = get_node_positions(u, layout)           # view of the whole z block (ALE only)
thickness  = get_thickness_from_state(u, layout)     # z_top − z_bottom  [m]; 0 for Eulerian
```

You can build a `StateLayout` directly from the component count, mesh size, and mode:

```julia
layout = StateLayout(NC, n_cells(mesh), n_nodes(mesh), mode)   # mode :: SimulationMode
```

In practice the layout is created for you inside the solver workspace; when you write a callback
you receive it (or can read it from the workspace). See Technical Reference §15 for the full
state-layout descriptor and the residual that consumes it, and Chapter 7 for how
`SimulationMode` is selected.

---

## 6.7 Eulerian versus ALE geometry: when `z` becomes state

There are two ways Pyrolysis.jl handles the mesh during a solve, controlled by the `use_ale`
keyword of `solve` (default `false`):

### 6.7.1 Eulerian (fixed mesh) — the default

With `use_ale = false`, the node positions are **fixed** for the whole simulation. Geometry
(cell centers, volumes, face positions) is constant and is mirrored once when the solver builds
its workspace. The state vector contains only the `T` and `ξ` blocks (plus `χ` if you have a
variable cross-section). This is the right choice for almost all cases: TGA-like idealizations,
cone-calorimeter runs without explicit surface recession, and any problem where you do not need
the mesh to physically move with the receding/swelling material. Surface recession in the
Eulerian picture is represented through the density/concentration fields and optional adaptive
surface-cell depletion (Chapter 10), not by moving nodes. A fixed mesh parks only its top live
cell, so a depletion run stops when a cell beneath the surface empties first; a material that
regresses belongs on the ALE mesh (§10.2.1).

```julia
sol = solve(problem; progress = false)            # Eulerian (use_ale = false, the default)
```

### 6.7.2 ALE (moving mesh) — `use_ale = true`

With `use_ale = true`, the node positions `z_i` become **part of the solved state**. The mesh
moves with the material: every component, pore gas included, rides the mesh (the gas moves
relative to it only through its transport flux), so the mesh contracts as the solid pyrolyzes and shrinks (or
expands for an intumescent/swelling material). The state vector then includes the `z` block
(see §6.6), and the residual computes a mesh velocity `w` cumulatively from the pinned bottom
node upward, advancing nodes as `dz_i/dt = w_i`. The cell volumes are recomputed geometrically
from the moving nodes each evaluation, which keeps the discrete Geometric Conservation Law
satisfied to machine precision (Technical Reference §11).

```julia
sol = solve(problem; use_ale = true, progress = false)   # moving (Lagrangian) mesh
```

In ALE mode the solution carries the node-position history in `solution.z` (shape
`(n_nodes, n_steps)`; empty for Eulerian runs), so you can plot mesh motion and surface
recession over time. Use ALE when surface regression or volume change is a first-class output
you want resolved geometrically; otherwise prefer the cheaper Eulerian default. ALE is also a
prerequisite for the variable-cross-section (lateral shrinkage) physics, since the column
strain and the area law are coupled through the mesh-velocity computation.

> **Caveat — vanishing cross-section / mesh tangling.** ALE meshes can in principle tangle
> (nodes crossing) or shrink toward zero thickness for very aggressive recession. The solver
> guards against this with the `min_thickness` keyword of `solve` (default `1e-6 m`), which
> terminates the run, and with depleted-cell parking (`handle_depletion = true`, Chapter 10);
> if a run terminates early or reports a thickness floor, see the
> troubleshooting notes in Chapter 15 and the ALE quality discussion in Technical Reference §11.

### 6.7.3 Geometry derived from state (advanced)

When you do need geometry from a raw ALE state vector (for a custom diagnostic), the
state-layout helpers compute it without touching the mutable mesh:

```julia
thickness = get_thickness_from_state(u, layout)   # current z_top − z_bottom  [m]
```

`compute_volumes_from_state(u, layout, A)` (in `Pyrolysis.Problem`) returns the cell volumes
`A·|z_{i+1} − z_i|` directly from the `z` block, and `update_state_geometry!` fills the
residual's `StateCache` (volumes, centers, face positions) from the state; the latter is what
the residual uses internally (Technical Reference §13). For variable-cross-section runs,
`cross_section_area_from_state` and `column_chi_bar_from_state` (in `Pyrolysis.Residual`)
reconstruct `A(t)` and `χ̄` from the state — advanced hooks you will rarely need directly.

---

## 6.8 Summary and cross-references

- Build meshes with `create_uniform_mesh(L, n_cells, Val(NC); ...)` or
  `create_stretched_mesh(L, n_cells, Val(NC), stretch; stretch_from = :top, ...)`. Both return
  `Mesh1D{NC}`; `NC` must match your `Material`.
- `z = 0` is the bottom/substrate, `z = L` is the top/exposed surface; cells number from bottom
  to top.
- Choose `n_cells` by resolving the near-surface front; confirm with a convergence study
  (Technical Reference §18; `convergence_douglas_fir` example). Stretch toward `:top` to
  concentrate cells at the exposed surface.
- `cross_section_area` (`A_0`, default `1.0 m²`) is a gauge constant; it cancels from intensive
  results and only matters as `A(t) = A_0·law(χ̄)` when a lateral-shrinkage law is active
  (Technical Reference §10).
- Inspect meshes with `n_cells`/`n_nodes`/`n_faces`, `cell_size` and `domain_bounds` —
  qualify with `Pyrolysis.Geometry.` or bring them in via `using Pyrolysis.Internal` —
  and read the cached geometry from the `Mesh1D` fields (`cell_centers`, `cell_volumes`,
  `face_areas`, `boundary_faces`, `cross_section_area`).
- The state vector is block-major `[T | ξ_1..ξ_NC | z (ALE) | χ (WithChi)]`; use
  `get_temperature`, `get_concentration`, `block_view`, and the ALE accessors instead of raw
  offsets.
- `use_ale = false` (default) keeps the mesh fixed (Eulerian); `use_ale = true` makes node
  positions part of the state so the mesh moves with the material (Technical Reference §11).

For the finite-volume discretization built on this mesh, see Technical Reference §13; for
moving-mesh mechanics and the Geometric Conservation Law, Technical Reference §11; for how the
mesh plugs into a problem and how modes are selected, Chapter 7; for boundary conditions on the
top/bottom faces, Chapter 5; and for surface depletion and the mesh-adaptation building blocks,
Chapter 10.
