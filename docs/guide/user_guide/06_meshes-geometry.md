# 6. Meshes and Geometry

The mesh is the spatial backbone of every Pyrolysis.jl simulation: it discretizes the
through-thickness coordinate `z` into the finite-volume cells whose component masses and
energy the solver integrates. A full solve runs on a **line**: a one-dimensional mesh of
cells stacked along `z`, together with a coordinate system and a cross-section measure
that turn node positions into face areas and cell volumes. The mesh holds geometry and
nothing else. It carries no component count and no initial state: the problem's
`T_initial` and `ξ_initial` are evaluated at the cell centres when `setup` builds the
problem the integrator runs (Chapter 7). This chapter shows how to create meshes (uniform,
stretched and from a list of nodes), how to choose the resolution, what the coordinate
system and the cross section mean in a one-dimensional model, how a material's lateral
shrinkage law makes the cross section vary, how to inspect a mesh, where the mesh appears
in the solved state, and how the nodes move when the mesh moves with the solid. For the
finite-volume discretization see Technical Reference §13; for the moving-mesh (ALE)
formulation, Technical Reference §11.

> **Coordinate convention (read this first).** Throughout Pyrolysis.jl, `z = 0` is the
> **back** face (boundary tag `:bottom`, the substrate side) and `z = L` is the **exposed
> surface** (tag `:top`). Heat enters from the top; a receding material recedes toward the
> back, whose node is pinned. Cell indices increase from the back (cell 1, next to `z = 0`)
> to the surface (cell `n`, next to `z = L`). The normal of an interior face points toward
> increasing `z`, and the normal of a boundary face points out of the domain.

---

## 6.1 The mesh in one dimension

A line mesh consists of:

- **Nodes**, the points that bound the cells. A mesh of `n` cells has `n + 1` nodes at
  strictly increasing positions; node 1 sits at `z = 0`, node `n + 1` at `z = L`.
- **Cells**, the finite-volume control volumes. Cell `i` lies between nodes `i` and `i + 1`,
  and its centre is the midpoint of the two. The solver carries, per cell, the mass of every
  component and an energy state (the temperature by default); the concentrations
  `ξ_{j,i} = m_{j,i}/V_i` and the temperature are derived from them in each evaluation.
- **Faces**, the points through which fluxes pass. In one dimension a face sits at a node,
  so there are also `n + 1` faces, face `f` at node `f`. Each face records its two cells,
  `(owner, neighbour)`: interior face `f` has owner `f − 1` and neighbour `f`; boundary
  face 1 is owned by cell 1 and tagged `:bottom`; boundary face `n + 1` is owned by cell `n`
  and tagged `:top`. Conduction and gas fluxes are evaluated at the faces, the boundary
  conditions on the tagged faces.

The mesh type is `Mesh{1}`. It holds two things:

- `mesh.topology`, a `MeshTopology`: the faces and cells above, the node count, the
  boundary tags (tag → face ids), the coordinate system, the cross-section measure `A₀`,
  the **pinned** nodes (whose velocity is zero; node 1 on every mesh the constructors build)
  and the **Lagrangian** nodes (which always move with the solid; none by default);
- `mesh.positions`, the reference node positions (one `SVector{1}` per node): the positions
  at which the initial conditions are evaluated and from which a moving mesh starts.

Every metric (face areas, face normals, cell volumes, cell centres) is computed from the
node positions by `Pyrolysis.Meshes.geometry!`, according to the coordinate system (§6.4).
For the default planar coordinates the volume of cell `i` is

```
V_i = A · (z_{i+1} − z_i)
```

with `A` the cross section. You never construct `Mesh` directly; you call one of the three
constructors of §6.2. `solve` and `setup` never modify the mesh you pass: they build their
own topology with the nodes the solve moves (§6.7) and carry the positions in the state.
See Technical Reference §13 for the cell, face and node layout and the operators built on it.

---

## 6.2 Creating a mesh

Three mesh constructors are exported at the top level (`using Pyrolysis` makes them
available):

| Function | Purpose |
|---|---|
| `create_uniform_mesh(L, n_cells; ...)` | Equal cell sizes `Δz = L / n_cells`. |
| `create_stretched_mesh(L, n_cells, stretch; ...)` | Geometrically graded cells, finest at one end. |
| `create_mesh_from_nodes(z; ...)` | Nodes at positions you give, from zero. |

All three return a `Mesh{1}` with the tags `:bottom` (at `z = 0`) and `:top` (at the last
node), node 1 pinned and no Lagrangian nodes, and all three take the same two keywords:

| Keyword | Default | Meaning |
|---|---|---|
| `coordinates` | `Planar()` | The coordinate system: `Planar()`, `Cylindrical()` or `Spherical()` (§6.4.1). |
| `cross_section` | `1.0`, `1.0`, `4π` | The cross-section measure `A₀`: the face area of a planar column [m²], the axial length of a cylinder [m], the solid angle of a sphere [sr] (§6.4.2). Must be finite and positive. |

### 6.2.1 Uniform meshes

```julia
create_uniform_mesh(L::Real, n_cells::Integer;
                    coordinates = Planar(),
                    cross_section = <default for the coordinates>) -> Mesh{1}
```

Positional arguments:

- `L`, the thickness of the domain [m]: the nodes run from `z = 0` to `z = L`.
- `n_cells`, the number of cells, at least one.

A 1 cm slab with 100 cells:

```julia
using Pyrolysis

mesh = create_uniform_mesh(0.01, 100)    # L = 1 cm, 100 cells
```

A typical cone-calorimeter sample, 5.8 mm with 72 cells (0.08 mm per cell):

```julia
sample_thickness = 5.8e-3      # [m]
n_cells          = 72

mesh = create_uniform_mesh(sample_thickness, n_cells)
```

The node positions are `z_i = (i − 1)·Δz` with `Δz = L / n_cells`, so cell `i` occupies
`[(i − 1)Δz, iΔz]` and, in planar coordinates, every cell has the volume `V_i = A·Δz`. The
composition and temperature of the sample are not part of the mesh: they are the problem's
`ξ_initial` and `T_initial` (Chapter 7; Chapter 5, §5.10).

### 6.2.2 Stretched (graded) meshes

In a cone-calorimeter or fire-exposure problem the steepest temperature and reaction-rate
gradients are at the exposed surface, where the pyrolysis front lives. A graded mesh puts
small cells there and large cells deep in the cold substrate, giving the same accuracy for
fewer cells.

```julia
create_stretched_mesh(L::Real, n_cells::Integer, stretch::Real;
                      stretch_from::Symbol = :top,
                      coordinates = Planar(),
                      cross_section = <default for the coordinates>) -> Mesh{1}
```

The arguments beyond the uniform case:

| Argument | Default | Meaning |
|---|---|---|
| `stretch` (positional) | — | The ratio of successive cell sizes away from the refined end. `stretch = 1` is uniform (any value within `1e-10` of one is treated as exactly uniform); `stretch > 1` makes each cell that factor larger than its neighbour on the refined side. A value below one throws an `ArgumentError`: choose the refined end with `stretch_from`, never with `stretch < 1`. |
| `stretch_from` (keyword) | `:top` | The end that gets the finest cell: `:top`, the exposed surface at `z = L`, or `:bottom`, the back at `z = 0`. Any other value throws an `ArgumentError`. |

The cell sizes follow a geometric series from the refined end, the finest cell's size
chosen so that the cells sum exactly to `L`:

```
Δz_0 = L · (stretch − 1) / (stretch^n − 1),     Δz_k = Δz_0 · stretch^k
```

The end nodes are placed at exactly `0` and `L`. Fine cells at the exposed surface,
coarser toward the back:

```julia
# 50 cells over 1 cm, the finest (about 8.6 µm) at the exposed surface
mesh = create_stretched_mesh(0.01, 50, 1.1)

# the same grading with the finest cell at the back face
mesh_back = create_stretched_mesh(0.01, 50, 1.1; stretch_from = :bottom)
```

A modest ratio (`1.05` to `1.15`) is usually enough. Strong grading produces large
cell-to-cell size jumps, which degrade the accuracy of the face gradients; the conduction
and gas-transport face coefficients are distance-weighted harmonic means, which keeps graded
meshes consistent (Technical Reference §7 and §13), but they do not remove the truncation
error of a large jump. Note how small the surface cell becomes: at `stretch = 1.1` and 50
cells it is about a hundredth of the uniform cell.

> **Tip.** When in doubt, start uniform, run a convergence study (§6.3), and switch to a
> stretched mesh only if a uniform mesh with the same near-surface resolution is too
> expensive.

### 6.2.3 Meshes from a list of nodes

```julia
create_mesh_from_nodes(z::AbstractVector{<:Real};
                       coordinates = Planar(),
                       cross_section = <default for the coordinates>) -> Mesh{1}
```

builds a line with nodes at the positions `z` [m], which must start at zero, be finite and
increase strictly (at least two nodes); the cells lie between consecutive nodes and are
sized independently. Use it to put a node on every interface of a layered sample, or to
refine one region by hand. A 6 mm substrate under a 1 mm coating, with a node on the
interface at 6 mm:

```julia
z = [range(0.0, 6e-3; length = 25); range(6e-3, 7e-3; length = 11)[2:end]]
mesh = create_mesh_from_nodes(z)        # 24 substrate cells, 10 coating cells
```

The initial conditions then give each layer its composition; for a material whose first
component is the substrate and whose second is the coating,
`ξ_initial = [z -> z < 6e-3 ? 2000.0 : 0.0, z -> z < 6e-3 ? 0.0 : 900.0]`. Since they are
evaluated at the cell centres, a node on the interface keeps every cell within one layer. When the burnout correction moves a layered mesh, the interface nodes are held to the
solid (§6.7.2).

The constructors of this section are built on `Pyrolysis.Meshes.line_mesh(x; coordinates,
cross_section, bottom = :bottom, top = :top, fixed = false, lagrangian = Int[])`, which
accepts any strictly increasing positions (a hollow cylinder from an inner radius, for
instance) and lets you name the boundary tags, pin every node (`fixed`) or declare Lagrangian
nodes. A solve needs the tags `:top` and `:bottom` (§7.9).

---

## 6.3 Choosing `n_cells`: resolution versus cost

The number of cells controls both accuracy and run time. There is no universal correct
`n_cells`; it depends on how thin the reacting and thermal layers are relative to the
sample.

- **Resolve the thermal and reaction front.** Under a high cone flux the front can be a
  fraction of a millimetre thick, and you want several cells across it. The cone example
  above uses 72 cells over 5.8 mm. The Douglas-fir convergence example
  (`examples/convergence_douglas_fir.jl`) solves its deck at 15 to 480 cells and tabulates
  the peak mass-loss rate, the ignition time and the solve time against the cell size, so
  you can see where the answer stops moving.
- **Cost grows with the state.** Each cell carries `NC + 1` unknowns (its `NC` component
  masses and its energy state); a moving mesh adds two per node (position and velocity),
  in-depth radiation one per face, a variable cross section two per cell and two for the
  column (§6.6, §7.8). The state is ordered so that the Jacobian is banded, its bandwidth
  set by the component count rather than the cell count (a variable cross section adds two
  dense rows), so the work of each step grows roughly in proportion to the cells. Time a
  short run at two resolutions rather than guessing.
- **Burnout needs its own resolution.** A solve that deactivates cells through a burnout is
  verified at stated cell counts across the vanishing layer; Chapter 10 gives them.
- **Run a convergence study.** Solve the same problem at increasing `n_cells` and watch a
  scalar output of interest (peak mass-loss rate, ignition time, back-face temperature)
  converge. Stop refining when successive results agree to your tolerance. This is the
  recommended way to pick a production resolution; Technical Reference §18 describes the
  verification and convergence methodology.

A complete convergence loop: the same problem at each resolution, two scalars per solve.

```julia
using Pyrolysis

virgin = SolidComponent(:virgin; ρ = 1190.0, c = 1500.0, k = 0.20, λ = 1e-6)
gas    = GaseousComponent(:MMA; M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)
pmma   = Material(name = :PMMA, components = (virgin, gas),
                  reactions = (Reaction(:pyrolysis, :virgin => :MMA;
                                        A = 8.5e12, E = 188e3, h = 870e3),))
bcs = BoundaryConditionSet(
    thermal = (top = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300.0),
               bottom = AdiabaticBC()),
    mass    = (top = DiffusiveMassBC(0.01), bottom = ImpermeableBC()),
)

for n in (25, 50, 100, 200)
    problem = PyrolysisProblem(
        mesh      = create_uniform_mesh(0.01, n),
        material  = pmma,
        bc_set    = bcs,
        tspan     = (0.0, 200.0),
    )
    sol = solve(problem; use_ale = true)
    lost = 1 - sol.extras.total_mass[end] / sol.extras.total_mass[1]
    println("n = $n: peak MLR = $(maximum(sol.mass_loss_rate)) kg/(m²·s), ",
            "mass lost = $(round(100lost, digits = 3)) %, T_s = $(sol.surface_T[end]) K")
end
```

When successive lines agree to your target tolerance, you have enough cells.

---

## 6.4 Coordinate systems and the cross section

### 6.4.1 Coordinate systems

The coordinate system is a property of the mesh, chosen with the `coordinates` keyword of
the constructor. It decides how face areas and cell volumes follow from the node positions:

| `coordinates` | Position | Face area at `r` | Volume between `r₋` and `r₊` | `cross_section` | Default |
|---|---|---|---|---|---|
| `Planar()` | depth `z` | `A` | `A (r₊ − r₋)` | face area [m²] | `1.0` |
| `Cylindrical()` | radius `r` | `2π r A` | `π (r₊² − r₋²) A` | axial length [m] | `1.0` |
| `Spherical()` | radius `r` | `A r²` | `A (r₊³ − r₋³)/3` | solid angle [sr] | `4π` |

In cylindrical and spherical coordinates the positions are radii: a mesh from
`create_uniform_mesh(L, n; coordinates = Cylindrical())` runs from the axis (the `:bottom`
tag, a face of zero area) to the exposed surface at radius `L` (`:top`), so the cylinder or
sphere is heated from outside and recedes toward its axis or centre.

```julia
rod  = create_uniform_mesh(0.01, 40; coordinates = Cylindrical(), cross_section = 0.1)
ball = create_uniform_mesh(0.01, 40; coordinates = Spherical())   # the whole sphere, 4π
```

Conduction, reactions, gas transport, surface absorption and the moving mesh are
implemented in every coordinate system. Deactivation is not refused on a radial mesh, but
the burnout path of Chapter 10 is verified on planar decks only. Two models need planar
coordinates and are refused with an `ArgumentError` otherwise: in-depth radiation
(`radiation_model = BEER_LAMBERT`) and the variable cross section (§6.4.3).

### 6.4.2 The cross-section measure `A₀`

A planar one-dimensional model has no real lateral extent, so its cross-sectional area is a
**gauge constant** rather than a physical width. It enters only through the face areas and
cell volumes:

```
V_i = A (z_{i+1} − z_i)          A_face = A   (the same on every face)
```

Because the balances are written per unit volume and the fluxes per unit area, `A` cancels
from the intensive results (temperature profiles, mass-loss rate per unit area, surface
temperature, the per-area surface fluxes). The default `cross_section = 1.0` m² is the
natural choice for almost every planar simulation. The extensive outputs scale with it: the
cell masses `sol.m` [kg], the live inventories `sol.extras.live_mass` [kg] and
`sol.extras.live_energy` [J], and the conservation audit's terms. Set a different area only
when you want those in the units of a physical sample:

```julia
# A 100 cm² sample face; the per-area results are unchanged
mesh = create_uniform_mesh(0.01, 100; cross_section = 1e-2)
```

In cylindrical and spherical coordinates the measure is not a gauge in the same sense: the
face area varies with the radius, and `A` is the axial length or the solid angle that scales
it. The per-area outputs divide by the current area of the `:top` faces.

### 6.4.3 The variable cross section (lateral shrinkage)

A charring sample in a cone calorimeter shrinks (or, if intumescent, swells) laterally as well
as through its thickness. A one-dimensional model represents this with a **variable cross
section**: the area follows the column's pyrolysis progress,

```
A(t) = A₀ · law(χ̄),          χ̄ = Σ_i G_i / Σ_i M_i
```

with `A₀` the mesh's `cross_section`, `G_i` the cumulative dry pyrolysis gas cell `i` has
produced (every gas the reactions make except water vapour) and `M_i` the cell's initial dry
solid mass. The law is not a mesh setting: it is the material's `lateral_shrinkage_law`
(Chapter 7, §7.6), and `setup` adds the cross-section rows whenever the material carries one.
The law must be a smooth function (its derivative is taken by automatic differentiation),
positive and finite on `[−0.001, 1.001]`, and equal to one at zero conversion; a variable
cross section needs a moving, planar mesh (`use_ale = true`).

The area history is reported in `sol.extras.cross_section_area_history` [m²]. With a
variable area, `sol.mass_loss_rate` is per unit **current** area of the `:top` faces and
`sol.extras.MLR_gauge` per unit **initial** area, which keeps the latter comparable with a
fixed-area run and with a load-cell measurement on the original sample face. See Technical
Reference §10 for the model.

---

## 6.5 Inspecting a mesh

The mesh constructors are exported at the top level; the mesh's own functions live in the
`Pyrolysis.Meshes` submodule and are not part of the curated top-level surface. Import
them by name from `Pyrolysis.Meshes`, or bring in every submodule's exports with
`using Pyrolysis.Internal` (Chapter 1, §1.9). They take the topology, `mesh.topology`.

### 6.5.1 Size and tags

| Call | Returns |
|---|---|
| `n_cells(topo)` | Number of cells `n`. |
| `n_nodes(topo)` | Number of nodes (`n + 1`). |
| `n_faces(topo)` | Number of faces (`n + 1`). |
| `boundary_tags(topo)` | The boundary tags, sorted: `[:bottom, :top]`. |
| `boundary_faces(topo, tag)` | The face ids carrying `tag` (an `ArgumentError` for an unknown tag). |
| `is_fixed(topo)` | Whether every node is pinned. |

```julia
using Pyrolysis
using Pyrolysis.Meshes: n_cells, n_nodes, n_faces, boundary_tags, boundary_faces

mesh = create_uniform_mesh(0.01, 100)
topo = mesh.topology
n_cells(topo), n_nodes(topo), n_faces(topo)    # (100, 101, 101)
boundary_tags(topo)                            # [:bottom, :top]
boundary_faces(topo, :top)                     # [101], the face at z = L
boundary_faces(topo, :bottom)                  # [1], the face at z = 0
topo.coordinates, topo.cross_section           # (Planar(), 1.0)
topo.pinned                                    # [1], the back node
```

The boundary conditions are keyed by these tags (Chapter 5): the `top` entries of a
`BoundaryConditionSet` apply on the faces tagged `:top`, the `bottom` entries on `:bottom`.

### 6.5.2 Node positions, cell centres and volumes

```julia
using Pyrolysis.Problem: cell_centers, cell_volumes

z  = first.(mesh.positions)      # node positions [m]: z[1] = 0.0, z[end] = 0.01
xc = cell_centers(mesh)          # cell centres [m]: 5.0e-5 … 0.00995
V  = cell_volumes(mesh)          # cell volumes [m³]; sum(V) == 0.01 × 1.0
```

`cell_centers` returns the points at which the problem's initial conditions are evaluated:
the midpoints of the cells at the reference positions.

### 6.5.3 Every metric: `MeshGeometry`

The full set of metrics at one set of node positions is a `MeshGeometry`, filled by
`geometry!` (the only function that writes one):

```julia
using Pyrolysis.Meshes: MeshGeometry, geometry!

g = MeshGeometry{Float64}(topo)
geometry!(g, first.(mesh.positions), topo.cross_section, topo)

g.cell_volume       # V_i [m³]
g.cell_center       # cell centres, one SVector{1} per cell [m]
g.face_area         # face areas [m²]
g.face_normal       # unit normals: −1 on the bottom face, +1 elsewhere
g.face_center       # face positions [m]
```

The element type follows the positions, so the same function evaluates the geometry of a
dual-number state; this is how the Jacobian differentiates through a moving mesh. To get
depth profiles from a solution, use the solution's own positions instead (§6.7.3), which
are those of the solve at each saved time.

---

## 6.6 Where the mesh lives in the state

Most users never touch the raw state vector: `solve` returns a `PyrolysisSolution` with
named fields (Chapter 9). For a custom loop over `setup` and `integrate`, or a callback that
reads the integrator's state, it helps to know how the mesh enters the state. The state is a
sequence of node blocks and cell blocks: on a line, node 1, cell 1, node 2, cell 2, …,
cell `n`, node `n + 1`, so every row couples only to its neighbours and the Jacobian is
banded.

| Block | Holds | Present when |
|---|---|---|
| node `k` | position `x_k` [m], velocity `w_k` [m/s] | moving mesh (`use_ale = true`) |
| node `k` | intensity `I_k` of the face at that node [W/m²] | `radiation_model = BEER_LAMBERT` |
| cell `i` | masses `m_{i,1} … m_{i,NC}` [kg], energy state (`T_i` [K] or `H_i` [J]) | always |
| cell `i` | conversion history `G_i`, `M_i` [kg] | variable cross section |
| global | cross section `A` [m²], lateral strain rate `θ_A` [1/s] | variable cross section |

The node velocities, the face intensities, `A` and `θ_A` are algebraic unknowns (rows of the
system without a time derivative); the others are differential. A fixed mesh carries no node
unknowns at all: its geometry is constant. Never compute an index by hand; the layout of a
set-up problem gives every index through functions of `Pyrolysis.States` (reachable through
`Pyrolysis.Internal`):

```julia
using Pyrolysis.Internal: mass_index, energy_index, position_index, velocity_index,
                          state_length

cp = setup(problem; use_ale = true)          # `problem` as in §6.3, 50 cells
layout = cp.def.layout
state_length(layout)                         # 252 = 50 × (2 + 1) + 51 × 2
cp.u0[energy_index(layout, 50)]              # temperature of the surface cell [K]
cp.u0[mass_index(layout, 50, 1)]             # its mass of component 1 [kg]
cp.u0[position_index(layout, 51)]            # the surface node's position [m]
cp.u0[velocity_index(layout, 51)]            # and its velocity [m/s]
```

Chapter 11 describes the layout, the residual and the other building blocks `Internal`
exposes; Technical Reference §13 the state mapping from these unknowns to concentrations,
temperatures and volumes.

---

## 6.7 Fixed and moving meshes

The `use_ale` keyword of `solve` (default `false`) decides whether the nodes move.

### 6.7.1 Fixed mesh: the default

With `use_ale = false`, every node is pinned for the whole solve. The geometry (cell
centres, volumes, face positions) is constant, and the state holds only the cell blocks (and
the face intensities under `BEER_LAMBERT`). This is the right choice for problems whose
thickness does not change appreciably: TGA-like idealizations, inert or weakly reacting
slabs, and conduction studies. On a fixed mesh a reacting material that loses mass becomes
less dense in place; its surface does not recede.

```julia
sol = solve(problem)                         # fixed mesh (use_ale = false)
```

A fixed mesh can deactivate burned-out cells (`deactivation = true`, Chapter 10): from the
exposed side, each cell whose condensed mass falls below the threshold goes, its condensed
inventories pass to the live cell beneath it, and the gas it held is vented and recorded in
the deactivation event, since the receiver's pores could hold it only at a higher pressure.
A cell below the threshold beneath a live cell ends the solve with outcome `:interior`. A
material that regresses belongs on the moving mesh.

### 6.7.2 Moving mesh: `use_ale = true`

With `use_ale = true`, the node positions and velocities become part of the solved state.
The back node stays pinned at rest, and every other node moves with the velocity

```
w = v_s + c
```

the solid velocity `v_s` plus a correction `c`. The solid velocity comes from the cells'
volume balances: one algebraic row per cell, on the velocity of its node farther from the
pinned node, states that the solid velocity across the cell's faces accounts for the volume
change of its reactions (the strain rate `θ`, Technical Reference §10), so the rows are
swept from the back node outward. The positions follow `dx_k/dt = w_k`, and the cell volumes
are always computed from the positions, so the discrete geometric conservation law holds by
construction (Technical Reference §11). The velocities are consistent with the rows at the
start of the solve and after every deactivation event.

Without a correction (`c = 0`) the mesh moves with the solid: every component, pore gas
included, rides the mesh, and the gas moves relative to it only through its transport flux.
The mesh contracts as a pyrolyzing solid shrinks and expands for a swelling or intumescent
one. A correction moves the mesh relative to the solid, and the inventories then cross the
faces at the relative velocity `−c`. The correction a solve uses is decided by `setup`: on a
moving mesh with `deactivation = true` and a layer that can burn out, the equidistribution
correction, which thins a vanishing layer's cells together so that the layer leaves in one
deactivation event; otherwise none (§7.3.3, Chapter 10). Under the equidistribution
correction the interior nodes between cells that start with different condensed components,
the layer interfaces, are Lagrangian nodes that always move with the solid; the `lagrangian`
keyword of `solve` names them explicitly.

```julia
sol = solve(problem; use_ale = true)         # the mesh follows the solid
```

The solve ends with return code `:Terminated` and outcome `:min_thickness` when the
thickness, from the back node to the outermost live node, falls below `min_thickness`
(default `1e-6` m). Use a moving mesh when surface regression or volume change is an output
you want resolved geometrically, and whenever a material recedes; a variable cross section
requires it (§6.4.3).

### 6.7.3 Node positions in the solution

A moving-mesh solution records the geometry at every saved time:

| Field or accessor | Contents |
|---|---|
| `sol.z` | Node positions, `n_nodes × n_times` [m]; `NaN` beyond the outermost live node. `0 × 0` on a fixed mesh. |
| `get_profile(sol, :z, k)` | The node positions at saved time `k`; on a fixed mesh the mesh's positions. |
| `sol.cell_centers` | Cell centres, `n_cells × n_times` [m]; `NaN` in a deactivated cell. |
| `get_cell_positions(sol; at = k)` | The cell centres at saved time `k` (default the last); `at = 1` gives the initial centres. |
| `sol.correction_velocity` | The node corrections `c = w − v_s`, `n_nodes × n_times` [m/s]; zero without a correction, `0 × 0` on a fixed mesh. |
| `sol.extras.thickness_history` | The distance from the first node to the outermost live one [m]. |
| `sol.extras.min_cell_thickness` | The thickness of the thinnest live cell [m]. |

```julia
sol = solve(problem; use_ale = true)
z_end  = sol.z[:, end]                    # node positions at the last saved time [m]
x_end  = get_cell_positions(sol)          # cell centres there, for a T(z) profile
x_0    = get_cell_positions(sol; at = 1)  # the initial centres
L_t    = sol.extras.thickness_history     # the sample thickness over time [m]
fields = sol(50.0)                        # the fields at t = 50 s, positions included
fields.z[end]                             # the surface position at 50 s [m]
```

Plot a temperature profile against `get_cell_positions(sol; at = k)` with `sol.T[:, k]`, so
the depth axis spans the thickness the sample had at that time. Chapter 9 describes the
solution object and its interpolation.

---

## 6.8 Meshes of more than one dimension

`Pyrolysis.Meshes` also builds planar meshes of quadrilaterals and hexahedra:
`box_mesh(xs, ys[, zs]; transform, pinned = (:bottom,), fixed, lagrangian)` on the tensor
product of strictly increasing node coordinates, each node optionally moved by
`transform(x::SVector) -> SVector`, with the sides tagged `:left`/`:right` along x,
`:bottom`/`:top` along the last axis and, in 3D, `:front`/`:back` along y; and
`cell_mesh(positions, cells; tag, pinned, fixed, lagrangian)` from node positions and each
cell's corner nodes, the faces found by matching node sets. They return a `Mesh{2}` or
`Mesh{3}` whose metrics `geometry!` computes exactly for the multilinear cells.

```julia
using Pyrolysis.Meshes: box_mesh, boundary_tags

slab = box_mesh(range(0.0, 0.02; length = 5), range(0.0, 0.01; length = 4))
boundary_tags(slab.topology)        # [:bottom, :left, :right, :top]
```

These meshes are building blocks at the level of the residual, not a supported path for
full solves. `setup` accepts them for conduction (with an anisotropic conductivity tensor),
reactions and per-tag thermal conditions (a `BoundaryConditionSet` with a thermal entry for
every tag), on a fixed mesh or a moving one closed by one multiplier per cell (or, for a
single column of cells, by the balance rows); cells leave such a mesh by element death
(`Pyrolysis.Deactivation.element_death!`) applied to the set-up layout. `setup` refuses what
is a line's alone: gas transport, in-depth radiation, the variable cross section, the
equidistribution correction and threshold deactivation. A full simulation runs on a line mesh
from §6.2. Technical Reference §13 describes the multi-dimensional geometry and the closure.

---

## 6.9 Summary and cross-references

- Build meshes with `create_uniform_mesh(L, n_cells; ...)`,
  `create_stretched_mesh(L, n_cells, stretch; stretch_from = :top, ...)` or
  `create_mesh_from_nodes(z; ...)`. Each returns a `Mesh{1}`: a topology and the reference
  node positions, with no component count and no initial state.
- `z = 0` is the back face (`:bottom`, node 1, pinned), `z = L` the exposed surface
  (`:top`); cells number from the back to the surface.
- Choose `n_cells` by resolving the near-surface front and confirm with a convergence study
  (§6.3; `examples/convergence_douglas_fir.jl`; Technical Reference §18). Stretch toward
  `:top` to concentrate cells at the exposed surface.
- `coordinates` (`Planar()`, `Cylindrical()`, `Spherical()`) and `cross_section` (`A₀`,
  default `1.0`, `1.0`, `4π`) are mesh properties. In planar geometry `A₀` is a gauge that
  cancels from intensive results; a material's lateral shrinkage law makes it vary as
  `A₀·law(χ̄)` on a moving planar mesh.
- Inspect a mesh through `mesh.topology` with the functions of `Pyrolysis.Meshes`
  (`n_cells`, `n_nodes`, `n_faces`, `boundary_tags`, `boundary_faces`), its positions through
  `mesh.positions`, and its metrics with `MeshGeometry` and `geometry!`.
- The state interleaves node and cell blocks; a moving mesh adds each node's position and
  velocity. Read it through the layout's index functions (Chapter 11).
- `use_ale = false` (default) keeps every node pinned; `use_ale = true` moves the nodes with
  `w = v_s + c`, the correction installed by `setup` for burnout (§7.3.3, Chapter 10).
- Meshes of quadrilaterals and hexahedra exist at the residual level only.

For the finite-volume discretization built on the mesh, see Technical Reference §13; for the
moving mesh and the geometric conservation law, Technical Reference §11; for how the mesh
enters a problem and what `setup` decides, Chapter 7; for the boundary conditions on the
tagged faces, Chapter 5; for deactivation and the correction, Chapter 10.
