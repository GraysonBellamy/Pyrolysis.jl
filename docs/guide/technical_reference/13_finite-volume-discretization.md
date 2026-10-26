# 13. Finite-Volume Spatial Discretization

This chapter derives the cell-centered finite-volume discretization that turns the
continuum balances of §3 (Governing Equations) into the semi-discrete mass-matrix
system integrated in §15 (Temporal Integration, the DAE, and the Colored Jacobian). It
describes the face-based mesh (nodes, faces, cells and boundary tags,
`src/Meshes`), the metrics computed from the node positions, the scatter of face
fluxes into extensive cell rows, the kernels that compute the fluxes and the cell
sources (`src/Kernels`), the order in which the residual runs them
(`src/Assembly`), the face values and gradients the fluxes use, mesh generation,
and the mapping from the state vector to the quantities the kernels read
(`src/States`). Every formula below is the one the code evaluates; where the
implementation adds a guard, a floor or a sign choice, it is stated.

The organising principle is that **every transport process is a flux through a
face, and every cell row is extensive**. The integrated unknowns are the
component masses of every cell, $m_{j,i}$ [kg], and an energy state per cell. A
kernel writes a flux per unit area along each face normal; the assembler adds
$-F_fA_f$ to the face's owner and $+F_fA_f$ to its neighbour. Interior faces
therefore cancel in any sum over cells, mass is conserved exactly by the
semi-discrete system, and the energy ledger closes against the boundary fluxes,
the sources and the deactivation events (§16). The coupling of every row to its
neighbours is local, which gives the banded sparsity the colored Jacobian of §15
exploits.

Conventions (binding, from §2, Nomenclature):

- The spatial coordinate is $z$; $z=0$ is the back (tag `:bottom`) and $z=L$ the
  exposed surface (tag `:top`). Heat enters from the top and the material recedes
  downward.
- $i$ is a cell index, $j$ a component index, $f$ a face index, $n$ a node index
  (and, where it is a count, the number of cells), $r$ a reaction index (overload
  I1). $N_C$ is the number of components.
- A face's normal points from its **owner** cell to its **neighbour** cell, out of
  the domain on a boundary face. $s_{i,f}=+1$ when cell $i$ owns face $f$ and $-1$
  when it is the neighbour. On a face quantity the subscripts $o$ and $n$ denote
  the face's owner and neighbour cells ($d_o$, $T_n$, …).
- A face flux $F_f$ (energy, W m⁻²) or $F_{j,f}$ (mass of component $j$,
  kg m⁻² s⁻¹) is the transport per unit face area along the normal: along $+z$ on
  every interior face and on the `:top` face, along $-z$ (outward) on the
  `:bottom` face.
- Heat of reaction: $h>0$ endothermic, $Q_\mathrm{rxn}=-\sum_r h_r r_r$.
  $h_j(T)$ is the sensible enthalpy of component $j$ on the $T_\mathrm{ref}$ datum.
- A tilde marks a guarded trial-state quantity (§13.10.3): $\tilde T_i$ the trial
  temperature, $\tilde V_i$ the floored cell volume.

---

## 13.1 Cell-Centered Finite-Volume Layout in 1D

![Nonuniform column with cell centres, node faces, outward boundary normals and equal-and-opposite owner/neighbour flux contributions.](../../assets/figures/01_mesh_anatomy.svg)

**Figure 1. Finite-volume mesh anatomy (schematic).** The pinned back is at `z = 0`; heat enters at the exposed face `z = L(t)`. In one dimension each node is a face. An interior flux is oriented from owner to neighbour and contributes equal and opposite extensive rates, so it cancels from the global balance. Widths are deliberately nonuniform. See Technical Reference §13.1 for the face convention and Figure 2 in §13.4 for the framework.

### 13.1.1 Topology: nodes, faces, cells

Pyrolysis.jl stores its primary unknowns, the component masses and the energy
state, per **cell**: the classic cell-centered (collocated) finite-volume
arrangement. A line of $n$ cells has

$$
n_\text{nodes}=n_\text{faces}=n+1 .
$$

- **Node** $n$ sits at position $z_n$. Node 1 ($z=0$) is pinned. On a moving mesh
  (`use_ale = true`) the node positions are differential unknowns and the node
  velocities algebraic ones (§11.4, §13.10); on a fixed mesh every node is pinned
  and the geometry is computed once.
- **Face** $f$ (`Face{1,1}`) sits at node $f$ and carries its nodes, its two cells
  `cells = (owner, neighbour)`, the neighbour `0` on a boundary face, and a tag:
  face 1 has cells `(1, 0)` and tag `:bottom`; interior face $f$ ($2\le f\le n$)
  has cells `(f − 1, f)` and tag `:interior`; face $n+1$ has cells `(n, 0)` and
  tag `:top`.
- **Cell** $i$ (`Cell{1,2,2}`) lies between nodes $i$ and $i+1$ and is bounded by
  faces $i$ and $i+1$.

Every interior face's owner is the cell below it and its neighbour the cell above,
so its normal points along $+z$; the boundary faces' normals point out of the
domain, $-z$ at face 1 and $+z$ at face $n+1$. A positive interior flux therefore
carries its quantity upward, toward the exposed surface; a positive flux on either
boundary face carries it out of the domain.

> **Implementation note.** `Face`, `Cell`, `MeshTopology` and `Mesh`
> (`src/Meshes/topology.jl`). `MeshTopology{D,C,FT,CT}` holds the coordinate
> system, the faces, the cells, the node count, the boundary tags (tag → face
> ids), the pinned nodes, the Lagrangian nodes (layer interfaces that always move
> with the solid, §11.7.5) and the cross-section measure $A_0$. `Mesh{D}` is the
> topology plus the reference positions of the nodes, the positions at which
> `setup` evaluates the problem's `T_initial` and `ξ_initial` (at the cell
> centres, §12.8). The mesh carries no component count, no initial state and no
> per-cell inventory. Its topology never changes during a solve: deactivation
> (§14) changes which cells are live, recorded in the state layout's mask, and
> never resizes an array.

### 13.1.2 Boundary tags

Boundary conditions are keyed by tag; there are no numeric boundary
identifiers. `line_mesh` tags face 1 `:bottom` and face $n+1$ `:top` (other names
through its `bottom` and `top` keywords); `boundary_tags(topology)` returns the
sorted tags and `boundary_faces(topology, tag)` the face ids that carry one, an
unknown tag being an `ArgumentError`. The thermal boundary kernel and the
gas-transport kernel hold one patch per tag (§13.7.5), and `setup` refuses a
boundary-condition set that names a tag the mesh does not have.

A patch is applied on the face that bounds the live cells on its side,
`live_face(topology, active, f)`: $f$ while its owner is live; otherwise, walking
inward across the block of dead cells, the first face whose owner is live; zero
when no live cell remains on that side. After cells are deactivated from the
exposed side, the `:top` condition therefore applies at the upper face of the
outermost live cell, an interior face whose $+z$ normal points out of the live
region (§14.3.1).

### 13.1.3 What a cell-centered scheme means for accuracy

The unknowns are cell inventories, but fluxes are needed at faces, so every flux
needs a reconstruction from cell values to the face. Diffusion uses a face
coefficient (the distance-weighted harmonic mean) and the two-point difference
along the normal (§13.6, §13.7): for smooth fields on a mesh whose non-uniformity
is smooth, the face gradient is second-order. The resistance-in-series mean keeps
the flux continuous across a material interface that lies on a face; V-HC-04
reproduces the steady piecewise-linear profile of a two-layer slab to integrator
precision for conductivity ratios 10:1 and 100:1 (§18). Because the initial
conditions are sampled at the cell centres, a composition interface that does not
lie on a face is moved to the nearest face (there are no mixed cells); a node
placed on each interface (`create_mesh_from_nodes`, §13.9.4) avoids that
quantization.

Two places are first-order by construction. A **boundary face** is reached from
its cell centre by a one-sided difference over the centre-to-face distance
(§13.2.3), so a face value the condition fixes (a prescribed temperature, the
solved surface temperature) enters with a first-order error at the face; this is
the surface-temperature closure error the coarse-mesh verification cases measure
(§13.12). **Darcy advection** is first-order upwind (§13.7.3). The fluxes relative
to a corrected moving mesh are second-order upwind with a smooth slope (§11.6).

---

## 13.2 Geometric Quantities

Every metric is computed from the node positions and the cross-section measure by
`geometry!` (`src/Meshes/geometry.jl`), the only writer of a `MeshGeometry`. The
element type of the metrics follows the positions, so a dual-number state carries
geometry sensitivities and automatic differentiation passes through the geometry.

### 13.2.1 Face areas and cell volumes

With $r$ the node coordinate (the position $z$ in planar geometry, the radius in
cylindrical and spherical geometry) and $A$ the mesh's cross-section measure (the
planar face area [m²], the axial length of a cylinder [m], the solid angle of a
sphere), the coordinate system fixes the face areas and the cell volumes:

| Coordinates | Face area at $r$ | Volume of cell $i$ | $d$ (`geometric_dimension`) | Default $A$ |
|---|---|---|---|---|
| `Planar()` | $A$ | $A\,(r_{i+1}-r_i)$ | 1 | 1 |
| `Cylindrical()` | $2\pi rA$ | $\pi\,(r_{i+1}^2-r_i^2)\,A$ | 2 | 1 |
| `Spherical()` | $Ar^2$ | $A\,(r_{i+1}^3-r_i^3)/3$ | 3 | $4\pi$ |

In every coordinate system $\partial V_i/\partial r_{i+1}=A_{i+1}$ and
$\partial V_i/\partial r_i=-A_i$: the volume changes by the volume its faces sweep,
which is the discrete Geometric Conservation Law (§11.5.2; `test/meshes/test_meshes.jl`,
"volume derivative by AD equals the face sweep operator"), and every cell satisfies
$V_i=\tfrac1d\sum_f A_f\,(x_f\cdot n_f)$. A cylindrical or spherical line that
starts at the axis or the centre has a zero face area there, so that face carries
no flux whatever its condition. On a line from the axis to $r=1$ cm in four
cylindrical cells the face areas are $0$, $0.0157$, $0.0314$, $0.0471$ and
$0.0628$ m² per metre of axial length and the volumes sum to $\pi r^2$ exactly
(`geometry!` on `create_uniform_mesh(0.01, 4; coordinates = Cylindrical())`).

The volume is **signed**: an inverted trial mesh gives a non-positive volume,
which the state mapping floors for evaluation (§13.2.4) and the admissibility test
reads unfloored (§15.5.4). Under a variable cross section $A$ is an algebraic
unknown of the state (§11.3.3) and enters every area and volume.

### 13.2.2 Cell centres and distances along the normal

The centre of a cell is the midpoint of its two nodes,
$z_i^c=\tfrac12(z_i+z_{i+1})$, in every coordinate system (a radial cell's centre
is the midpoint of its radii, not its volume centroid). Distances are measured
along the face normal:

$$
d_{f,i}=\pm\,(x_f-x_i^c)\cdot n_f\quad(\text{positive out of cell }i),\qquad
d_o+d_n=(x_n^c-x_o^c)\cdot n_f ,
$$

`center_to_face_distance(g, topology, f, i)` and `cell_distance(g, topology, f)`.
On an interior face of a line, $d_o=z_f-z_o^c$ and $d_n=z_n^c-z_f$ are the
half-widths of the two cells and $d_o+d_n$ is the centre-to-centre distance, the
denominator of every two-point gradient. On a uniform mesh $d_o=d_n=\Delta z/2$.
A gradient formed over the cell width instead of the centre distance would carry a
zeroth-order error wherever the mesh's non-uniformity does not vanish under
refinement.

### 13.2.3 Boundary distances

At a boundary face the neighbour is absent, and every boundary flux is formed over
the centre-to-face distance $d$ of the face's cell, its half-width
$\tfrac12\Delta z_i$ on a line: the surface energy balance
$q_\mathrm{bc}(T_s)=k\,(T_s-T_\mathrm{cell})/d$ (§12.2), the film-law mass
conditions (§12.4) and the pressure conditions (§12.5). On the face that bounds
the live cells after a deactivation, $d$ is the distance from the outermost live
cell's centre to that face.

The interior kernels (conduction, gas transport) visit only interior faces between
two live cells; the boundary faces, and the face that bounds the live cells after
a deactivation, belong to the boundary kernels (§13.7.5). No flux is therefore
computed twice, and no distance other than $d$ enters a boundary flux.

### 13.2.4 Distance and volume floors

Each kernel holds, per face, floors of the two centre-to-face distances:
`DISTANCE_FLOOR_FRACTION = 1e-8` times their values on the reference mesh. Every
distance a flux divides by passes through `smooth_floor(d, d_floor)`, and the
volume used for the concentrations through `smooth_floor(V, V_floor)` with
$V_\mathrm{floor}=10^{-8}\,V_i(0)$ (`VOLUME_FLOOR_FRACTION`). The floors act only
on inverted or collapsed trial meshes: a cell is deactivated long before its
volume approaches $10^{-8}$ of its initial value, and `smooth_floor` is the
identity above twice the floor (§13.10.3), so on every admissible state the
floored value is the geometric one.

### 13.2.5 Meshes of more than one dimension (residual level)

`MeshGeometry` also stores the face normals and centres, the **sweep vectors**
`face_sweep[q, f]` $=\int_f N_q\,\mathbf n\,dA$ (the volume face $f$ sweeps per unit
velocity of its $q$-th node; they sum to the area vector $A_f\mathbf n_f$) and the
lumped corner volumes `cell_node_volume[a, i]` $=\int_i N_a\,dV$. On a line the
sweep vector of a face is its area vector and the corner volumes integrate the two
linear shape functions over the shell by two-point Gauss–Legendre quadrature.

On planar meshes of quadrilaterals and hexahedra (dimension $D=2$ or 3; `box_mesh`,
`cell_mesh`, §13.9.6; 2D cells per unit depth) `geometry!` integrates the sweep vectors exactly over
straight or bilinear faces, takes each face's area vector as their sum, its area
as the vector's length and its normal as the vector over the area, computes the
volumes from the divergence identity
$V_i=\tfrac1D\sum_f s_{i,f}\sum_q\mathbf x_q\cdot\mathbf s_{q,f}$, and the
centroids and corner volumes by Gauss quadrature of the multilinear map;
`cell_volume_quadrature` is an independent path the tests compare against
(`test/multid/test_meshes.jl`). These meshes are assembled and differentiated at
the residual level (the reactions, conduction with conductivity tensors, the
thermal boundary, element death §14.5, the potential closure §11.3.5); they are
not a supported path for full solves. In-depth radiation, gas transport, the
variable cross section, deactivation and the equidistribution correction are a
line's: `setup` refuses them on a mesh of more dimensions.

---

## 13.3 The Divergence: Face Fluxes Scattered into Extensive Rows

### 13.3.1 The scatter

By the divergence theorem the rate of change of a cell inventory is its source
minus the net flux out through its faces. With $F_f$ the flux per unit area along
the normal of face $f$,

$$
-\oint_{\partial\Omega_i}\mathbf F\cdot\hat{\mathbf n}\,dS
\;\longrightarrow\;
-\sum_{f\in\partial i}s_{i,f}\,A_f\,F_f .
$$

The assembler loops over the faces once: with $(o,n)$ the face's cells it
subtracts $F_fA_f$ from the row of $o$ and, on an interior face, adds it to the
row of $n$. For cell $i$ of a line, face $i$ (below, cell $i$ its neighbour)
contributes $+A_iF_i$ and face $i+1$ (above, cell $i$ its owner) contributes
$-A_{i+1}F_{i+1}$: the inflow from below minus the outflow above. On the bottom face
cell 1 is the owner and the normal points out through the back, so a positive
$F_1$ removes $A_1F_1$ from cell 1.

The rows are **extensive** (kg s⁻¹ for masses, W for energy) and the code never
divides by a volume: $(1/V_i)\sum_f s_{i,f}A_fF_f$ is the discrete divergence of §2,
positive for net outflow, but the volume enters only through the concentrations
$\xi_{j,i}=m_{j,i}/\tilde V_i$ the kernels read.

> **Implementation note.** `assemble_rows!` (`src/Assembly/residual.jl`) is the
> only writer of `du`. It scatters `F_energy` into the energy rows and `F_mass`
> into the face-flux part of the mass rates (`M_flux`), writes the mass, energy,
> intensity, position, velocity, history, cross-section and deformation rows, and
> writes the dead rows last.

### 13.3.2 The mass rows

$$
\frac{dm_{j,i}}{dt}=\omega_{j,i}\tilde V_i+\dot M_{j,i},\qquad
\dot M_{j,i}=-\sum_{f\in\partial i}s_{i,f}\,A_f\,F_{j,f},
$$

with $\omega_{j,i}\tilde V_i$ [kg s⁻¹] the cell's reaction source (plus an external
mass source, §13.5.3). The face fluxes $F_{j,f}$ are those of the gas transport
(gas components only, §13.7) and, under a correction, the relative fluxes (every
component, §13.8). On a mesh that rides the skeleton, a condensed component has
no face flux and its row is its reaction source. Because the masses are the
states, a shrinking cell needs no dilation term: its concentrations rise as its
volume falls (§11.3.4).

### 13.3.3 Discrete conservation

Summing the face-flux parts over the cells, every interior face appears once with
$s_{i,f}=+1$ and once with $s_{i,f}=-1$ and the same $A_fF_{j,f}$, so

$$
\sum_{i}\dot M_{j,i}=-\sum_{f\ \text{boundary}}A_f\,F_{j,f},
$$

the discrete analogue of Gauss's theorem. The stoichiometry closes per reaction
(§6.5.3), so $\sum_j\omega_{j,i}=0$ and the column's total mass changes only by the
boundary outflow. The relation is linear in the state, so every step of a
Runge–Kutta or Rosenbrock method preserves it to round-off; the mass ledger of §16
closes to the integrator's tolerance on that basis. `test/assembly/test_residual_jacobian.jl`
("row assembly conserves mass and scatters fluxes antisymmetrically") asserts it on
a closed slab: the mass rows sum to zero within $10^{-12}$ of the largest mass row,
and in the enthalpy form the energy rows sum to the reaction heats and formation
terms within $10^{-12}$ relative, the conduction fluxes cancelling pairwise.

### 13.3.4 The energy row

With $S_{E,i}$ [W] the cell's energy sources (the reaction heat
$Q_{\mathrm{rxn},i}\tilde V_i$ and the absorbed in-depth radiation, §13.5) and $F_f$ the face energy flux (the
conduction, the enthalpy the mass fluxes carry, and on a boundary face the flux of
the condition, §13.7), the cell's energy right-hand side is

$$
E_i=S_{E,i}-\sum_{f\in\partial i}s_{i,f}\,A_f\,F_f .
$$

**Temperature form** (`energy_state = :temperature`, the default):

$$
\frac{dT_i}{dt}=\frac{E_i-\sum_j h_j(T_i)\,\dot M_{j,i}}{C_i}+S_{T,i},
\qquad C_i=\mathrm{smooth\_floor}\Bigl(\sum_j m_{j,i}\,c_{p,j}(\tilde T_i),\ \rho c_{p,\mathrm{floor}}V_i(0)\Bigr),
$$

with $\rho c_{p,\mathrm{floor}}=10^3$ J m⁻³ K⁻¹. Substituting $E_i$ and the face
flux of §13.7,

$$
C_i\frac{dT_i}{dt}=S_{E,i}-\sum_{f\in\partial i}s_{i,f}A_f\Bigl[q_f+\sum_jF_{j,f}\bigl(h_j(T_f)-h_j(T_i)\bigr)\Bigr]+C_iS_{T,i}:
$$

the gas enthalpy enters through its difference from the cell enthalpy. At an
incoming boundary, use the specified incoming temperature and reservoir heat;
an LTE film uses the equivalent converted wall-film representation (§9.9). A mass flux carried at the cell's own temperature neither heats
nor cools it; gas that leaves through a hotter face takes the enthalpy difference
with it. This is the semi-discrete form of the gas-advection term
$-\sum_j c_{p,j}\,\mathbf J_j\cdot\nabla T$ (§7.4), written as a face flux so that
it telescopes with the rest.

**Enthalpy form** (`energy_state = :enthalpy`):

$$
\frac{dH_i}{dt}=E_i+\sum_j\omega_{j,i}\tilde V_i\,h_j(T_i)+C_iS_{T,i},
\qquad H_i=\sum_j m_{j,i}h_j(T_i)\ \ (\text{regularized, §13.10.4}).
$$

The formation term $\sum_j\omega_{j,i}\tilde V_ih_j(T_i)$ appears because the
reaction heats are quoted at the local temperature (§6.6.3). The two forms are the
same balance, $dH_i/dt=C_i\,dT_i/dt+\sum_jh_j(T_i)\,dm_{j,i}/dt$. The enthalpy form is
a linear invariant every Runge–Kutta step conserves to round-off; the temperature
form conserves energy to the integrator's order (§16.5.2). $S_{T,i}$ is an external
temperature-rate source (§13.5.3), zero in a production solve.

### 13.3.5 Dead rows

After a deactivation (§14.3) the rows of a dead cell, and of the nodes and faces
beyond the live region, are identity rows: a zero derivative on a differential row
and $u_r$ itself on an algebraic row (`_dead_rows!`). The state length never
changes, the integrator's matrix stays nonsingular, and the live block is the
compacted problem's (§15.12). The interior kernels skip every face with a dead
cell, and the boundary kernels apply their conditions at `live_face` (§13.1.2).

---

## 13.4 Kernels and the Assembly Order

![Flow from fixed topology and reference data through setup and evolving state to geometry, derived fields, correction, kernels, closures and residual assembly.](../../assets/figures/02_mesh_framework.svg)

**Figure 2. Mesh framework (schematic).** Read left to right for one residual evaluation. Setup fixes connectivity and row layout; positions, inventories, cross section and the activity mask evolve. Geometry and derived fields feed correction, kernels and closure before assembly. Correction precedes face fluxes, and gas transport precedes the thermal boundary calculation. External sources and masking enter the final assembly; §13.4 gives the implementation details.

### 13.4.1 The kernel contract

A kernel (a subtype of `Kernel`, `src/Kernels`) implements
`contribute!(c, kernel, ds, topology, material, t)`: it adds its face fluxes and cell
sources to a `Contributions` buffer `c`, reading only the derived state `ds`
(§13.10), the topology and metrics, the material and the boundary conditions it
holds. It never reads the state vector, the row derivatives, or another kernel's
output except through the buffer, in the order the assembler fixes. The buffer's
fields:

| Field | Size | Units | Written by |
|---|---|---|---|
| `F_mass` | $N_C\times n_\text{faces}$ | kg m⁻² s⁻¹ | gas transport, relative fluxes |
| `F_energy` | $n_\text{faces}$ | W m⁻² | conduction, gas transport (carried enthalpy), thermal boundary, relative fluxes |
| `S_mass` | $N_C\times n$ | kg s⁻¹ | reactions, external sources |
| `S_energy` | $n$ | W | reactions, in-depth radiation |
| `θ` | $n$ | s⁻¹ | mesh closure (strain rates) |
| `M_flux` | $N_C\times n$ | kg s⁻¹ | assembler ($\dot M_{j,i}$) |
| `T_face` | $n_\text{faces}$ | K | thermal boundary (surface temperatures) |
| `r_intensity` | $n_\text{faces}$ | W m⁻² | in-depth radiation (algebraic residuals) |
| `r_velocity`, `r_multiplier` | $D\,n_\text{nodes}$, $n$ | — | mesh closure (algebraic residuals) |
| `F_carried`, `S_carried` | $2\times n_\text{faces}$, $2\times n$ | kg m⁻² s⁻¹, kg s⁻¹ | relative fluxes, mesh closure (conversion-history pair) |
| `gradient` | $D\times n$ | K m⁻¹ | conduction ($D\ge2$) |

`reset!` zeroes the fluxes, the sources and the algebraic residuals before the
kernels run; the assembler zeroes `M_flux`; `T_face` and `gradient` are overwritten.

### 13.4.2 The kernels and their order

`setup` builds the kernel tuple in assembler order (`_kernels`,
`src/Integration/setup.jl`):

1. `ReactionKernel`, when the material has reactions (§13.5.1).
2. `GasTransportKernel{M}`, for mobile gas, with the declared `LumpedTransfer`
   or `DustyGas`. Gas-bearing sealed materials still evaluate their EOS;
   `NoGasTransport` is reserved for materials without gas components.
3. `BeerLambertKernel`, under `BEER_LAMBERT` (§13.5.2).
4. `ConductionKernel`, always (§13.7.1).
5. `ThermalBoundaryKernel`, always (§13.7.5).
6. `RelativeFluxKernel`, when the correction is not `NoCorrection`: the
   equidistribution correction or a prescribed mesh velocity (§13.8).

One residual evaluation (§15.3) runs `derive!` (the state mapping, §13.10), then
`correction!` (the node corrections $c_n$, §11.7), `reset!`, the kernels in this
order, `closure_rows!` (the strain rates $\theta_i$ from the reaction mass sources
at the trial temperature, the history sources and the velocity-row residuals,
§11.3), `add_sources!` (§13.5.3), and `assemble_rows!` (§13.3). The order matters in
three places:

- **Coupled boundary exchange.** The thermal boundary kernel solves and stores
  boundary species and carried-energy fluxes together. Interior gas transport
  uses its own shared face state.
- **The relative fluxes after the thermal boundary.** The surface balance and its
  blowing models see the gas-transport fluxes alone.
- **The closure after the kernels and before the external sources.** The strain
  rates are the reactions' alone; a manufactured mass source does not move the mesh.

Every physics option is a type parameter of the definition (the kernel tuple's
type, the layout's flags, the closure), so `residual!` has one method and
`find_kernel` resolves a kernel at compile time.

### 13.4.3 Buffers and element types

The residual is a pure function of the state, the time and an immutable
`CoreDefinition` (the layout, the state setup, the kernels, the closure and the
external sources), writing into a `CoreWorkspace` of buffers: the derived state and
the contributions with `Float64` elements, and their twins for every dual-number
element type the colored Jacobian evaluates the residual at, built on first use
(`src/Assembly/definition.jl`; §15.2). The sizes are fixed for the whole solve; a
dead cell keeps its slots and the kernels skip it. After warm-up a residual and a
colored Jacobian allocate nothing (`test/jacobian/test_allocations.jl`, on 100-cell
columns in both energy forms).

---

## 13.5 Cell Sources: Reactions, In-Depth Radiation, External Sources

### 13.5.1 Reaction sources

`ReactionKernel` evaluates, cell by cell, at the concentrations
$\xi_{j,i}=m_{j,i}/\tilde V_i$ as stored (not their positive part) and at the trial
temperature $\tilde T_i$: the reaction rates $r_r$ (Arrhenius with the smooth
temperature gates and the material's `DepletionLimiter`, the one-sided limiter whose
rate is zero at and below zero reactant concentration, §6.1–§6.3), the species
sources $\omega_{j,i}=\sum_r\nu_{r,j}r_r$ and the reaction heat
$Q_{\mathrm{rxn},i}=-\sum_rh_rr_r$ with the heats quoted at the local temperature
(§6.6). It adds

$$
S_{\text{mass},j,i}\mathrel{+}=\omega_{j,i}\tilde V_i,\qquad
S_{E,i}\mathrel{+}=Q_{\mathrm{rxn},i}\tilde V_i .
$$

With the limiter enabled, a cell in which every reaction has a reactant at or
below zero concentration is skipped: its rates and every derivative of them are
exactly zero. The mass sources are also what the mesh closure reads for the strain
rates (§11.3.1).

### 13.5.2 In-depth radiation: the face intensities

Under `BEER_LAMBERT` the intensity of every face is an algebraic unknown of the
state. `BeerLambertKernel` orders the cells from the exposed face (the single
`:top` face) into the depth and, for each cell $i$, records its upstream face
(toward the exposed face) and its downstream face. Its rows are

$$
0=I_e-a\,q_\mathrm{inc}(t)\quad\text{(the face the beam enters)},\qquad
0=I_\text{down}-I_\text{up}\,e^{-\alpha_i\Delta_i}\quad\text{(the downstream face of each live cell)},
$$

with $q_\mathrm{inc}$ the incident flux of the `RadiativeFluxBC`, $a$ its explicit
absorptivity or else the effective absorptivity of the first live cell (§5.8.2),
$\alpha_i$ the cell's absorption coefficient at the positive-part concentrations
and $\Delta_i$ its thickness (the sum of its two floored centre-to-face distances).
The cell absorbs

$$
S_{E,i}\mathrel{+}=I_\text{up}A_\text{up}-I_\text{down}A_\text{down},
$$

the exact cell average of the attenuated beam, so the absorbed energy telescopes
to the entering intensity minus the transmitted one
(`test/kernels/test_radiation.jl`, "the beam's energy telescopes"). Each row reads
one cell and its two faces, so the Jacobian stays banded where the eliminated
recurrence would couple every cell to every cell above it. After cells are
deactivated from the exposed side the beam enters the first live cell's upstream
face, and the faces upstream of it carry identity rows. Under `BEER_LAMBERT` the
incident flux contributes nothing to the surface balance; the condition's other
members (re-radiation, convection) stay at the surface. `setup` refuses in-depth
radiation off planar coordinates, where the attenuated intensity and the face
areas would disagree. The model and its limits are §8.3.

### 13.5.3 External sources

`VolumetricSources(mass, energy)` (`src/Assembly/definition.jl`) adds external
sources for verification by manufactured solutions (§18.2): `mass(z, t)` returns
the volumetric mass source of every component [kg m⁻³ s⁻¹], added to `S_mass` as
$S\,\tilde V_i$ after the closure has read the reaction sources; `energy(z, t)`
returns a temperature-rate source $S_T$ [K s⁻¹], added to the temperature row as
$S_T$ and to the enthalpy row as $C_iS_T$. Both are pure functions of the cell-centre
position and the time, so the colored Jacobian differentiates the residual with
them exactly. The default is `NoSources()`. The ledger of §16 does not include them.

---

## 13.6 Gradient and Face-Value Reconstruction

### 13.6.1 Two-point normal gradients

On a line every gradient a flux uses is the two-point difference along the face
normal, for a cell field $\Phi$,

$$
\Bigl(\frac{\partial\Phi}{\partial n}\Bigr)_f=\frac{\Phi_n-\Phi_o}{d_o+d_n},
$$

second-order at the face for a smooth field on a uniform mesh and on a mesh whose
non-uniformity is smooth. No cell-centred gradient is formed on a line: neither
conduction nor gas transport needs one, and the enthalpy the gas carries is a face
flux (§13.7.4), not a cell source built on a cell gradient. At a boundary face the
gradient is one-sided over the centre-to-face distance, $(\Phi_s-\Phi_i)/d$,
with $\Phi_s$ the face value the condition supplies (the surface temperature of
the surface balance, the film law's face concentration); it is first-order at the
face (§13.1.3).

### 13.6.2 Face temperatures

Four face temperatures appear, each where it is consistent:

| Where | Face temperature |
|---|---|
| enthalpy carried by the gas fluxes, interior face | $T_f=(d_nT_o+d_oT_n)/(d_o+d_n)$, the linear interpolation to the face at the cell temperatures (exact for a linear profile) |
| inside the lumped wall flux | $\bar T_f=\tfrac12(\tilde T_o+\tilde T_n)$, the mean of the trial temperatures (§13.7.2) |
| boundary face | the surface temperature $T_s$ of the condition: closed form or surface Newton (§13.7.5) |
| relative fluxes | $\tilde T_f$, the upwind van Albada reconstruction on the stencil of the mass fluxes it carries (§11.6.2) |

### 13.6.3 Upwind face values

Two fluxes are advective. The Darcy flux carries the upwind cell's intrinsic gas
density $\xi/\varphi_f$, first-order (§13.7.3). The relative fluxes of a corrected
moving mesh carry a second-order upwind reconstruction with a smooth van Albada
slope and a positivity-preserving bound on every inventory (§11.6.2–§11.6.3); a
smooth slope is used because a limiter with a kink would put the kink into the
Jacobian and into finite-difference sensitivities. Both read the upwind side from
the value of the velocity alone, so at zero velocity every derivative pass takes
the derivative of the $v\ge0$ side, whatever its dual seeds
(`test/kernels/test_gas_transport.jl`, "at zero Darcy velocity every seed takes the
v ≥ 0 side").

### 13.6.4 Least-squares cell gradients (more than one dimension)

On a mesh of quadrilaterals or hexahedra, the conduction kernel computes for every
live cell the least-squares temperature gradient over its live face neighbours $k$,
each difference weighted by the inverse squared distance,

$$
\bigl(\mathsf N_i+\epsilon_\mathrm{LS}\,\mathrm{tr}(\mathsf N_i)\,I\bigr)\nabla T_i=\sum_k w_k\,(T_k-T_i)\,\mathbf d_k,
\qquad \mathsf N_i=\sum_k w_k\,\mathbf d_k\mathbf d_k^{\mathsf T},\quad w_k=|\mathbf d_k|^{-2},
$$

with $\mathbf d_k$ the vector between the cell centres and $\epsilon_\mathrm{LS}=10^{-8}$
(`GRADIENT_REGULARIZATION`) where the neighbours do not span every direction (the
normal matrix has $\det\mathsf N_i\le10^{-6}(\mathrm{tr}\,\mathsf N_i/D)^D$, `GRADIENT_SPAN`, decided
on its plain values; zero otherwise): a column one cell wide, or a cell beside dead ones, then has
a zero gradient component in the unresolved direction instead of an undetermined
one. The gradients feed the cross-diffusion part of the face flux (§13.7.1) and the
tangential correction of the thermal boundary kernel.
`test/multid/test_conduction.jl` asserts exact fluxes for linear fields with rotated
anisotropic tensors on distorted meshes and second-order convergence of a steady
manufactured solution with a rotated tensor on a distorted mesh.

---

## 13.7 Conductive, Gas-Diffusive, and gas Face Fluxes

This section states the discrete face fluxes. The physical derivations are in §7
(conduction, the gas enthalpy flux), §9 (gas transport and pressure) and §12
(boundary conditions).

### 13.7.1 Conductive heat flux

`ConductionKernel` visits the interior faces between two live cells:

$$
q_f=-k_f\,\frac{T_n-T_o}{d_o+d_n},\qquad
k_f=\frac{(d_o+d_n)\,k_ok_n}{d_ok_n+d_nk_o+\eta\,(d_o+d_n)},\quad\eta=10^{-20}\ \text{W m}^{-1}\text{K}^{-1}.
$$

$k_f$ (`guarded_harmonic_mean`, guard `FACE_MEAN_GUARD`) is the distance-weighted
**resistance-in-series** mean $(d_o+d_n)/(d_o/k_o+d_n/k_n)$ in product form: equal
to it to a relative $\eta/k_f$, it goes to zero continuously, with a finite
derivative, as either cell's conductivity vanishes. On a uniform mesh it is the
harmonic mean $2k_ok_n/(k_o+k_n)$. The cell conductivities $k_i$ are the effective
conductivities at the positive-part concentrations and the trial temperature (the
mixing rule of §5), smoothly floored at $10^{-6}$ W m⁻¹ K⁻¹ (`CONDUCTIVITY_FLOOR`);
the temperatures in the difference are the cell temperatures themselves. With
constant $k$ on a uniform mesh this is the textbook central difference; with a
conductivity jump on a face the steady profile is nodally exact (V-HC-04); with
$k(T)$ it carries a genuine second-order discretization error (V-HC-05).

On a mesh of more than one dimension the face's conductivity tensor is
$K_f=k_f$ times the distance-weighted mean of the cells' dimensionless tensors, and
the flux along the normal splits into a two-point part and a cross-diffusion part,

$$
q_f=-(\mathbf n\cdot K_f\mathbf n)\,\frac{T_n-T_o}{d_o+d_n}-\mathbf t\cdot\nabla T_f,\qquad
\mathbf t=K_f\mathbf n-(\mathbf n\cdot K_f\mathbf n)\,\frac{\mathbf d}{d_o+d_n},
$$

with $\nabla T_f$ the distance-weighted mean of the two cells' least-squares
gradients (§13.6.4), so a linear temperature field's flux is exact on any mesh whose
neighbours span every direction. On a line $\mathbf d$ is along the normal, the cross
part vanishes identically and is not evaluated.

### 13.7.2 Gas wall and molecular flux

The pore state $p_j=\xi_jRT/(M_j\phi)$ supplies the transported potential.
A lumped half path gives $n_j^D=\Lambda(p_{j,L}-p_{j,R})/(RT_{path}d)$.
Dusty-gas paths solve the coupled friction matrix using a consistent mean
face pressure and collective donor composition (§9.4–9.5). The conditioned
rank-one evaluator is linear in gas-species count and has explicit vacuum
and pure-gas limits.

---

### 13.7.3 Total-flux interface matching

Each half path adds Darcy flux $v_D p_{j,up}/(RT_{up})$. The interface partial
pressures satisfy equality of total species flux through the two halves.
This preserves a wall-only half connected to a Darcy-only half. The lumped
and pure-gas systems have a positive analytic quadratic solution; general
dusty-gas systems use a qualified local root and its first implicit derivative.
Pressure increments retain resolved gradients near equilibrium without
subtracting reconstructed absolute pressures. Chapter 9 gives degeneracies,
donor assumptions and root acceptance criteria.

---

### 13.7.4 Carried enthalpy

The single transmitted molar vector becomes $F_{j,f}=M_j n_{j,f}$ and carries
$\sum_jF_{j,f}h_j(T_f)$, where
$T_f=(d_RT_L+d_LT_R)/(d_L+d_R)$. Both species and energy rows scatter this
shared face value conservatively. Temperature and enthalpy formulations
use the same fluxes.

---

### 13.7.5 Boundary faces

`ThermalBoundaryKernel` assembles the complete boundary exchange. Named
prescribed fluxes, reservoir transport, external films and hydraulic
resistances provide the species closure. Outward species carry surface
enthalpy; inward species carry their specified incoming enthalpy.
For an LTE wall film, $q_r=q_s+H_r-H_s$ converts the wall heat to the
reservoir control surface and $E_{out}=H_r-q_r$. For adiabatic or prescribed
heat, donor enthalpy and imposed heat enter directly (§9.9).

The mass and thermal residuals are coupled only when their constitutive
closures require it. The solve reports its residual and iterations, and
implicit differentiation includes all local couplings. The boundary state
is a diagnostic, not a hidden solution-history input.

---

### 13.7.6 Solid and liquid components carry no transport flux

Condensed components are not transported by diffusion or Darcy flow: the
gas-transport kernel writes fluxes for the gas components only. On a mesh that rides
the skeleton, a condensed component's row is its reaction source, its mass stays in
its cell, and the mesh motion is what carries it (§11.1.2). The gas-transport fluxes
are relative to the skeleton; on a moving mesh the gas is mesh-bound, and these
fluxes are the whole of its transport (§11.2.3). Under a correction the relative
fluxes carry every component, condensed and gas, across the moving faces (§13.8).

---

## 13.8 Fluxes Relative to a Moving Mesh (Recap)

On a moving mesh the node velocities are $w=v_s+c$, the solid velocity from the
balance rows plus a correction (§11.3, §11.7). Where $c\ne0$, `RelativeFluxKernel`
carries every inventory across the faces at the relative velocity
$v_{\mathrm{rel},f}=-c_f\,n_f$:

$$
F_{j,f}\mathrel{+}=v_{\mathrm{rel},f}\,\xi_{j,f},\qquad
F_f\mathrel{+}=\sum_jF_{j,f}\,h_j(\tilde T_f),
$$

and the conversion-history pair at the solid mass flux times its per-solid-mass
values. The face states are second-order upwind reconstructions with the smooth van
Albada slope, an inventory's face state is bounded by the upwind cell's inventory,
a boundary face takes the one-sided extrapolation of its cell, and a face with an
exactly zero relative velocity is skipped (§11.6). The energy is never reconstructed
separately, so an isothermal composition profile advected across the mesh stays
isothermal. A prescribed mesh velocity (`w_override`, §11.7.8) moves the mesh over a
solid at rest, and the same kernel carries every inventory at the mesh-relative
velocity. Under `NoCorrection` the kernel is absent.

No dilation term is assembled: the masses are the states, so a cell's concentrations
follow from its masses and its volume (§11.3.4), and the volume is recomputed from
the node positions at every evaluation, so the discrete GCL holds identically
(§11.5).

---

## 13.9 Mesh Generation

### 13.9.1 Line meshes

`line_mesh(x; coordinates = Planar(), cross_section, bottom = :bottom, top = :top,
fixed = false, lagrangian = Int[])` (`src/Meshes/generators.jl`) builds a line with
nodes at the strictly increasing, finite positions `x` (non-negative radii in
cylindrical and spherical coordinates), cell $i$ between nodes $i$ and $i+1$, face
$f$ at node $f$, and the topology of §13.1. It pins node 1 (every node with
`fixed = true`) and records the Lagrangian nodes. The cross section must be finite
and positive and defaults to 1, 1 and $4\pi$ in planar, cylindrical and spherical
coordinates. The three user constructors below build on it; each returns a
`Mesh{1}`.

### 13.9.2 Uniform mesh

`create_uniform_mesh(L, n; coordinates = Planar(), cross_section)` divides $[0,L]$
into $n\ge1$ equal cells:

$$
\Delta z=\frac Ln,\qquad z_k=(k-1)\,\Delta z\quad(k=1,\dots,n+1),
$$

with the volumes of §13.2.1 ($V_i=A\,\Delta z$ in planar coordinates). In
cylindrical and spherical coordinates the positions are radii, so the line runs from
the axis or the centre to the surface at radius $L$.

### 13.9.3 Geometric-stretch mesh

`create_stretched_mesh(L, n, s; stretch_from = :top, coordinates, cross_section)`
concentrates cells near one boundary to resolve the steep thermal and reaction fronts
there. Neighbouring cells differ by the ratio $s\ge1$, and the cell widths form a
geometric series summing to $L$:

$$
L=\Delta z_0\,\frac{s^n-1}{s-1}\quad\Longrightarrow\quad\Delta z_0=L\,\frac{s-1}{s^n-1}.
$$

The series is built from the finest cell at $z=0$, $z_k=z_{k-1}+\Delta z_0\,s^{k-2}$
($k=2,\dots,n+1$), and scaled by $L/z_{n+1}$ to end at $L$. `stretch_from` names the
boundary that receives the finest cell: with the default `:top` the positions are
reflected, $\tilde z_k=L-z_{n+2-k}$, so the finest cell lies at the exposed surface;
`:bottom` keeps the finest cell at the back. The end positions are then set to $0$
and $L$ exactly. For $|s-1|<10^{-10}$ the mesh is the uniform range. The constructor
refuses $s<1$ and any `stretch_from` other than `:top` and `:bottom`: the
orientation is chosen by `stretch_from` alone, never by a ratio below one.

### 13.9.4 Mesh from a node list

`create_mesh_from_nodes(z; coordinates, cross_section)` builds a line on strictly
increasing positions that start at zero. It puts a node on every layer interface so
that no cell straddles a composition change (§13.1.3), and builds graded layered
meshes such as the one of V-BO-06 (§18).

### 13.9.5 Example

```julia
using Pyrolysis
using Pyrolysis.Internal: geometry!, center_to_face_distance, boundary_faces

mesh = create_stretched_mesh(0.01, 5, 1.2)      # 1 cm, 5 cells, finest at z = L
topo = mesh.topology
x = [p[1] for p in mesh.positions]
println(diff(x))                                 # cell widths [m]
println(topo.faces[1], "\n", topo.faces[2], "\n", topo.faces[6])
println(boundary_faces(topo, :top), " ", topo.pinned)

g = MeshGeometry{Float64}(topo)
geometry!(g, x, topo.cross_section, topo)
println([n[1] for n in g.face_normal])           # face normals
println(center_to_face_distance(g, topo, 2, 1), " ",
        center_to_face_distance(g, topo, 2, 2))  # d_o, d_n of face 2
```

prints (the type prefix `Pyrolysis.Meshes.` shortened)

```
[0.00278649752741346, 0.002322081272844549, 0.0019350677273704585, 0.0016125564394753808, 0.0013437970328961518]
Face{1, 1}((1,), (1, 0), :bottom)
Face{1, 1}((2,), (1, 2), :interior)
Face{1, 1}((6,), (5, 0), :top)
[6] [1]
[-1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
0.00139324876370673 0.0011610406364222746
```

The widths fall by the factor 1.2 toward the surface, the finest equals
$L(s-1)/(s^n-1)=1.3438$ mm, the bottom face's normal points out through the back,
and the two distances of face 2 are the half-widths of cells 1 and 2.

### 13.9.6 Quadrilaterals and hexahedra (residual level)

`box_mesh(xs, ys[, zs]; transform = identity, pinned = (:bottom,), fixed = false,
lagrangian = Int[])` builds a structured planar mesh of quadrilaterals or hexahedra
on the tensor product of the node coordinates, each node then moved by `transform`;
its sides are tagged `:left` and `:right` along $x$, `:front` and `:back` along $y$
(3D), and `:bottom` and `:top` along the last axis. `cell_mesh(positions, cells; tag,
…)` builds one from node positions and per-cell node lists, matching faces by their
node sets. These meshes serve the residual-level machinery of §13.2.5.

---

## 13.10 The State and its Mapping

### 13.10.1 The state vector

The state is a sequence of node and cell blocks followed by global scalars, laid out
by `StateLayout` (`src/States/layout.jl`; the full layout and its index formulas are
§15.1):

- **cell block** $i$: the component masses $m_{1,i},\dots,m_{N_C,i}$ [kg], the energy
  state ($T_i$ [K], or the sensible-enthalpy inventory $H_i$ [J] with
  `energy_state = :enthalpy`), and under a variable cross section the
  conversion-history pair $(G_i, m_{\mathrm{dry},i})$;
- **node block** $n$: on a moving mesh the position $z_n$ (differential) and the
  velocity $w_n$ (algebraic); under `BEER_LAMBERT` the intensity of the face at
  node $n$ (algebraic);
- **globals**, under a variable cross section: the area $A$ and its lateral strain
  rate $\theta_A$ (algebraic).

On a line the order is node 1, cell 1, node 2, …, cell $n$, node $n+1$, so every row's
unknowns and its neighbours' lie within a fixed distance and the Jacobian is banded.
No kernel computes an index: every access goes through `mass_index`,
`energy_index`, `position_index`, `velocity_index`, `intensity_index`,
`history_index`, `area_index` and `lateral_rate_index`. `setup` fills the initial
state from the problem: $m_{j,i}=\xi_j(z_i^c)\,V_i(0)$ and the energy state from
$T(z_i^c)$ at the cell centres of the reference mesh, then makes the algebraic
unknowns consistent (§15.4).

### 13.10.2 The state mapping

`derive!(ds, u, layout, setup)` (`src/States/derived.jl`) maps the state vector to
the `DerivedState` every kernel reads; it is the only reader of $\boldsymbol u$. In
order:

1. the live cells, from the layout's mask;
2. on a moving mesh, the node positions and velocities (or, under the potential
   closure of more dimensions, the cells' multipliers), under a variable cross section
   $A$ and $\theta_A$, and the
   geometry from the positions (`geometry!`, §13.2); a fixed mesh keeps the
   reference geometry computed when its buffers are built;
3. the carried fields (the history pair; in more dimensions the deformation
   gradient) and the face intensities;
4. per cell, the signed volume $V_i$ and $\tilde V_i=\mathrm{smooth\_floor}(V_i,10^{-8}V_i(0))$,
   the masses, the concentrations $\xi_{j,i}=m_{j,i}/\tilde V_i$ and the positive-part
   concentrations $\xi^+_{j,i}=(m_{j,i})^+/\tilde V_i$ at the cell's scale $\delta_i$;
5. the temperatures: the stored temperature, or the inversion of the enthalpy
   inventory (§13.10.4);
6. per cell, the trial temperature $\tilde T_i$, the floored heat-capacity sum $C_i$,
   the sensible enthalpies $h_j(T_i)$, the floored conductivity $k_i$, the surface
   emissivity $\varepsilon_i$ and absorptivity $a_i$, the pore-gas properties the gas
   transport needs (the wall coefficients, molecular resistance, porosity,
   partial pressures and mobility), and under in-depth
   radiation the absorption coefficient $\alpha_i$.

The node corrections $c_n$ are written next, by the closure (`correction!`, §11.7).
Which value of each field a kernel reads is fixed:

| Field | Read by |
|---|---|
| $\xi_{j,i}=m_{j,i}/\tilde V_i$ (signed, as stored) | reactions and relative solid-motion fluxes |
| $\xi^+_{j,i}$ | properties, pore pressure, species transport and reservoir exchange |
| $\tilde T_i$ | reactions, properties, the lumped wall flux, pressure, the surface balance, the strain rate |
| $T_i$ | the conductive difference, the sensible enthalpies, the face temperatures $T_f$ |

### 13.10.3 Trial extensions and physical admissibility

The shared C¹ guards are defined in `src/Materials/guards.jl`:
`smooth_clamp` for trial temperature, `positive_part` for inventories, and
`smooth_floor` for positive divisors. The porosity divisor floor is $10^{-12}$
and becomes the identity above twice that value. Such extensions let
intermediate integrator states be evaluated within their finite envelope;
they do not make a closed pore or an invalid constitutive root admissible.

Condensed caloric scales retain their reference-mass definition. The gas
transport positive-part scale is recomputed from current temperature and
pore volume, sharing the pressure-error budget among gas species. Its total
pressure bias is below 0.26 times that budget (§9.10). Accepted physical states and post-event states are checked
using raw quantities and a mass/thermal/geometry uncertainty margin. Local
transport failures are diagnosed explicitly (Chapter 9).

---

### 13.10.4 The energy state and the caloric extension

The sensible enthalpies $h_j(T)=\int_{T_\mathrm{ref}}^Tc_{p,j}\,dT'$ have closed forms
for constant, linear, polynomial and tabulated heat-capacity laws (each law constant
outside the property bracket, so the enthalpy continues linearly there), and a
`TabulatedEnthalpy` for any other law: the integral at 1 K knots on the bracket by
five-point Gauss–Legendre quadrature and the cubic Hermite interpolant between.

In the temperature form the cell temperature is read from the state. In the
enthalpy form it is the root of the cell's **extended caloric relation**
(`CaloricExtension`, `src/States/caloric.jl`)

$$
H_i=\sum_j(m_{j,i})^+\,h_j(T)+C_{\mathrm{reg},i}(\boldsymbol m)\,(T-T_\mathrm{ref}),
$$

where $C_{\mathrm{reg},i}$ [J K⁻¹] is a $C^1$ regularizing capacity that is zero while the guaranteed heat capacity
$\sum_j(m_{j,i})^+c_{p,j,\min}$ is at least twice the cell's capacity floor
$C_{\mathrm{floor},i}=\rho c_{p,\mathrm{floor}}V_i(0)$ and rises to $C_{\mathrm{floor},i}$ as
it vanishes ($c_{p,j,\min}$ is half the sampled minimum of the component's law). The
relation is strictly increasing in $T$ with slope at least $C_{\mathrm{floor},i}$, so it
has exactly one root for every real $H_i$, and it agrees with the physical relation
on every state in which each component is absent or holds at least $\delta_i$ and the
guaranteed capacity is at least $2C_{\mathrm{floor},i}$. `invert_caloric` solves it by a
bracketed Newton iteration started at the cell's temperature from the latest plain
evaluation, replacing a step that does not contract the bracket by bisection. A
dual-number inventory is differentiated through the defining relation: the primal
root is corrected by one Newton step per nesting level, which carries the
implicit-function derivatives exactly. A non-finite inventory has temperature `NaN`,
so the residual of such a trial state is `NaN` and the integrator rejects the step.

### 13.10.5 Reference scales and pressure-aware error control

`CellReference` holds initial volume, condensed mass, geometry, capacity
floors and the pressure-error budget. Condensed-mass and geometry absolute
tolerances retain their reference scales. Gas-row reference tolerances are
$\phi_i(0)V_i(0)M_jp_{tol}/(N_gRT_i(0))$.

The internal error norm also reconstructs each unscaled local error at the
proposed state and differentiates $p_j=m_jRT/(M_jW)$, with
$W=V-\sum_{condensed}m_j/\rho_{i,j}(T)$. Temperature, condensed masses,
positions, cross section and caloric inversion therefore contribute to the
pressure-error estimate. Each species shares the absolute budget, and the
sum of absolute partial-pressure errors is controlled by the total budget.
This is a local estimate; accumulated global error is qualified against
tighter trajectories. No residual evaluation mutates solver tolerances.

---

### 13.10.6 Example

The mapping of a set-up problem, read through `Pyrolysis.Internal`:

```julia
using Pyrolysis
using Pyrolysis.Internal: residual!, buffers, mass_index, energy_index, position_index,
    velocity_index, state_length, is_differential

solid = SolidComponent(:solid, ρ = 1190.0, φ=.3, c = 1420.0, k = 0.19, λ = 1e-6)
gas = GaseousComponent(:gas, M = 0.1, c = 1500.0, k = 0.02,)
material = Material(
    gas_transport=LumpedTransfer(endpoint=LumpedEndpoint(λ=1e-6,κ=0.,max_skeletal_occupancy=.001)),

    name = :slab,
    components = (solid, gas),
    reactions = (Reaction(:gasification, 1 => 2, A = 8.5e12, E = 188e3, h = 8.7e5),),
)
problem = PyrolysisProblem(
    mesh = create_uniform_mesh(0.002, 20),
    material = material,
    bc_set = BoundaryConditionSet(
        thermal = (
            top = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
        mass = (top = AmbientGasBC(GasReservoir(P=0., T=300.)), bottom = ImpermeableBC()),
    ),
    T_initial = z -> 300.0 + 2e5 * z,      # 300 K at the back, 700 K at the surface
    ξ_initial = [z -> 1190.0, z -> 0.0],
    tspan = (0.0, 60.0),
)
cp = setup(problem; use_ale = true)
layout = cp.def.layout
println(state_length(layout), " rows, ", count(!, is_differential(layout)), " algebraic")
println((position_index(layout, 1), velocity_index(layout, 1), mass_index(layout, 1, 1),
         mass_index(layout, 1, 2), energy_index(layout, 1), position_index(layout, 2)))

du = similar(cp.u0)
residual!(du, cp.u0, cp.def, 0.0, cp.ws)
ds, c = buffers(cp.ws, cp.def, Float64)
println((ds.m[1, 20], ds.V[20], ds.ξ[1, 20], ds.T[20], ds.cp_sum[20]))
println((c.F_energy[21], c.T_face[21]))
println((du[energy_index(layout, 20)], du[mass_index(layout, 20, 1)]))
```

prints

```
102 rows, 21 algebraic
(1, 2, 3, 4, 5, 6)
(0.11900000000000005, 0.00010000000000000005, 1190.0, 690.0, 168.98000000000008)
(-38157.35850998746, 700.0414101347891)
(-29.609607022784147, -0.0059319424191924845)
```

The state has $21\times2+20\times3=102$ rows: a position and a velocity per node and
two masses and a temperature per cell, the 21 velocities algebraic; node 1's block
comes first, then cell 1's, then node 2's. The top cell (centre at 1.95 mm, so
690 K) holds $1190\times10^{-4}$ kg in $10^{-4}$ m³ and $C=0.119\times1420$ J K⁻¹. The
surface balance of the `:top` face gives $T_s=700.04$ K, and its outward flux is
$-38.16$ kW m⁻² (into the domain: 50 kW m⁻² less the re-radiation
$0.9\,\sigma(T_s^4-300^4)$). The conduction into the cell and out of it nearly
cancel on this linear profile; the endothermic gasification at 690 K (a rate of
$0.05$ s⁻¹ on 1190 kg m⁻³) cools the cell at about 30 K s⁻¹ and removes its solid
at $59.3$ kg m⁻³ s⁻¹ times its volume.

---

## 13.11 Comparison to Gpyro and ThermaKin2Ds Discretizations

### 13.11.1 Gpyro

The documented 1D Gpyro method uses fully implicit control-volume balances,
with changing cell thicknesses, harmonic interface transport coefficients,
source linearization and tridiagonal solves in an iterative coupled solution
(GP §§4.1, 6). Boundary nodes lie on half-cell boundaries.

GP's convection–diffusion scheme is **selected by an input option**; central,
upwind, hybrid, power-law and exponential forms have different Péclet-number
coefficient functions (GP Table 2; 0.8200 `GET_A`/`GET_A_SINGLE`). It does
not automatically choose among all schemes as the Péclet number changes.
Pyrolysis.jl uses first-order upwind Darcy advection, a separate diffusive
law and distance-weighted harmonic face coefficients. Its semi-discrete
system is advanced with a stiff integrator and a coupled colored Jacobian.
Those shared finite-volume features do not make the diffusion laws or time
errors identical (§9.10).

### 13.11.2 ThermaKin2Ds

TK §3 calls its spatial elements **finite volumes**, carrying component
masses and temperature. It uses Crank–Nicolson integration in the primary
coordinate and explicit updates in the transverse coordinate, with changing
element sizes and conservative size adjustment. It should not be contrasted
with Pyrolysis.jl as a fixed-grid finite-difference method with no volume
inventory conservation.

Pyrolysis.jl carries extensive mass states, optional enthalpy states, and
the moving geometry in one mass-matrix system. Its source code establishes
these choices; TK's manual establishes its own method. Neither description
alone establishes relative accuracy, computational cost, or stability on a
particular ignition transient.

---

## 13.12 Limitations and Open Issues

1. **First-order gas advection.** The pressure-driven advective face flux is
   first-order upwind (§13.7.3); the second-order van Albada reconstruction exists
   only for the relative fluxes of a corrected moving mesh (§11.6). Sharp
   pressure-driven composition fronts are smeared at first order unless the mesh is
   refined.

2. **First-order boundary faces.** A boundary face is reached from its cell centre
   over the half-width (§13.2.3), so a face value fixed by a condition (a prescribed
   temperature, the solved surface temperature) carries a first-order error at the
   face. It is visible where a steep profile sits inside the boundary half-cell: in
   V-CPL-02 the surface-temperature closure, with the reaction zone inside that
   half-cell, shifts the steady-wave identity by about 1.3 % at $\Delta z\approx2\delta_r$
   and $7\times10^{-4}$ at $\Delta z\approx\delta_r$ ($\delta_r$ the reaction-zone
   thickness; the case's comments). The manufactured-solution cases impose boundary
   fluxes, not face values, for that reason (§18.2).

3. **No mixed cells.** Initial compositions are sampled at the cell centres, so an
   interface that does not lie on a face is moved to the nearest face. V-HC-04c
   measures the effect on a two-layer slab: the numerics match the moved-interface
   resistance network to $10^{-11}$, and the whole error against the true interface,
   1.64 % in flux and 1.47 K in the profile at $\Delta z=0.5$ mm, is the quantization.
   Put a node on every interface (§13.9.4).

4. **In-depth radiation on planar lines only.** `setup` refuses `BEER_LAMBERT` off
   planar coordinates and on meshes of more dimensions. Only the first
   `RadiativeFluxBC` member of a combined `:top` condition feeds the beam; a second
   radiative member is dropped without a warning (§8).

5. **More dimensions at the residual level.** Meshes of quadrilaterals and
   hexahedra are assembled and differentiated (§13.2.5, §13.6.4) but are not a
   supported path for full solves: gas transport, in-depth radiation, the variable
   cross section, deactivation and the equidistribution correction are a line's, and
   the relative fluxes are first-order there (§11.6.5).

For the temporal integration and the Jacobian that consume the rows of this chapter,
see §15 (Temporal Integration, the DAE, and the Colored Jacobian); for deactivation
and burnout, §14 (Cell Deactivation and Burnout); for the ledgers that exploit the
antisymmetric scatter of §13.3.3, §16 (Conservation Ledgers and the Audit).
