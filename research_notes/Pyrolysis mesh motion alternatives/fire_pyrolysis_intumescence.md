# Mesh and coordinate strategies in fire pyrolysis and intumescence models (shrinkage, strong swelling, burnout)

Scope note: this file covers the fire-science condensed-phase literature and codes. It focuses on how each discretizes a solid that changes thickness. Local copies of code manuals are cited as `file://` links: FDS 6.11.0 guides, the ThermaKin2Ds v2 reference, the Gpyro technical reference and user guide, and the Abdelhafiz et al. 2026 Gpyro V2 paper. Where I read source code directly (FDS 6.11.0 `wall.f90`/`read.f90`, OpenFOAM-6/7 `reactingOneDim.C`/`regionModel1D.C`), the GitHub link is given. Some items rest only on an abstract or a search-engine summary; those are flagged.

---

## Q1. How do FDS (incl. HT3D), ThermaKin2Ds, Gpyro (incl. V2/3D), Gpyro-derived codes, FireFOAM/OpenFOAM, Rein-group models, Di Blasi and Lautenberger handle shrinkage, swelling and burnout? What exact remeshing rules do they use?

### Takeaway
Every mainstream fire code uses a *mass-carrying, solid-Lagrangian 1D cell*: solid mass never crosses a cell face, and cell thickness follows from mass divided by density. The codes differ only in what they do when cells become too large or too small:
- **FDS** does a threshold-triggered *conservative remap* (renoding) and removes thin layers.
- **ThermaKin** *rezones every time step* to a preset Δx, splitting and topping up elements.
- **Gpyro and OpenFOAM** never rezone. Cells keep a fixed count, and below a minimum thickness they are frozen inert.

None of them integrates node positions inside an adaptive, error-controlled implicit solver. None extends shrink/swell cleanly to 3D: FDS HT3D forbids significant shrink/swell and burn-away, and Gpyro 2D/3D had no volume change as of 2017.

### Cited Findings

**FDS 6.11.0 (1D solid, default)**
- **Thickness update rule.** Each time step, every condensed-phase cell is multiplied by a factor δ built from ratios ρ_s,α/ρ_α. The documented rule is: δ = max_α(ρ_s,α/ρ_α) if that max is ≥1, otherwise δ = Σ_α(ρ_s,α/ρ_α). Densities are divided by δ to conserve mass — [FDS Tech Ref §7.2.5, Eq. 7.50](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/FDS/fds_technical_reference_guide.md); [FDS User Guide §9.3.4, Eq. 9.6](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/FDS/fds_user_guide.md).
  - The 6.11.0 source computes `REGRID_FACTOR` as the *sum* Σ_n ρ_n/ρ_adjusted,n. It clamps the factor to ≥1 if any present material has `ALLOW_SHRINKING=F`, and to ≤1 if any has `ALLOW_SWELLING=F`, then divides the densities by it. I did not find the max branch in the excerpt I inspected, so there may be doc/code drift — [FDS wall.f90 (FDS-6.11.0), ~L2380–2427](https://github.com/firemodels/fds/blob/FDS-6.11.0/Source/wall.f90).
- **Which face moves.** New node coordinates are rebuilt from the back toward the front: R_new(i) = (R_new(i+1)^I + ΔV_i·REGRID_FACTOR)^(1/I). With x=0 kept at the exposed face, the exposed face stays fixed and the solid "shrinks from the bottom up" — [FDS wall.f90 L2432–2440](https://github.com/firemodels/fds/blob/FDS-6.11.0/Source/wall.f90); [McDermott et al. 2018, FEMTC](https://files.thunderheadeng.com/femtc/2018_d2-07-mcdermott-paper.pdf).
- **Remesh triggers** (6.11.0 source). A layer is flagged for remeshing when any of these occur:
  1. A cell's new thickness drops below 20·ε.
  2. |Δx_new/Δx_at_last_remesh − 1| > `REMESH_RATIO`.
  3. For a shrinking cell, linear blending with a neighbour would change T by more than `RENODE_DELTA_T`.
  4. The layer's thickness has dropped by more than 1.5× its smallest cell since the last remesh.

  Source: [FDS wall.f90 L2443–2537](https://github.com/firemodels/fds/blob/FDS-6.11.0/Source/wall.f90).
- **`RENODE_DELTA_T`** is computed per material as `REAC_RATE_DELTA` (default 0.15) divided by d ln r/dT. For Arrhenius kinetics that is 0.15/(E/(R·T_ref²)), i.e. the temperature change that would alter the reaction rate by 15% — [FDS read.f90 L7829–7853](https://github.com/firemodels/fds/blob/FDS-6.11.0/Source/read.f90).
- **The user guide describes the same safeguards.** "If nodes change size during pyrolysis, FDS will try and renode layers… This renoding process is costly if it is done frequently." Renoding "will combine cells causing changes in temperature… this can cause large swings in the pyrolysis rate," which is why `REAC_RATE_DELTA` exists — [FDS User Guide §8.3.8](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/FDS/fds_user_guide.md).
- **`REMESH_RATIO` default is inconsistent across the docs.**
  - User-guide text: 0.25 (first renode check of a 0.1 m node happens at 0.125 or 0.075 m).
  - User-guide parameter table: 0.15.
  - 6.11.0 source: `REMESH_RATIO = 0.25`.

  Sources: [FDS User Guide §8.3.8 and SURF table](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/FDS/fds_user_guide.md); [FDS read.f90 L9136](https://github.com/firemodels/fds/blob/FDS-6.11.0/Source/read.f90).
- **Remap procedure** (6.11.0 source):
  1. Regenerate the layer's nodes with the original gridding rules (`GET_WALL_NODE_COORDINATES`: smallest cell size, stretch factor).
  2. Interpolate volume-integrated quantities (ρ·V per material, sensible enthalpy ρh·V, T·V, Q·V) with overlap weights (`GET_INTERPOLATION_WEIGHTS` / `INTERPOLATE_WALL_ARRAY`).
  3. Recover temperature from enthalpy by Newton iteration (tolerance 1e-4 K, ≤20 iterations).
  4. Reset the ghost temperatures so the front and back surface temperatures are unchanged.

  This remap conserves mass and sensible enthalpy — [FDS wall.f90 L2701–2803](https://github.com/firemodels/fds/blob/FDS-6.11.0/Source/wall.f90).
- **Default solid grid.** First cell ≤ √(τk/ρc) with τ = 1 s. `STRETCH_FACTOR` defaults to 2 (cells double toward mid-depth, then shrink). `N_LAYER_CELLS_MAX` defaults to 1000. The user guide notes pyrolysis can cause "spurious fluctuations in temperature and burning rate" and recommends uniform spacing (`STRETCH_FACTOR=1`) and smaller cells — [FDS User Guide §8.3.8](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/FDS/fds_user_guide.md); [FDS Tech Ref App. M, Eq. M.3](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/FDS/fds_technical_reference_guide.md).
- **Time integration.**
  - Crank–Nicolson in the solid ([FDS Tech Ref Eq. M.4](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/FDS/fds_technical_reference_guide.md)).
  - Updated every 2 gas steps by default (`WALL_INCREMENT`).
  - Sub-stepped up to `TIME_STEP_FACTOR` = 10 when interior T would change by more than `DELTA_TMP_MAX` = 25 °C ([FDS User Guide §8.3.8](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/FDS/fds_user_guide.md)).
- **Burnout handling.**
  - A layer is removed when its thickness falls below `MINIMUM_LAYER_THICKNESS` (default the lesser of 1e-6 m and 0.1×`THICKNESS`) or its mass fraction falls below `MINIMUM_LAYER_MASS_FRACTION` (default 1e-12). The source adds the remaining layer mass to the released mass flux (`M_DOT_LAYER_PP`).
  - Delamination can also drop a layer.
  - At zero total thickness without `BURN_AWAY`, surface T is set to `TMP_GAS_BACK`, convective flux to 0 and burning rate to 0.

  Sources: [FDS User Guide §8.3.8](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/FDS/fds_user_guide.md); [FDS wall.f90 L2501–2536, 2679–2697](https://github.com/firemodels/fds/blob/FDS-6.11.0/Source/wall.f90).
- **`BURN_AWAY=T`** removes gas-phase obstruction cells as their mass is consumed. Heat conduction in such an OBST is limited to 1D — [FDS User Guide §9.3.9](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/FDS/fds_user_guide.md).

**FDS HT3D (3D heat transfer)**
- **Current 6.11.0 restrictions.** HT3D cannot be applied to an OBST that is to `BURN_AWAY`. The guide also says "if the solid undergoes significant shrinking or swelling, do not use HT3D. 'Significant' means that the number of internal cells changes, in which case the 3-D nodal structure breaks down." The 3D solver reuses the 1D solver in an ADI scheme — [FDS User Guide §8.4.2](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/FDS/fds_user_guide.md); [FDS Tech Ref App. M](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/FDS/fds_technical_reference_guide.md).
- **Earlier 3D pyrolysis prototype** (McDermott, Zhang, Bruns, Benkorichi, FEMTC 2018):
  - **Grid.** Solid mass is tied to fixed Eulerian cells: "The grid cell itself cannot disappear." The method tracks a solid volume fraction φ_s = V_s/V_cell = Σρ_s,α/ρ_α.
  - **Deformation.** Conduction distances are rescaled: isotropic δx̃ = φ_s^(1/2)δx (2D), or unidirectional δx̃ = φ_s δx (e.g. cone test).
  - **Burn-away.** Cells burn away when their mass drops below a threshold, after the mass is shifted to a neighbour.
  - **Authors' caveat.** "the numerical schemes needed for shrinking and swelling materials in multi-dimensions are not straight-forward, and the choices made in how to deal with intumescence and burn away can have leading-order effects on the model results."

  Source: [McDermott et al. 2018](https://files.thunderheadeng.com/femtc/2018_d2-07-mcdermott-paper.pdf). The 2018 prototype's burn-away capability conflicts with the 6.11.0 user-guide prohibition, so it was apparently not retained in current HT3D. That reading is mine, from comparing the two documents.

**ThermaKin2Ds (Stoliarov group)**
- **Element descriptors.** Elements are described by component masses and temperature, which are the "primary object descriptors."
- **Deformation in the equations.** The conservation equations carry a deformation-advection term ∂/∂x(ξ_j ∫₀ˣ (1/ρ)(∂ρ/∂t)dx). This is a material velocity integrated from the back boundary (x = 0, "the reference point for the expansion or contraction").
- **Swelling factor γ** (0–1) sets how much generated gas adds volume through the mixture-density rule.
- **Integration.** Crank–Nicolson in x with fixed user Δt; explicit in y/z.
- **Rezoning.** "element sizes are adjusted before every timestep. If an element is larger than a pre-set value of Δx, it is split in two. If it is smaller, a fraction of the following element is added to bring it to the pre-set value. The temperature and composition of the mixed element are recalculated to ensure the conservation of energy and species."
- **Recommended settings.** Δx = 5×10⁻⁵ m and Δt = 0.01 s are suggested as a starting point, with sensitivity checks.
- **Radiation absorption** uses either a "maximum absorption" or a Monte-Carlo "random absorption" element-selection algorithm at every time step.

  Source for all of the above: [ThermaKin2Ds Reference v2, §1, §3](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/ThermaKin/ThermaKin2Ds%20Reference%20v2.md).
- **Deformation is along x only.** x carries radiative transport and object deformation. The front boundary may become non-uniform in 2D. For 2Dax objects the non-uniform front profile's effect on surface area and exposure is accounted for — [ThermaKin2Ds Reference v2, §2–3](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/ThermaKin/ThermaKin2Ds%20Reference%20v2.md).
- **Practice in Swann's thesis** (UMD 2019, Stoliarov advisor):
  - Same split/top-up rezoning, applied "after every time-step."
  - Δx = 5×10⁻⁵ m, Δt = 0.01 s. Convergence was verified by insensitivity to a factor-of-four change in the integration parameters.
  - Swelling of intumescent PVC and PC was fitted by prescribing component densities, with radially varying densities in 2D.

  Sources: [Swann 2019 PhD thesis, UMD](https://drum.lib.umd.edu/bitstream/1903/21906/2/Swann_umd_0117E_19758.pdf); [Swann et al. 2019, IJHMT (PVC)](https://doi.org/10.1016/j.ijheatmasstransfer.2018.12.011).

**Gpyro (Lautenberger; 1D/2D/3D)**
- **Cell model.** Each cell "is permeable to gaseous mass transfer but is impermeable to condensed phase mass transfer." Δz is an unknown solved from (ρ̄Δz)_P^{n+1} − (ρ̄Δz)_P^n = −ω̇‴_fg Δz Δt plus the species equations. "Equations 3, 56, and 57 are M+2 equations for M+2 unknowns in each grid cell: the M species' mass fractions, ρ̄, and Δz" — [Gpyro Tech Ref §4.1.1–4.1.3](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/Gpyro/gpyro_technical_reference.md).
- **Swelling and shrinkage via χ_k.**
  - χ_k = 1: no volume change (charring).
  - 0 ≤ χ_k < 1 with ρ_B < ρ_A: intumescence.
  - 1 < χ_k < ρ_A/(ρ_A − ρ_B): shrinkage.
  - Noncharring reactions force surface regression.

  Source: [Gpyro Tech Ref §3.3](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/Gpyro/gpyro_technical_reference.md).
- **Burnout.**
  - A positivity-preserving update: (ρ̄Δz)_P = (ρ̄Δz)°/Δt / (1/Δt + ω̇Δz/max(ρ̄Δz, ε)).
  - "if a cell's thickness decreases below a user-specified value, then all reaction rates in that cell are set to zero and that cell behaves as if it was inert."
  - "Although no numerical problems are encountered as the cell size gets large, the solution is not necessarily accurate if gradients are not well-resolved due to large grid sizes."

  Source: [Gpyro Tech Ref §6.4, Eq. 141](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/Gpyro/gpyro_technical_reference.md).
- **The threshold** is keyword `EPS`, default 1×10⁻¹⁰ m ("found to work well"). `NOCONSUMPTION=.TRUE.` disables regression; it was used in 2D FDS-coupled flame spread over thermoplastics — [Gpyro User's Guide](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/Gpyro/gpyro_users_guide.md).
- **Numerics.** Segregated, fully implicit, iterative (relaxation α, TDMA, residual-based convergence). The 3D formulation is labelled "Three-dimensional (no volume change)" — [Gpyro Tech Ref §4.3, §6](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/Gpyro/gpyro_technical_reference.md).
- **Lautenberger, MaCFP workshop slides, June 2017.**
  - Solid momentum in 3D is "N/A (no movement/shrinkage of solid phase in 3D)."
  - Listed limitations: "Particle shrinkage / swelling not accounted for in 2D/3D formulation (only in 1D)", and in FDS coupling "Geometry is static – no shrinkage; Objects can't 'burn away'."
  - The grid is Cartesian, with curved surfaces "stair stepped."

  Sources: [Lautenberger, MaCFP 2017 slides](https://iafss.org/wp-content/uploads/MaCFP_Lautenberger.pdf); [Lautenberger 2014, Gpyro3D, FSS 11](https://publications.iafss.org/publications/fss/11/193).
- **Gpyro V2** (optimized; Abdelhafiz, Amokrane, Bourbigot, Debenest; Fire Safety Journal 2026, IAFSS 2026 special issue):
  - **Grid.** "The computational grid is initially orthogonal and regular, but cells are allowed to deform dynamically in order to capture… swelling or shrinkage."
  - **Mass update.** The iterative implicit mass update is replaced by a direct explicit update with a correction step. Negative masses are clipped by rescaling the destruction rate so the species mass is exactly zero, and the rescaling is propagated to all co-participating species.
  - **Consumed cells.** "If the ratio of the current to the initial cell volume falls below a user-defined threshold (typically 0.01%), the cell is considered consumed… all reaction rates are set to zero and the cell becomes inert."
  - **Thermal solver.** Semi-implicit: implicit in one chosen (depth) direction and explicit laterally, under a Fourier limit. The authors report order-of-magnitude speed-ups.

  Source: [Abdelhafiz et al. 2026, FSJ 104840, doi:10.1016/j.firesaf.2026.104840 (local copy)](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/Gpyro/GPYRO_optimized.md). The paper does not show a multi-D swelling case.
- **Gpyro intumescence validation.** Lautenberger & Fernandez-Pello (2009) applied Gpyro to "gasification and swelling of an intumescent coating" alongside PMMA, white pine and PU-foam smolder — [Lautenberger & Fernandez-Pello 2009, FSJ](https://doi.org/10.1016/j.firesaf.2009.03.011).
- **A Lille group added a Gpyro swelling model** for EVA/ATH that "doesn't need the initial and final material densities." This is known only from a search-engine summary of the HAL record; the full text was blocked by bot protection and is not verified — [HAL hal-04123109](https://hal.univ-lille.fr/hal-04123109).

**Rein-group / Gpyro-derived smouldering models**
- Huang & Rein built 1D peat smouldering models "in the open-source code Gpyro," including ignition, spread and extinction — [Huang & Rein 2014, PROCI](https://doi.org/10.1016/j.proci.2014.05.048); [Huang et al. 2015, IJWF](https://doi.org/10.1071/wf14178). They also note that peat "volume… expands when absorbing water," which changes organic density — [Huang & Rein 2017, IJWF](https://doi.org/10.1071/wf16198).
- A later 2D Gpyro smouldering model exists, but its abstract does not mention volume change — [Yuan et al. 2023, Fuel](https://doi.org/10.1016/j.fuel.2023.128064). This is consistent with Gpyro's 2D/3D having no shrinkage (see above).
- A separate 2D axisymmetric peat model explicitly "considers volume shrinkage" — [Yang et al. 2016, Energy & Fuels](https://doi.org/10.1021/acs.energyfuels.6b02293). Its discretization is not given in the abstract.

**OpenFOAM / FireFOAM `reactingOneDim` pyrolysis region model**
- **Header note.** "The moveMesh option can only be applied to solid reaction such as: PMMA -> gas at the moment" — [OpenFOAM reactingOneDim.H API](https://www.openfoam.com/documentation/guides/latest/api/reactingOneDim_8H_source.html).
- **Implementation** (OpenFOAM-6/7 Foundation source):
  1. Solve ∂ρ/∂t = −RRg, with a mesh-flux term if the mesh moves.
  2. Compute newV = mass0/ρ per cell.
  3. Call `moveMesh(V − newV, minimumDelta)`.
  4. Inside `moveMesh`, the displacement deltaV/|S_f|·n accumulates cell by cell along each 1D column, starting from the gas-coupled patch. Points move only if the resulting cell thickness stays above `minimumDelta` (default 1e-4 m). Cells that did not move are flagged non-reacting via `setCellReacting(i,false)`.
  5. The gas-coupled face stays fixed and the opposite face moves; no cell is ever removed or remapped.

  Sources: [OpenFOAM-7 reactingOneDim.C L214–266, L405](https://github.com/OpenFOAM/OpenFOAM-7/blob/master/src/regionModels/pyrolysisModels/reactingOneDim/reactingOneDim.C); [OpenFOAM-6 regionModel1D.C L180–280](https://github.com/OpenFOAM/OpenFOAM-6/blob/master/src/regionModels/regionModel/regionModel1D/regionModel1D.C).

**Di Blasi (wood and intumescent)**
- **Shrinking wood (1996).** A biomass model in which "variation of the reacting medium properties and particle shrinkage are also described," with shrinkage parameters varied — [Di Blasi 1996, Chem. Eng. Sci.](https://doi.org/10.1016/s0009-2509(96)80011-x). The abstract does not state the grid/coordinate treatment.
- **Intumescent coatings (2001, 2004).** Di Blasi & Branca: "Mass and heat transfer take place across a variable volume medium, including a simplified mechanism for bubble dynamics and material swelling" — [Di Blasi & Branca 2001, AIChE J.](https://doi.org/10.1002/aic.690471020). The model was later run for 20–200 kW/m² — [Di Blasi 2004, JAAP](https://doi.org/10.1016/j.jaap.2003.10.003).

**Other shrinking-solid treatments in the biomass/fire literature (comparison points)**
- **Bellais et al. 2003.** Modelled shrinkage three ways (uniform, shrinking shell, shrinking cylinders). Uniform shrinkage slows pyrolysis while shell/cylinder models speed it, but "the effect was sufficiently small to be neglected given the uncertainty about some wood physical properties" — [Bellais et al. 2003, Fuel](https://doi.org/10.1016/s0016-2361(03)00062-0).
- **Hagge & Bryden 2002.** Shrinkage is negligible for Bi < 10 but matters in the thermal-wave regime: a 2-cm poplar particle at 2000 K with shrinkage factor 0.4 had pyrolysis time cut by 43% — [Hagge & Bryden 2002, Chem. Eng. Sci.](https://doi.org/10.1016/s0009-2509(02)00167-7).
- **Newer multi-D biomass approaches** use different machinery:
  - Eulerian VOF interface tracking for shrinking particles (Basilisk, with a reported "10% average error in the final particle shape") — [Caraccio et al. 2026, Chem. Eng. J.](https://doi.org/10.1016/j.cej.2026.174297).
  - A volume-averaged "pyromechanics" solid-mechanics deformation model in PATO/OpenFOAM — [Lahouze et al. 2025, Fuel](https://doi.org/10.1016/j.fuel.2025.134557).
- **PFEM for melting/dripping thermoplastics in fire** uses "an updated Lagrangian description… Nodes are viewed as material points that can freely move and even separate from the main analysis domain" — [Oñate, Rossi, Idelsohn & Butler 2009, IJNME](https://doi.org/10.1002/nme.2731).

### Inferences
- **The cell concept is shared.** FDS, ThermaKin, Gpyro and OpenFOAM all hold solid mass fixed per cell and derive thickness from density. That is the same solid-Lagrangian concept Pyrolysis.jl uses. The real design axis is the *remediation policy*:
  - FDS: threshold remap with conservative overlap weights.
  - ThermaKin: per-step rezone to a preset Δx (split/top-up).
  - Gpyro/OpenFOAM: never remap, freeze below ε or `minimumDelta`.
- **ThermaKin's integral term matches Pyrolysis.jl's mesh velocity.** ThermaKin's material velocity v(x) = −∫₀ˣ(1/ρ)∂ρ/∂t dx, referenced to the back boundary, is the same kinematics as Pyrolysis.jl's cumulative per-cell volumetric strain from the back face. FDS and OpenFOAM instead hold the exposed face fixed and move the back face.
- **Pyrolysis.jl's parking is cheaper and gentler than the alternatives.** Pyrolysis.jl "parks" burnt-out cells in a fixed-size state vector, while Gpyro and OpenFOAM keep zero-thickness inert cells. Pyrolysis.jl's merge-and-park is the only scheme of the four that both removes the degenerate cell from the physics and keeps a fixed state size. The FDS layer/cell removal changes array length. ThermaKin's split/merge also changes element count every step.
- **None of these codes has error control on the mesh.** All use fixed or heuristically sub-stepped time steps (FDS CN with a ΔT_max sub-step criterion; ThermaKin fixed Δt; Gpyro fixed Δt with iterations). None couples node motion to an adaptive local-error estimate as Pyrolysis.jl does. On error control, Pyrolysis.jl's monolithic approach has no peer in the fire codes.
- **Unidirectional deformation is the standard escape route to 2D.** FDS HT3D's explicit statement that a change in cell count breaks the 3D nodal structure, plus Gpyro's 2D/3D lacking volume change, suggest that cell-count-changing remeshing does not generalize to multi-D in fire codes. ThermaKin2Ds and McDermott's prototype both get to 2D by allowing deformation only along x, column by column.

### Gaps
- I could not confirm how ThermaKin2Ds computes lateral (y/z) fluxes between columns whose element boundaries no longer align after independent x-rezoning. The reference is silent.
- No primary source states whether Gpyro V2 implements cell deformation in 2D/3D; the paper's wording is generic.
- Di Blasi's 1996 shrinking-wood grid treatment and the 2001/2004 intumescent numerics could not be verified: full texts were not accessible, and abstracts do not specify a moving grid or a transformation.
- No published account of Rein-group-specific mesh rules beyond Gpyro defaults was found.

---

## Q2. Intumescent coating models: how do they represent 10–50× expansion, and what resolution/accuracy issues are reported?

### Takeaway
Intumescent models span the full range of approaches:
- Prescribed/empirical expansion with constant-thickness effective properties.
- Front (Stefan) models with coordinate transformation.
- Moving or added-element grids (Lucherini; Kim/Makarov/Molkov).
- ALE/deformed geometry in commercial FE (COMSOL).
- Material (undeformed) coordinates on a fixed mesh (Butler 1995; Ogrin 2017; Sandia/Larsen 2017).
- General pyrolysis codes using density-driven cell growth (Gpyro χ; ThermaKin densities/γ; FDS δ).

The reported accuracy problems are dominated by the *physics of expansion* (the expansion law, its heating-rate dependence, cracking), not by numerics. One quantitative numerical finding stands out: assuming an initially expanded, constant-volume layer underestimates pre-expansion temperatures.

### Cited Findings
- **Expansion magnitudes.**
  - Intumescent coatings "can make the intumescent coatings swell up to 100 times of the initial thickness" ([Cirpici, Wang, Rogers & Bourbigot 2016, Polym. Eng. Sci.](https://doi.org/10.1002/pen.24308)).
  - "expand by a factor of 40 or more" ([Omrane et al. 2007, FSJ](https://doi.org/10.1016/j.firesaf.2006.08.006)).
  - A Sandia-tested commercial coating expanded "to approximately twenty times its initial volume and one-third of its initial mass" ([Larsen et al. 2017, TFEC, SAND2016-9107C](https://www.osti.gov/servlets/purl/1394042)).
- **Classical models with prescribed expansion.** Cagliostro (1975), Anderson (1984/85) and Henderson & Wiecek assumed the maximum expansion and then treated the coating as non-reactive. Di Blasi & Branca and Yuan & Wang also "assumed the final expansions of intumescent coatings were known a priori" — [Cirpici et al. 2016](https://doi.org/10.1002/pen.24308).
  - Anderson & Wauters (1984): control-volume model where "Expansion is explicitly accounted for by assuming it to be a function of mass loss" ([Anderson & Wauters 1984, IJES](https://doi.org/10.1016/0020-7225(84)90036-3)).
- **Front model with transformation (Buckmaster et al. 1986).** All outgassing and swelling occur at a thin front where density, velocity and temperature gradient jump. "Through an elementary mathematical transformation, this problem is reduced to one in which the outer surface is fixed and only the inner front moves," yielding a Stefan problem — [Buckmaster, Anderson & Nachman 1986, IJES](https://doi.org/10.1016/0020-7225(86)90084-4).
- **Butler, Baum & Kashiwagi (NIST, 1995): 3D bubble model in Lagrangian coordinates.**
  - About 10,000 bubble sites nucleate, swell and migrate; flow fields are summed analytically.
  - "As the bubbles nucleate and swell… the model geometry becomes increasingly distorted from the initial rectangular solid… This can be arranged by tracing the locations of convected gridpoints back in time to their original positions through the introduction of a Lagrangian coordinate system. The Lagrangian formulation eliminates the convection term from the energy equation, the complex geometry is incorporated into the diffusion coefficient."
  - The energy equation is solved by finite differences in the original rectangular geometry.

  Source: [Butler, Baum & Kashiwagi 1995, Int. Conf. Fire Research & Engineering](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=911858).
- **Di Blasi & Branca 2001.** Three-component coating with a "variable volume medium" and a simplified bubble-dynamics swelling mechanism. Predicted substrate temperatures agree well with measurements — [Di Blasi & Branca 2001](https://doi.org/10.1002/aic.690471020). Numerics are not given in the abstract.
- **Staggs et al. 2012.** The model is "fully continuous (rather than semi-discrete as in other approaches)." A fraction of evolved gas is trapped and causes expansion. "results from MLC experiments demonstrate a larger than expected range in coating expansion than predicted by the model… the basic expansion mechanism… requires more research" — [Staggs et al. 2012, Chem. Eng. Sci.](https://doi.org/10.1016/j.ces.2011.12.010).
- **Zhang, Wang et al. 2012/2013.** A local expansion law relates "the local rate of expansion… to the local rate of mass loss, rate of temperature change and temperature."
  - Steel temperatures are "typically within 10% and the final expanded thickness within 20%" across furnace curves ([Zhang et al. 2012, FSJ](https://doi.org/10.1016/j.firesaf.2012.02.004); [Zhang et al. 2013, J. Fire Sci.](https://doi.org/10.1177/0734904112453566)).
  - Cirpici et al. 2016 replaced the empirical law with a bubble-growth (Amon–Denson) analytical expansion model ([Cirpici et al. 2016](https://doi.org/10.1002/pen.24308)).
- **Griffin 2010.** Includes "the increase in thermal resistance as gases are formed and the coating expands," gas convection, and radiation across pores; validated in the cone — [Griffin 2010, J. Fire Sci.](https://doi.org/10.1177/0734904109346396).
  - Griffin's five limitations of existing models are summarized by Nørgaard et al.: models are 1D only; data are incomplete; furnace radiation/turbulence is absent; and the expansion-factor basis is questionable — [Nørgaard et al. 2016, AIChE J.](https://doi.org/10.1002/aic.15291).
- **Gillet et al. 2007.** Multi-layer swollen structure under 45 kW solar-furnace fluxes. They note "a better identification of each layer's thermal properties is needed" — [Gillet et al. 2007, J. Phys. D](https://doi.org/10.1088/0022-3727/40/3/030).
- **Ogrin, Saje & Hozjan 2017: 2D FE with material coordinates and element-wise expansion.**
  - Each element "instantaneously expands to its final thickness by the same constant factor of expansion" once its mean Gauss-point temperature first exceeds T_exp. Expansion therefore proceeds element by element.
  - "To alleviate remeshing and updating nodal coordinates due to extension of the coating, the formulation was also adapted for use with material coordinates… Thus, the finite element mesh remains unchanged during the thermal analysis, while the extension of the coating is fully accounted for."
  - Integrals use the undeformed configuration with the deformation-gradient Jacobian J. Validation used 30 coating elements and Δt = 2 s with implicit backward Euler. Convergence was checked by refining elements and the time step.
  - Asaro et al.'s constant-volume (pre-expanded) assumption "showed a substantial underestimation of temperatures prior to the actual expansion."

  Source: [Ogrin et al. 2017, Fire Mater. 42:145](https://doi.org/10.1002/fam.2466).
- **Sandia (Larsen et al. 2017): fixed FE mesh with a local expansion ratio M.** "The FEM implementation of the model uses a fixed mesh to represent the expanding IM by defining a local variable, expansion ratio (M). The effect of dilating the original material space is implemented by adjusting the local effective thermal conductivity" (k*_eff = k_eff/M). M ramps with decomposition to M_max. The authors note it "is not proper when the IM is applied to surfaces with curvature or at intersections of plane surfaces" — [Larsen et al. 2017](https://www.osti.gov/servlets/purl/1394042).
- **Wu et al. 2016.** "The numerical treatment for the expanding coating layer is by considering the expansion effects on the property parameters, instead of dealing with the variation of grid system" — [Wu et al. 2016 (Consensus record)](https://consensus.app/papers/details/48df3ecba6ff5900a74afab07dae087e/).
- **Lucherini, Hidalgo, Torero & Maluk 2022/2023 (IJTS).**
  - Crank–Nicolson FD. "The coating swelling is implemented… by… expanding the mesh representing the physical domain in proximity to the substrate-coating interface."
  - The UQ manuscript version states swelling is implemented "by adding finite elements in the proximity of the substrate-coating interface following empirical correlations," with constant char properties.
  - "accuracy is highly influenced by the coating thickness ahead of the reaction zone… the coating swelling rate plays a key role, while the thermo-physical properties… have a secondary effect." The model is quasi-steady and less accurate before swelling onset and at low heat flux.

  Sources: [Lucherini et al. 2022, IJTS 107922](https://doi.org/10.1016/j.ijthermalsci.2022.107922); [Lucherini et al. 2020, UQ eSpace](https://doi.org/10.14264/9051b1e).
- **Lucherini et al. 2024 review (JSFE).** Compares analytical, FD and FE approaches. It finds "the tendency to simplify the problem into effective thermo-physical properties" and a "lack of homogenisation of modelling approaches," and calls for a universal method — [Lucherini et al. 2024, JSFE](https://doi.org/10.1108/jsfe-10-2023-0038). The full text was not accessible (SSL/Cloudflare blocks), so its per-model mesh classification was not extracted.
- **Kim, Makarov & Molkov 2017.** 1D intumescent paint on a hydrogen tank: "A moving mesh, that is employed for a better prediction of the expansion process… is based on the local changes of heat and mass." Validated against Cagliostro et al. 1975 — [Kim et al. 2017, IJHE](https://doi.org/10.1016/j.ijhydene.2016.02.157).
- **COMSOL ALE (Nyazika et al. 2018).** A silicone/expandable-graphite coating was handled "using the Arbitrary Lagrangian-Eulerian Method (ALE) implemented in COMSOL," assuming homogeneous 1D expansion. It failed to capture the temperature rise after 450 s "because of cracking" — [Nyazika et al. 2018, J. Phys. Conf. Ser.](https://doi.org/10.1088/1742-6596/1107/3/032012).
- **Other FE models with moving boundaries.** Kang et al. 2017 is an FE model of an inorganic intumescent with "expansion of its volume" and "moving boundaries" — [Kang et al. 2017, J. Fire Sci.](https://doi.org/10.1177/0734904117701765). Dedey et al. 2025 is a 3D pyrolysis model using "a moving boundary condition (expansion rate) to mimic the thermal expansion" — [Dedey et al. 2025, J. Fire Sci.](https://doi.org/10.1177/07349041251321408).
- **Adaptive FD (Piloto, Mesquita & Vaz 2014).** A conference paper with keywords "Adaptive mesh," "Moving boundary." A search summary claims finite differences with adaptive space refinement and the method of lines, but the repository page does not show this — [Piloto et al. 2014](https://portaberta.ipb.pt/display/cv-t-8367f056c39be801f61ccc0b2adf6238).
- **ThermaKin for intumescent coatings.**
  - A 9-reaction IFRC model in ThermaKin: "A swelling ratio of 1 combined with the optimized mass transportation coefficient at upper boundary well represents the expansion." Coating-thickness evolution was predicted with R² = 0.94 ([Yi et al. 2024, J. Build. Eng.](https://doi.org/10.1016/j.jobe.2024.110127)).
  - Swann's ThermaKin2Ds 2D models "reproduce the experimental sample shape profiles" for intumescent polymers, but under-predicted rapid early swelling in some cases ([Swann 2019 thesis](https://drum.lib.umd.edu/bitstream/1903/21906/2/Swann_umd_0117E_19758.pdf)).
- **Simplified engineering models avoid resolving expansion altogether.** Tian et al. treat the virgin layer as thermally thin and the char as having negligible heat capacity, expressing insulation as a flux ratio versus char height — [Tian et al. 2018, FSJ](https://doi.org/10.1016/j.firesaf.2018.04.010).
- **Pyropolis (Snegirev et al. 2013).** For charring polymers "material intumescence is assumed to be controlled by the amount of char produced… given the intrinsic char porosity." Expansion propensity is a key sensitivity parameter — [Snegirev et al. 2013, FSJ](https://doi.org/10.1016/j.firesaf.2013.03.012).
- **Aerospace analogue.** Thermo-mechanical swelling of charring ablators was co-simulated with an ALE method in SAMCEF/AMARYLLIS — [Pinaud, swelling ablators (Consensus lists 2026; DOI suggests 2012 Ablation Workshop)](https://doi.org/10.13023/ablat.2012.26).

### Inferences
- **Why material coordinates fit 10–50× expansion.** They avoid the Lagrangian failure mode Pyrolysis.jl suffers, where cells stretch and resolution degrades in physical space. The mesh stays fixed in the reference configuration, and expansion enters only as a metric factor (k/M, J). The trade-off is that resolution per unit *physical* length still falls as 1/M in the char. That is acceptable if the char is near quasi-steady conduction (Lucherini's finding), but not for resolving a moving reaction front in physical space.
- **Pyrolysis.jl's 1D stretching problem is in principle the same as Ogrin's/Sandia's fixed-mesh material-coordinate scheme.** Both keep one cell per material parcel. Material coordinates do not add resolution where the char has expanded. They change only the bookkeeping, unless combined with h-refinement in reference space.
- **Lucherini's sensitivity finding sets the resolution target.** Accuracy hinges on the thickness ahead of the reaction zone and the swelling rate. The fine resolution matters at the reaction/expansion front and at the virgin layer, not inside the fully expanded char. A Lagrangian mesh concentrated in the virgin material, which Pyrolysis.jl's mesh naturally is, already resolves the right place.
- **Threshold-triggered expansion is non-smooth.** Element-wise instantaneous expansion (Ogrin) and front models (Buckmaster) introduce switching on T_exp, so outputs are not smooth in parameters. Continuous kinetically driven swelling (Gpyro χ, ThermaKin γ/densities, Pyrolysis.jl per-component swelling) is smoother.

### Gaps
- I could not obtain the numerical details (grid motion rules, node insertion, energy conservation on insertion, convergence) for these:
  - Lucherini 2022 (UCL/UQ full text blocked).
  - Kim/Makarov/Molkov 2017.
  - Di Blasi & Branca 2001 and Di Blasi 2004.
  - Griffin 2010.
  - Gillet 2007.
  - Zhang 2012.
  - Kang 2017.
  - Staggs 2012.
- No primary source was found for "Mamleev," "Hehnen" or "Elliott" intumescent *models*.
  - The Hehnen hit is an FDS PMMA inverse-modelling paper ([Hehnen et al. 2023, FSJ](https://doi.org/10.1016/j.firesaf.2023.103926)), not intumescence.
  - Elliott/Maluk-group items I recall are experimental, but this is not verified.
  - Mamleev's work appears to be kinetics, but this is not verified.
- No intumescent study found reports grid-convergence specifically for 20–100× expansion on a Lagrangian/moving mesh.

---

## Q3. Does anyone use mass-based (Lagrangian mass / von Mises) coordinates for charring or intumescent solids, and what did it buy them?

### Takeaway
Yes, in two forms:
- **Explicit material (undeformed) coordinates on a fixed mesh** for intumescent coatings: Butler et al. 1995 (3D FD), Ogrin et al. 2017 (2D FE) and Sandia/Larsen et al. 2017 (FE, k/M scaling). Their stated benefit is *no remeshing or node updates*: "the finite element mesh remains unchanged… while the extension of the coating is fully accounted for." Butler adds that the convection term disappears from the energy equation.
- **Implicit mass-Lagrangian cells** in the general codes: Gpyro's impermeable-to-solid cells with Δz unknown, ThermaKin's mass-descriptor elements, and FDS's mass-conserving cell scaling. These give exact per-cell solid-mass conservation, but each code still needs a remedy for degenerate cells (freeze, remap or rezone).

I found no fire paper that uses a cumulative-mass (von Mises) coordinate as the independent variable and reports an accuracy/cost comparison against moving grids.

### Cited Findings
- **Butler et al. 1995.** Lagrangian coordinates trace grid points "back in time to their original positions." This "eliminates the convection term from the energy equation, the complex geometry is incorporated into the diffusion coefficient," and FD is done "in the original rectangular geometry" — [Butler, Baum & Kashiwagi 1995](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=911858).
- **Ogrin et al. 2017.**
  - "To perform the whole analysis on the initial undeformed finite element mesh, we must use the initial undeformed (ie, material) coordinates, rather than the deformed (ie, spatial) coordinates." The deformation map must be "one-to-one, locally invertible… with a positive Jacobian determinant."
  - Expansion was restricted to the direction normal to the substrate, with the J factor applied only in that direction.
  - Benefit: "to alleviate remeshing and updating nodal coordinates… the finite element mesh remains unchanged," which they describe as "particularly convenient for… variously shaped steel cross sections."

  Source: [Ogrin et al. 2017](https://doi.org/10.1002/fam.2466).
- **Sandia/Larsen et al. 2017.** A fixed mesh plus local expansion ratio M, with k*_eff = k_eff/M (the 1D material-coordinate flux form). Caveat: it is not proper at curvature or corners — [Larsen et al. 2017](https://www.osti.gov/servlets/purl/1394042).
- **Buckmaster et al. 1986.** A transformation that fixes the outer surface reduces the two-moving-boundary problem to a single-front Stefan problem — [Buckmaster et al. 1986](https://doi.org/10.1016/0020-7225(86)90084-4).
- **Gpyro.** A cell is "impermeable to condensed phase mass transfer," and Δz is solved per cell from (ρ̄Δz) conservation. This is mass-Lagrangian by construction, with a fixed cell count — [Gpyro Tech Ref §4.1](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/Gpyro/gpyro_technical_reference.md).
- **ThermaKin.** Elements are "characterized by masses of components and temperature… primary object descriptors," with a rezone after each step — [ThermaKin2Ds Reference v2 §3](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/ThermaKin/ThermaKin2Ds%20Reference%20v2.md).
- **Lagrangian/fixed-domain coordinates in regression and ablation analyses:**
  - Perturbation solutions for charring ablators "in the Lagrangian coordinate system because of surface recession and moving boundaries" ([Arabgol et al. 2018, Fire Mater.](https://doi.org/10.1002/fam.2641)).
  - A regressing thermoplastic "transformed to a fixed-domain coordinate system" ([Wichman et al. 2026, Fire Mater.](https://doi.org/10.1002/fam.70036)).
  - A "movable domain model… transformed to a fixed domain model" for charring cable insulation ([Lin et al. 2005, Mech. Adv. Mater. Struct.](https://doi.org/10.1080/15376490490491990)).
- **The same Eulerian→Lagrangian (stationary-mesh) transformation is used for strongly deforming porous biopolymers outside fire.** Solutions are obtained "using a stationary mesh" and converted back, with ∇_E = j^(−1/3)∇_L — [Zhao et al. 2022, J. Food Sci.](https://doi.org/10.1111/1750-3841.16279). It also appears in finite-strain soil consolidation with material coordinates — [Singh et al. 2015, IJNAMG](https://doi.org/10.1002/nag.2393).

### Inferences
- **Pyrolysis.jl already uses mass-Lagrangian cells.** In 1D a material coordinate X (or cumulative-mass coordinate m) is exactly what a solid-Lagrangian cell-centred FV with moving nodes discretizes, so "switching to mass coordinates" would mainly change the state variables:
  - **Option A:** stop integrating z_i as ODE states and recover z(X,t) = ∫J dX algebraically.
  - **Option B:** keep the reference-space mesh fixed and carry J (or per-cell thickness) as a state.

  Either way the Jacobian sparsity simplifies (no node-position coupling). Resolution behaviour in physical space is unchanged.
- **Mass coordinates do not fix burnout.** A cell whose mass goes to zero has zero width in m-space too. Burnout still needs either an inert floor (Gpyro, OpenFOAM) or removal or merging (FDS, Pyrolysis.jl). Strong swelling is likewise not cured: cells keep their mass but grow in physical size, and refinement would still require remapping in X.
- **The established fire-code alternative to remapping is a fixed reference mesh with metric scaling.** Ogrin, Sandia and McDermott's φ_s-scaled distances all avoid discrete remeshing events, which is relevant for smooth inverse problems. It extends to 2D/3D naturally when deformation is unidirectional, or isotropic through J, without changing the cell count. FDS's HT3D failure mode, a changing cell count, never arises.

### Gaps
- No fire-science paper was found that quantifies what material or mass coordinates buy in accuracy per cost versus a moving or rezoned grid.
- Henderson & Wiecek's (1987) expanding-composite coordinate treatment is cited by others but was not verified from the primary source.
- Whether Staggs 2012's "fully continuous" model uses material coordinates could not be confirmed.

---

## Q4. Are there studies comparing mesh strategies (fixed vs moving vs rezoned) for pyrolysis accuracy or cost, or reporting grid-convergence problems during burnout?

### Takeaway
There are only a few direct comparisons, and most are old:
- The Ghent group's moving-grid charring model: Theuns, Vierendeels & Vandevelde 2002, 2004 and 2005.
- Ogrin et al.'s spatial-versus-material formulation with a convergence check.
- McDermott et al.'s 1D-regressing versus 3D-burn-away verification.

The burnout regime is handled everywhere by ad hoc thresholds: FDS layer removal at 1e-6 m, Gpyro `EPS` = 1e-10 m, Gpyro V2 at V/V₀ < 0.01%, OpenFOAM `minimumDelta` = 1e-4 m. None of these is accompanied by a published grid-convergence study of the burnout transient. The code manuals themselves document renoding-induced rate swings (FDS) and loss of accuracy with large cells (Gpyro).

### Cited Findings
- **Ghent moving-grid model (Theuns et al. 2002).** A 1D model in which "the char and virgin grid move along with the pyrolysis front." Non-uniform virgin-layer grids ran about 3× faster than uniform grids, according to the abstract (via search summary). Compared with an integral model, the latter showed significant errors under sudden boundary changes — [Theuns, Vierendeels & Vandevelde 2002, IJNMHFF 12:541](https://doi.org/10.1108/09615530210434287).
- **Theuns et al. 2005.** "Unlike the integral model, the moving grid model does not require any assumption for the temperature profile… the integral model has problems with sudden variations of the external heat flux" — [Theuns et al. 2005, Fire Mater.](https://doi.org/10.1002/fam.878).
- **Jia, Galea & Patel 1999.** "Standard numerical techniques to solve moving front problems… encounter difficulties when dealing with nonlinear boundaries." They propose a decoupling technique verified against exact solutions — [Jia, Galea & Patel 1999, Fire Mater.](https://onlinelibrary.wiley.com/ai/10.1002/(SICI)1099-1018(199903/04)23:2%3C71::AID-FAM672%3E3.0.CO;2-5).
- **Ogrin et al. 2017.**
  - Convergence was studied "by increasing the number of finite elements or/and decreasing the time step," and adequate ranges were assessed.
  - Both deformed-configuration (spatial) and undeformed (material) formulations were implemented.
  - A constant-volume (pre-expanded) representation underestimates pre-expansion temperatures.

  Source: [Ogrin et al. 2017](https://doi.org/10.1002/fam.2466).
- **McDermott et al. 2018.**
  - Compared the FDS 1D solver (exposed face fixed, solid "shrinks from the bottom up," no burn-away) with the 3D solver (cells burn away; mass shifted to a neighbour before removal; final cell clipped at a numerical threshold) on PMMA at 50 kW/m² and 1 cm, with matching δz.
  - They note "Rate of pyrolysis varies between the methods due to the surface heat flux boundary condition and the difference in grid resolution."

  Source: [McDermott et al. 2018](https://files.thunderheadeng.com/femtc/2018_d2-07-mcdermott-paper.pdf).
- **FDS documentation of renoding artefacts.** Renoding "is costly if it is done frequently," and combining cells across large temperature differences "can cause large swings in the pyrolysis rate." The default solid grid can produce "spurious fluctuations in temperature and burning rate" — [FDS User Guide §8.3.8](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/FDS/fds_user_guide.md).
- **Gpyro documentation.** Burnout cells are made inert "to prevent numerical difficulties that occur as Δz → 0"; `EPS` "can reduce convergence problems that may be encountered as grid cells become very small"; and "the solution is not necessarily accurate if gradients are not well-resolved due to large grid sizes" — [Gpyro Tech Ref §6.4](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/Gpyro/gpyro_technical_reference.md); [Gpyro User's Guide](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/Gpyro/gpyro_users_guide.md).
- **Gpyro V2.** A nearly consumed-cell stabilization was needed "in strongly non charring polymers such as PMMA" (threshold 0.01% of initial volume) — [Abdelhafiz et al. 2026 (local)](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/Gpyro/GPYRO_optimized.md).
- **ThermaKin.** Convergence is demonstrated by sensitivity to Δx and Δt (factor-of-four change). "The accuracy of the results gradually decreases with increasing value of either of these parameters" — [ThermaKin2Ds Reference v2 §3](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/ThermaKin/ThermaKin2Ds%20Reference%20v2.md); [Swann 2019](https://drum.lib.umd.edu/bitstream/1903/21906/2/Swann_umd_0117E_19758.pdf).
- **Shrinkage-representation studies (physics, not mesh).** Different shrinkage idealizations (uniform, shell, cylinders) shift pyrolysis rate in opposite directions — [Bellais et al. 2003](https://doi.org/10.1016/s0016-2361(03)00062-0). Shrinkage matters only in the thermal-wave regime — [Hagge & Bryden 2002](https://doi.org/10.1016/s0009-2509(02)00167-7).
- **Inverse-fit dependence on discretization.** Bench-scale inverse-modelling setups for PMMA "differing in cell size and targets… all lead to similar and good results. Yet, they show significantly different performance in the real-scale" — [Hehnen et al. 2023, FSJ](https://doi.org/10.1016/j.firesaf.2023.103926). This likely refers to gas-phase cell size, since FDS solid gridding is separate; I did not verify which.

### Inferences
- **Pyrolysis.jl has little published competition here.** The literature offers no quantitative accuracy-per-cost benchmark for burnout across mesh strategies. Pyrolysis.jl's own verification suite, with Stefan and Landau ablation exact solutions and code-to-code comparisons, is already more systematic than anything published for fire codes on this regime.
- **Every code treats the final cell's disappearance as a singular event.** The universal use of thresholds (freeze, clip, remove) confirms this. Pyrolysis.jl's cost spike after each park is consistent with every other code needing special handling there. The difference is that adaptive-step codes pay for it in step-size collapse, while fixed-step codes pay in unquantified local error.

### Gaps
- No published grid-convergence or order-of-accuracy study of the burnout transient (last-cell disappearance, MLR tail) was found for FDS, Gpyro, ThermaKin or OpenFOAM.
- MaCFP condensed-phase comparisons (PMMA) might contain burnout-tail discrepancies attributable to numerics. I did not find a source that separates numerical from model-form differences.
- The Theuns 2002 full text was not accessed. The 3× speed-up claim comes from the abstract via search summary.

---

## Q5. Is there any treatment that keeps the discrete solution smooth in parameters for inverse analysis (no discrete merging events)?

### Takeaway
I found no fire-pyrolysis source that explicitly designs the shrink/swell/burnout discretization for parameter-smoothness. Inverse analysis in fire practice is dominated by gradient-free optimizers (genetic algorithms, SCE, manual or hierarchical fitting), which tolerate non-smooth forward maps. That is probably why the issue goes undiscussed. The non-smooth elements in existing codes are explicit:
- FDS threshold renoding and layer removal.
- ThermaKin per-step split/top-up and Monte-Carlo radiation absorption.
- Gpyro/OpenFOAM freeze thresholds.
- Gpyro V2's negative-mass clipping.
- Ogrin's element-wise expansion trigger.

Gradient-based or adjoint inference of pyrolysis parameters exists only in the ablation (re-entry) community.

### Cited Findings
- **Gradient-free inverse methods in fire pyrolysis.**
  - Gpyro's coupled genetic algorithm ([Lautenberger & Fernandez-Pello 2009](https://doi.org/10.1016/j.firesaf.2009.03.011); [Gpyro User's Guide](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/Gpyro/gpyro_users_guide.md)).
  - Rein-group peat kinetics by genetic algorithm ([Huang & Rein 2014, C&F](https://doi.org/10.1016/j.combustflame.2013.12.013)).
  - Shuffled Complex Evolution for an intumescent LDPE/APP/PER system ([Coimbra et al. 2023, FSJ](https://doi.org/10.1016/j.firesaf.2023.103813)).
  - A review notes that both gradient-based and direct-search methods are used ([Nyazika et al. 2019, J. Fire Sci.](https://doi.org/10.1177/0734904119852740)).
- **ThermaKin's radiation algorithm is stochastic by option.** The "random absorption algorithm" selects the absorbing element by Monte Carlo at every time step, which makes outputs noisy in parameters — [ThermaKin2Ds Reference v2 §3](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/ThermaKin/ThermaKin2Ds%20Reference%20v2.md).
- **FDS remeshing is explicitly event-based.** It is triggered by `REMESH_RATIO`, `RENODE_DELTA_T` and cell disappearance, and it is suppressed when it would shift pyrolysis rate by more than 15% (`REAC_RATE_DELTA`). Layers are deleted at `MINIMUM_LAYER_THICKNESS` — [FDS wall.f90](https://github.com/firemodels/fds/blob/FDS-6.11.0/Source/wall.f90); [FDS read.f90](https://github.com/firemodels/fds/blob/FDS-6.11.0/Source/read.f90).
- **Gpyro V2's negative-mass correction** rescales destruction rates only "whenever negative masses are predicted" (a switch) — [Abdelhafiz et al. 2026 (local)](file:///home/gbellamy/Documents/git/Pyrolysis.jl/docs/refs/Gpyro/GPYRO_optimized.md).
- **Fixed-mesh metric-scaling formulations** contain no remeshing events, but expansion triggers can be discrete:
  - Ogrin: element expands instantaneously when T̄ > T_exp ([Ogrin et al. 2017](https://doi.org/10.1002/fam.2466)).
  - Sandia: M is a continuous algebraic function of decomposition ([Larsen et al. 2017](https://www.osti.gov/servlets/purl/1394042)).
  - FDS 3D prototype: φ_s is a continuous function of density ([McDermott et al. 2018](https://files.thunderheadeng.com/femtc/2018_d2-07-mcdermott-paper.pdf)).
- **Gradient/adjoint inference in pyrolysis-type problems (ablation, not fire):**
  - A discrete adjoint of a pyrolysis (PICA heat-shield) model supplies gradients for Hamiltonian Monte Carlo ([Bosco 2019, MSc thesis, ULiège](https://matheo.uliege.be/bitstream/2268.2/7860/7/AbstractThesis.pdf)). The abstract does not describe the grid.
  - Automatic differentiation (Tapenade) on the Monopyro ablation/pyrolysis code was used for surface-flux estimation, with results described as "first encouraging." This is per search summary and not opened — [Alestra et al., ICIPE slides](https://www.imo.universite-paris-saclay.fr/~francois.dubois/travaux/vrac/alestra08/alestra-icipe-21janv.pdf).
  - Sandia reports kinetic calibration "can successfully be solved using classical gradient-based optimization," also per search summary — [Sandia 2021 publication record](https://www.sandia.gov/research/publications/details/efficient-kinetic-thermal-inverse-modeling-for-organic-material-decompositi-2021-07-01/).

### Inferences
- **Most promising smooth alternatives to merge-and-park, from the formulations surveyed:**
  - **(a) Fixed cell count with an inert floor** (Gpyro `EPS`, OpenFOAM `minimumDelta`), but implemented as a *smooth* rate switch-off, e.g. reaction rates multiplied by a smooth function of cell mass or thickness, instead of a hard cutoff. This keeps the state size constant with no discrete event. Its stiffness depends on how the vanishing cell's conduction coupling is regularized.
  - **(b) A fixed reference-space mesh with metric scaling** (material coordinates; φ_s-scaled distances). This removes node-position states and remeshing events entirely and extends to unidirectional 2D/3D. A burnt cell would still have to be handled by a smooth floor as in (a).
- **Pyrolysis.jl's discrete parks are the obstacle for gradient-based inversion.** The step-size restart after each park and the merge itself make the discrete map piecewise in parameters, because park times move with parameters. Pyrolysis.jl's monolithic adaptive implicit formulation is otherwise unusually well suited to sensitivity analysis, since forward sensitivities and adjoints through the ODE solver are available. No fire code offers an established smooth alternative to copy; option (a) is the closest existing practice.
- **The adaptive integrator itself can add noise.** Even without remeshing events, error-controlled step selection makes outputs only piecewise-smooth in parameters unless tolerances are tight or sensitivities are computed by differentiating the continuous ODE (forward or adjoint sensitivity). This is a general numerical-analysis point, not one sourced from the fire literature.

### Gaps
- I found no fire-science paper that evaluates objective-function smoothness or gradient accuracy as a function of the solid-mesh strategy.
- No fire pyrolysis code was found that uses AD or adjoints together with a deforming or burning-out solid mesh.
- Whether ThermaKin's default ("maximum absorption") mode is used in CAPA/cone inverse fits, and how much non-smoothness its rezoning introduces, is not documented.
