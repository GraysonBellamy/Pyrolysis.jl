# Keep the Lagrangian mesh, retire the merge

Pyrolysis.jl is not missing a better mesh. The survey covered fire codes (FDS, ThermaKin, Gpyro, OpenFOAM), sixty years of ablation codes (CMA through CHAR, PATO and Icarus), and the fixed-grid, level-set, diffuse-interface, cut-cell and particle literature. No family beats its solid-Lagrangian, monolithically integrated moving mesh on 1D accuracy per cost, error control or surface-boundary fidelity. Every mainstream fire code already uses the same mass-carrying cell. None of them puts node motion under adaptive error control. Converged ablation codes agree within 1–3% whatever their mesh strategy. What Pyrolysis.jl lacks is narrower. Burnout costs too much because of the merge itself, a conservative remap that knocks the keeper cell off its quasi-steady profiles, so no tuning of the step-size reset can remove that cost. The literature offers one form that removes such events without creating voids. It is a **continuous-ALE mesh velocity w = v_s + c**: a small, smooth correction that equidistributes the cells so that none collapses while its layer still has thickness. On the PMMA burnout benchmark this turns about 60 park events into **one terminal event per vanishing layer**, keeps all mesh motion inside the error-controlled ODE, and is the 1D precursor of the pseudo-elastic mesh smoothing that every multi-D ablation code uses. It has two costs: relative-velocity fluxes that smear composition wherever the correction acts, and no published precedent through burnout, so it must pass a prototype gate. A second, lower-risk change is to carry **per-cell thickness instead of node positions**. This is the material-coordinate form that Gpyro and the intumescent-coating models use, and it removes the prefix-sum coupling that makes the exact Jacobian 60–90× more expensive per step. For the 1D core, reject Lagrange-remap, Landau front-fixing, h-adaptivity and fixed-grid interface methods. For 10–50× intumescence, no published method resolves the expanded char better, so measure before building anything.

## Every production code shares Pyrolysis.jl's cell but none shares its integrator

The fire codes all use a *mass-carrying, solid-Lagrangian cell*: solid mass never crosses a face, and thickness follows from mass divided by density. They differ only in what they do when cells become too large or too small:

- **FDS** rescales every cell each step. It renodes when a cell deviates from its last remeshed size by more than `REMESH_RATIO = 0.25`, or when blending would move the temperature by more than the change that alters a reaction rate by 15% (`REAC_RATE_DELTA = 0.15`). The renode is an overlap-weighted conservative remap, and layers thinner than about 1 µm are deleted ([FDS wall.f90](https://github.com/firemodels/fds/blob/FDS-6.11.0/Source/wall.f90); [FDS read.f90](https://github.com/firemodels/fds/blob/FDS-6.11.0/Source/read.f90)). The user guide warns that renoding "is costly if it is done frequently" and "can cause large swings in the pyrolysis rate" ([FDS User Guide §8.3.8](../docs/refs/FDS/fds_user_guide.md)).
- **ThermaKin** rezones before every time step to a preset Δx (5×10⁻⁵ m suggested), splitting and topping up elements, with a fixed Δt of 0.01 s ([ThermaKin2Ds Reference v2](../docs/refs/ThermaKin/ThermaKin2Ds%20Reference%20v2.md)).
- **Gpyro** solves Δz as an unknown per cell. A cell thinner than `EPS = 1e-10` m is frozen inert ([Gpyro Tech Ref §4.1, §6.4](../docs/refs/Gpyro/gpyro_technical_reference.md)), and Gpyro V2 freezes a cell once it falls below 0.01% of its initial volume ([Abdelhafiz et al. 2026](../docs/refs/Gpyro/GPYRO_optimized.md)).
- **OpenFOAM's** `reactingOneDim` stops moving a cell below `minimumDelta = 1e-4` m and marks it non-reacting ([OpenFOAM-7 reactingOneDim.C](https://github.com/OpenFOAM/OpenFOAM-7/blob/master/src/regionModels/pyrolysisModels/reactingOneDim/reactingOneDim.C)).

ThermaKin's material velocity, integrated from the back face, is the same kinematics as Pyrolysis.jl's cumulative volumetric strain. In fire codes the real design choice is therefore the **remediation policy**, not the cell.

The ablation community has used five families:

- translate-and-drop grids (CMA);
- uniformly contracting or Landau grids (SODDIT, FIAT, Blackwell–Hogan, Amar, MOPAR);
- fixed grids with control-volume deactivation (ITRAC);
- multi-D ALE, where surface recession drives the boundary and a pseudo-elastic, Laplacian, spring or radial-basis-function solve places the interior nodes (CHAR, PATO, KATS, Icarus, Sandia);
- fixed-mesh volume-averaged, level-set or element-death methods (Argo, Sandia).

Sources: [Amar et al. 2022](https://ntrs.nasa.gov/api/citations/20210022987/downloads/AIAA_2022_Mesh_Deformation.pdf); [Sandia SAND2014-2970C](https://www.osti.gov/servlets/purl/1142845). In every production code surveyed, mesh nodes are driven by the surface and do not follow mass; CHAR states explicitly that its nodes are "not attached to a specific point of mass".

Coupling is loose almost everywhere:

- CHAR applies the recession computed at step n during step n+1 and claims "a significant speed-up over a fully implicit approach" ([CHAR overview 2016](https://ntrs.nasa.gov/api/citations/20160005889/downloads/20160005889.pdf)).
- Sandia solves mesh motion in a separate, segregated step.
- PATO solves every equation sequentially ([Lachaud & Mansour 2014](https://jeanlachaud.com/research/Lachaud2014_JTHT_preprint)).

The two data points on tighter coupling both favour it. Amar, Blackwell & Edwards found that putting recession-rate sensitivities into a full Newton "improves convergence and reduces total computation time" ([Sandia 2006](https://www.sandia.gov/research/publications/details/one-dimensional-ablation-using-a-full-newtons-method-and-finite-control-vol-2006-01-01/)). In 3dFIAT, time-lagged coupling produced "severe numerical oscillations" and failed near peak heating ([Chen & Gokcen 2013](https://ntrs.nasa.gov/api/citations/20140005552/downloads/20140005552.pdf)).

When codes are converged, the mesh strategy does not show up in the results:

- 13 of 14 codes in the 2011 Ablation Workshop case agreed with FIAT, with numerical differences of about ±1% ([Lachaud et al. 2011](https://ntrs.nasa.gov/api/citations/20110014340/downloads/20110014340.pdf)).
- PATO and FIAT agree within 2% ([NASA/CR-2015-218960](https://ntrs.nasa.gov/api/citations/20160006963/downloads/20160006963.pdf)).
- ITRAC's fixed grid and CMA's node dropping agree "excellent[ly]" ([Ewing et al.](https://ntrs.nasa.gov/api/citations/20140008557/downloads/20140008557.pdf)).

Strategies therefore differ in cost, in burnout, under strong swelling and in multi-D, and the published benchmarks exercise none of these.

| Family | Representative codes | Shrinkage and burnout | Fit for Pyrolysis.jl |
|---|---|---|---|
| Solid-Lagrangian cells with remediation | FDS, ThermaKin, Gpyro, OpenFOAM, **Pyrolysis.jl** | Threshold remap, per-step rezone, inert freeze, or merge-and-park | Current; exact transport of solid composition |
| Continuous ALE / MMPDE r-adaptation | MOVCOL, Huang–Russell; none in pyrolysis | Monitor floor stops collapse; one terminal event | **Missing; strongest candidate** |
| Material coordinates on a fixed reference mesh | Butler 1995, Ogrin 2017, Sandia 2017, Gpyro's Δz | Metric (J or Δz) as state; burnout still needs a floor | **Worth adopting as a state choice** |
| Landau or contracting grid | SODDIT, FIAT, Amar, MOPAR | Cannot reach zero thickness; spallation impossible | Reject: worse for volumetric pyrolysis |
| Indirect ALE (Lagrange-remap), h-adaptivity | Hydrocodes, ocean models, CHAR's h-AMR | Discrete remap or refine events | Reject: events outside error control |
| Fixed grid with sharp or diffuse interface | Level set, VOF, Argo, PATO-COACO, cut cell | Porosity → 1 or cut-cell agglomeration | 2D topology fallback only |
| Particle and mesh-free methods | Peridynamics, PFEM, SPH | Point deletion; dripping, cracking | Out of scope |

Pyrolysis.jl's distinctive pieces are its node positions held as error-controlled ODE states, its explicit discrete geometric conservation law, and its fixed-size park-and-mask burnout. None of these was found in any surveyed code. Against Gpyro V2 on the reference cone, the project's same-machine benchmark found Pyrolysis.jl 3–10× cheaper at equal accuracy. The current design is ahead of the field everywhere except at the burnout event.

## Burnout cost comes from the merge, not from the step-size reset

Today a cell is parked when it falls below 5% of its initial thickness. It is merged into its interior neighbour with enthalpy conserved, and the integrator is handed a proposed step of 10⁻⁶ s ([callbacks.jl](../src/Solver/callbacks.jl)). The code's own docstring gives the mechanism. The merged keeper cell sits off the quasi-steady gas and temperature profiles of its neighbours. The residual at the parked state is therefore orders of magnitude above its pre-park value, and the extrapolated predictor lands outside the Newton basin. The technical reference adds that `DirectSolve`, restarted at 5×10⁻² s after each park, **fails its Newton more than 500 times** on the 60-park PMMA burnout. The reduced strategy with the reset needs **858 steps and 8 failures** ([Technical Reference Guide](../docs/Technical_Reference_Guide.md)).

A step-control measurement on V-BO-01 attributed **61% of all steps and 47% of Newton iterations to post-park recovery**, roughly nine steps per park. In the same measurement, larger absolute resets nearly doubled the step count and resets relative to the pre-park step failed outright. The cheap experiment the numerics literature suggests, cutting to 0.1–0.5× the pre-event step, has therefore already been run on this code, and it does not help.

The literature shows the same mechanism. Furzeland, Verwer & Zegeling found that regridding between BDF steps caused **259 error-test failures, most of them right after a regrid**. Limiting regridding to about 100 events cut steps from 663 to 331 and Jacobians from 709 to 98 ([Furzeland et al. 1990](https://ir.cwi.nl/pub/1640/1640D.pdf)). Every discrete remap therefore injects the same off-manifold transient: FDS renoding, ThermaKin rezoning, indirect-ALE remap, h-coarsening, and Pyrolysis.jl's merge. Each remap also adds numerical diffusion that falls with reconstruction order but "does not vanish" ([Adcroft, ECMWF 2020](https://events.ecmwf.int/event/167/contributions/1377/attachments/796/1415/AS2020-Adcroft.pdf)). Lagrange-remap is mature in hydrocodes and ocean models, but for this solver it is the wrong direction. It would multiply the very events that already dominate the cost.

The two event-free alternatives that fire codes use also fail here:

- **Letting cells collapse.** Without parking, surface cells shrink to **10⁻¹⁰–10⁻¹³ m**, the step collapses, and the Newton fails on 20–40% of steps with either Jacobian ([Performance plan](../docs/Performance_Improvement_Plan_2026-09-26.md)). Gpyro and OpenFOAM tolerate collapsing or frozen cells only because their fixed-step segregated solvers never see the conditioning. They pay in local error that nobody measures.
- **A finite inert floor.** A floor, the smooth version of Gpyro's freeze, leaves void films where solid used to be. V-BO-05 shows what voids do to the energy path. A void cell conducts heat only at the gas conductivity and absorbs nothing, so a void column under a surviving skin sends the transmitted flux into the substrate "far outside the property bracket" ([V-BO-05](../docs/guide/verification/V-BO-05.md)). Sixty stacked films would build exactly that resistance.

For inverse work, park events make outputs only piecewise smooth in the parameters. Sensitivities exist between events (Galán, Feehery & Barton's jump conditions; [MIT record](https://yoric.mit.edu/?p=420)), but the parameter-to-output map develops kinks whenever a perturbation moves a park across an observation time or changes the park count ([Khan & Barton 2017](https://yoric.mit.edu/biblio/generalized-derivatives-hybrid-systems)). The finite-difference gradients the project now relies on are the most exposed to these kinks. Decks that never park, which includes most charring-material fits, are already event-free under the current mesh.

## A continuous-ALE floor is the one missing option that removes events

The r-adaptivity literature gives Pyrolysis.jl a clear framework. In 1D, simultaneous "quasi-Lagrangian" coupling of the mesh and the physics in a stiff method-of-lines solver "is generally the method of choice". Budd, Huang and Russell recommend that mesh points "move in a similar manner to the particles, but not to follow their motions too precisely" ([Budd, Huang & Russell 2009](https://people.bath.ac.uk/mascjb/Papers09/acta.pdf)). Huang and Kamenski proved that a semi-discrete moving-mesh PDE (MMPDE) mesh stays nonsingular for all time, with cell volumes bounded below, in any dimension ([Huang & Kamenski 2018](https://arxiv.org/pdf/1512.04971)). Large-strain poromechanics uses the same idea, updated-Lagrangian motion plus mesh smoothing, to stop elements inverting ([Sabetamal et al. 2015](https://doi.org/10.1002/nag.2386)). Pyrolysis.jl is already the degenerate case c = 0.

The concrete design has seven parts:

- **Node velocity.** Set ż_i = v_s(z_i) + c_i, with c the relaxation residual of an MMPDE5-type equidistribution (Budd et al.'s relaxed MMPDE; the τ ż_i term sets its response time).
- **Monitor.** Use a floor α, plus a solid-mass term that keeps nodes near-Lagrangian wherever composition varies, plus optionally a reaction-rate term.
- **Pinned nodes.** Pin the back node, move the surface node with the physical surface velocity as now, and keep layer-interface nodes Lagrangian.
- **Soft activation.** Make c ≈ 0 while cells are healthy, so charring decks keep today's exact composition transport. Switch it on smoothly: hard switches are nondifferentiable and stiff.
- **What burnout looks like.** Solid flows into a thinning cell through its faces, so the cell does not collapse. No voids form, because the skeleton velocity, not the mesh, still decides where solid goes, and the layer shrinks exactly as it does now.
- **The endpoint.** At L → 0 every body-fitted method is singular; Landau's transformed advective term diverges, for example ([Ewing et al.](https://ntrs.nasa.gov/api/citations/20140008557/downloads/20140008557.pdf)). With a substrate below, the vanishing layer collapses as a block, giving **one event per layer**. That event comes at the same minimum cell size the solver already handles at every park: 60 cells of 98 µm squeezed into 5% of the layer gives about 5 µm cells.
- **Residue.** Materials that leave char residue need no events at all.

The costs are real:

- **Relative fluxes for every quantity.** Every per-volume quantity needs a relative-velocity face flux (v_s − w)·ξ. The pore gas is carried at the skeleton velocity together with the solids, which is consistent with the project's gas-frame derivation and does not bring back the deleted gas sweep term.
- **Artificial smearing of solid species.** Solid species have zero physical diffusivity, so all the numerical diffusion this flux brings is artificial smearing of char/virgin fronts and layer interfaces. Temperature is far less sensitive: with a correction of about 2×10⁻⁵ m/s, 10⁻⁴ m cells and κ ≈ 10⁻⁷ m²/s, the cell Péclet number is about 0.02. First-order upwinding is therefore fine for T. Species need a limited second-order reconstruction with a *smooth* limiter, because minmod and van Leer reintroduce kinks.
- **Tuning.** The response time τ and the monitor smoothing must be tuned, and "the tuning, being rather problem-dependent, does not lend itself to automation" ([Furzeland et al. 1990](https://ir.cwi.nl/pub/1640/1640D.pdf)).
- **A documented failure.** One group's fully implicit MMPDE coupled to precipitation chemistry failed in ode15i and had to sub-step the mesh equation ([IJNMF 2024](https://doi.org/10.1002/fld.5294)).

Some hooks already exist: the `w_override` path with a limited advection kernel for T driven by the mesh velocity relative to the material ([ale_fluxes.jl](../src/Discretization/ale_fluxes.jl)), and the material-frame sources in the method-of-manufactured-solutions (MMS) driver. The species relative flux and its Jacobian entries would be new.

Turning the monitor's solution-adaptive terms up turns the same machinery into accuracy-driven r-adaptation. The expected gain is modest. Orders of magnitude are documented only for singular structures: **81 moving points beat 10⁵ uniform points** on self-similar blow-up ([Budd et al. 2009](https://people.bath.ac.uk/mascjb/Papers09/acta.pdf)). A moderately thin flame went from 293 to 150 steps ([Furzeland et al. 1990](https://ir.cwi.nl/pub/1640/1640D.pdf)). Pyrolysis fronts resolved by tens to hundreds of cells resemble the flame case.

The other family members fall away:

- **Landau front-fixing** is the same correction chosen to keep spacing uniform. It puts relative flux everywhere, interfaces "awkward[ly]" with substrates, and cannot model spallation ([Ewing et al.](https://ntrs.nasa.gov/api/citations/20140008557/downloads/20140008557.pdf)). It is strictly worse.
- **h-adaptivity** makes the state size variable and nonsmooth. ARKODE can resize one-step integrators without losing their adaptivity heuristics ([Reynolds et al.](https://arxiv.org/pdf/2205.14077)), but the variable size works against fixed sparsity and finite-difference gradients.

The decisive uncertainty is that **no published code documents continuous ALE through burnout to zero residue** under adaptive implicit integration.

## Carrying cell thickness instead of node positions makes the exact Jacobian local

The fire codes that solve shrinkage most cleanly store the metric, not positions:

- Gpyro solves (ρ̄Δz) and Δz per cell, "M+2 equations for M+2 unknowns in each grid cell" ([Gpyro Tech Ref §4.1](../docs/refs/Gpyro/gpyro_technical_reference.md)).
- The intumescent models of Butler, Ogrin and Sandia run on a fixed reference mesh in material coordinates, so "the finite element mesh remains unchanged… while the extension of the coating is fully accounted for" ([Ogrin et al. 2017](https://doi.org/10.1002/fam.2466); [Larsen et al. 2017](https://www.osti.gov/servlets/purl/1394042)). Butler adds that the Lagrangian form "eliminates the convection term from the energy equation" ([Butler, Baum & Kashiwagi 1995](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=911858)).

In 1D, Pyrolysis.jl's Lagrangian mesh already *is* a material-coordinate discretization; only its choice of state differs.

That choice matters for the Jacobian. The ALE residual has three nonlocal couplings ([Performance plan §5.1](../docs/Performance_Improvement_Plan_2026-09-26.md)):

- the node-velocity prefix w_k = Σ_{i<k}(θ_i − θ_A)Δz_i;
- the Beer–Lambert optical-depth suffix;
- the global scalars of the lateral-shrinkage (`WithChi`) mode.

The prefix couples every node to every cell below it, "which makes the factorized W effectively dense and the sparse direct solve 60–90× more expensive per step" ([Technical Reference Guide](../docs/Technical_Reference_Guide.md)). Production therefore uses `ApproxSolve`, an inexact Newton that drops these couplings. Carrying h_i = Δz_i as the state instead gives ḣ_i = (θ_i − θ_A)h_i, which is local. No production cell equation reads w (the residual writes it only into the z rows): T is advected by the relative mesh velocity only under the verification override, and species carry only the dilation term. The prefix therefore disappears entirely rather than being bordered or lifted. Only the physical Beer–Lambert suffix and the `WithChi` scalars remain.

This is a change of variables (z is recovered as a cumulative sum), not the mass-matrix DAE the performance plan rejected. The semi-discrete physics is unchanged. The discrete geometric conservation law holds automatically, a park becomes h_keeper += h_parked, and the hybrid's correction enters locally as c_{i+1} − c_i. A further option is to carry inventories m = ξh, as ThermaKin's and Gpyro's mass descriptors do. That deletes the dilation term and makes mass conservation exact by construction.

The costs are a large refactor and a change in the error-control norm. The refactor touches the state layout, the hand-derived Jacobian operators, parking and the MMS driver. Errors would be controlled on thicknesses rather than positions, so the surface-position error becomes a sum of per-cell errors. The wall-time payoff is also unproven: the reduced Newton has not been found to stall. The sure gains are an exact Jacobian at roughly the reduced cost and the deletion of the two-strategy machinery. None of this fixes burnout or swelling, because a cell whose mass goes to zero has zero width in material space too.

## Intumescence and 2D both point to column-wise material frames, not fixed grids

Intumescent coatings expand "up to 100 times" ([Cirpici et al. 2016](https://doi.org/10.1002/pen.24308)), "by a factor of 40 or more" ([Omrane et al. 2007](https://doi.org/10.1016/j.firesaf.2006.08.006)), or to about 20× their initial volume ([Larsen et al. 2017](https://www.osti.gov/servlets/purl/1394042)). The published models span four approaches:

- prescribed expansion;
- front transformations ([Buckmaster et al. 1986](https://doi.org/10.1016/0020-7225(86)90084-4));
- stretched or added-element grids ([Lucherini et al. 2022](https://doi.org/10.1016/j.ijthermalsci.2022.107922); [Kim et al. 2017](https://doi.org/10.1016/j.ijhydene.2016.02.157));
- COMSOL ALE ([Nyazika et al. 2018](https://doi.org/10.1088/1742-6596/1107/3/032012)) and material coordinates.

**No fixed-grid Eulerian, VOF or level-set intumescent model exists**, and smoothed-particle hydrodynamics (SPH) could not reach realistic 10–100× expansion ([Cirpici 2015](https://research.manchester.ac.uk/en/studentTheses/simulating-the-expansion-process-of-intumescent-coating-fire-prot)).

The reported errors come mostly from the physics of expansion. Lucherini found that accuracy "is highly influenced by the coating thickness ahead of the reaction zone" and by the swelling rate, with char properties secondary ([Lucherini et al. 2022](https://doi.org/10.1016/j.ijthermalsci.2022.107922)). Resolution matters at the front and in the virgin layer, which is where a Lagrangian mesh keeps its cells. Every material-coordinate scheme loses resolution in the char in proportion to 1/M (M being the local expansion ratio), just as Pyrolysis.jl does; nobody has done better.

Resolving expansion dynamically is worth what it costs. Ogrin et al. found that a constant-volume, pre-expanded representation "substantially underestimate[s] temperatures prior to the actual expansion" ([Ogrin et al. 2017](https://doi.org/10.1002/fam.2466)). Their threshold-triggered element expansion produced small non-physical oscillations, whereas kinetic swelling of the kind Pyrolysis.jl uses varies smoothly.

No published study reports grid convergence at 20–100× expansion. Pyrolysis.jl should therefore build that verification case before adding redistribution. Grading the initial mesh toward the swelling layer is a free mitigation. If stretched char cells still limit accuracy, the hybrid's material monitor can redistribute them, and relative flux through homogeneous char smears little.

For 2D, merging has no analogue ([Multi-D assessment](../docs/MultiD_Extension_Assessment_2026-09-26.md)):

- FDS warns that if a solid's cell count changes, "the 3-D nodal structure breaks down" ([FDS User Guide §8.4.2](../docs/refs/FDS/fds_user_guide.md)).
- Gpyro's 2D/3D formulations have no shrinkage ([Lautenberger 2017](https://iafss.org/wp-content/uploads/MaCFP_Lautenberger.pdf)).
- McDermott et al. warn that the choices for multi-D intumescence and burn-away "can have leading-order effects on the model results" ([McDermott et al. 2018](https://files.thunderheadeng.com/femtc/2018_d2-07-mcdermott-paper.pdf)).
- ThermaKin2Ds and McDermott's prototype reach 2D by deforming column by column.

In 2D the strain integral no longer fixes the solid velocity, so a mechanics or compatibility solve becomes necessary. KATS ([Fu 2018](https://uknowledge.uky.edu/me_etds/111)), PATO ([Oruganti et al. 2026](https://ntrs.nasa.gov/api/citations/20260007012/downloads/Oruganti_AblationWorkshopAbstract_2026_FP2.pdf)) and Amaryllis ([Pinaud 2012](https://uknowledge.uky.edu/ablation/2012/Coupling/1)) all couple one in for swelling. The interior placement used in multi-D practice (pseudo-elastic, Laplacian or radial-basis-function solves) is exactly the 1D hybrid's c generalized. CHAR's experience shows its limits: imposing all displacement components on curved or cornered boundaries "can become overly constraining and lead to poor quality or tangled meshes" ([Amar et al. 2022](https://ntrs.nasa.gov/api/citations/20210022987/downloads/AIAA_2022_Mesh_Deformation.pdf)).

Perforation and other topology changes defeat any body-fitted mesh. Sandia's element death "allows arbitrarily large surface deformations" at "modest to poor surface resolution" ([Sandia 2014](https://www.osti.gov/servlets/purl/1142845)). Diffuse surfaces carry measured penalties:

- Standard volume-of-solid formulations misjudge reactive surface area by a factor of 2–4 and miss the ALE reaction constant by up to 78% ([Maes et al. 2022](https://www.frontiersin.org/articles/10.3389/feart.2022.917931/full)).
- Smearing an interface causes "O(1) error in gradients" ([Bayat et al. 2022](https://doi.org/10.1016/j.jcp.2022.111627)).
- A VOF model of a shrinking biomass particle converges only to first order and misses the final shape by about 10% ([Caraccio et al. 2025](https://arxiv.org/html/2510.17588v2)).
- Only strictly conservative cut cells match body-fitted grids ([Başkaya et al. 2023](https://arxiv.org/abs/2309.09657)), and their small-cell problem simply relocates the burnout event.

Fixed grids are therefore a 2D topology fallback, not a 1D upgrade. A 1D-ray material network coupled to multi-D flow is a demonstrated intermediate step ([Zibitsker et al. 2024](https://doi.org/10.1016/j.ijheatmasstransfer.2023.124728)).

## A prioritized path that measures first and keeps exactness where it pays

| Priority | Action | Goals served | Effort | Go/no-go gate |
|---|---|---|---|---|
| 1 | Measure three things: finite-difference gradient smoothness of the V-BO-01 quantities of interest (QoIs) as park times shift; Newton iterations per stage, exact vs reduced Jacobian, on the audit decks; and a 20–50× intumescent grid ladder | Sizes priorities 2–5 | Days | None; this sets the evidence |
| 2 | Re-project the keeper's gas and temperature onto the neighbours' quasi-steady profiles at each park (burnout item D) | Burnout cost | Small | V-BO-01 steps fall and QoIs stay in tolerance; at most about 1.9× |
| 3 | Prototype the continuous-ALE floor: soft-activated MMPDE correction, relative fluxes with a smooth limiter, Lagrangian layer interfaces, parking kept as a fallback | Burnout robustness, smoothness, swelling, 2D precursor | Medium; research risk | V-BO-01 at ≤1 event per layer; peak MLR and its timing within current tolerances; ledgers closed; charring decks unchanged (c ≈ 0); smooth gradients across former park times; wall time no worse |
| 4 | Thickness-state (material-coordinate) reformulation, optionally with inventory states | Exact local Jacobian; deletes the two Jacobian strategies; column-ready | Medium–large | Do it if priority 1 shows exactness cuts Newton work, or bundle it with 3 |
| 5 | 2D as column-wise material frames, hybrid smoothing and deactivation for topology change | 2D path | Months | After 3 and 4; never port merging |

Do not pursue Lagrange-remap or renoding, Landau grids, h-adaptivity, level-set, phase-field, cut-cell or particle methods for the 1D core. Do not adopt the Gpyro-style inert floor either: on an adaptive solver it either collapses cells or builds void films.

## Conclusion

The meshing frontier in this field is not the mesh but the event. Every code surveyed, fire or ablation, fixed-step or implicit, handles the vanishing cell with a threshold. Pyrolysis.jl feels this more sharply only because its adaptive integrator measures what fixed-step codes absorb as unquantified local error. That reframes the burnout gap: it is the price of honest error control, and it is cured by never creating the off-manifold state, not by a cleverer reset. The two recommended changes are complementary. One changes what the mesh does (w = v_s + c); the other changes what the solver stores (h instead of z). Together they would make the 1D core local, light on events and smooth in its parameters, and they are the same pieces a column-wise 2D extension needs.

The open risks are concrete. No continuous-ALE pyrolysis code has burned through in print. Composition smearing has to stay confined to homogeneous vanishing layers. The 20–50× intumescence resolution question has no published answer. The project's verification suite is already more systematic on burnout than anything in the literature, so the evidence that settles these risks will come from the gates above, not from further reading.
