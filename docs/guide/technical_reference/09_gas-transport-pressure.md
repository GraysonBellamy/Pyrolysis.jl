# 9. Gas-Phase Mass Transport and Pressure

Gas mass is stored per bulk cell volume and transported relative to the condensed
skeleton. `LumpedTransfer()` is the default material model. `DustyGas(...)` adds
common-diffusivity molecular friction and species-dependent Knudsen resistance.
Both models use the same pore state, Darcy term, conservative face assembly,
reservoir data and boundary energy convention. Skeleton motion relative to the
mesh is a separate ALE flux. Mobile pore transport is supported in one dimension.

![Bulk and intrinsic pore gas densities, wall and molecular friction, Darcy flow, skeleton motion, and continuous species exchange through an optional film to a complete reservoir state.](../../assets/figures/12_gas_transport.svg)

## 9.1 Pore state, units and assumptions

For bulk concentration $\xi_j=m_j/V$, skeletal density $\rho_{i,j}$, molar mass
$M_j$, and local temperature $T$,

$$
\phi=1-\sum_{j\in S\cup L}\frac{\xi_j}{\rho_{i,j}(T)},\qquad
\rho_j^g=\frac{\xi_j}{\phi},\qquad
p_j=\frac{\xi_jRT}{M_j\phi},\qquad P=\sum_{j\in G}p_j.
$$

Gas mole and mass fractions are $x_j=p_j/P$ and
$Y_j=\rho_j^g/\sum_k\rho_k^g$. At vacuum these fractions are reported as zero.
The EOS is evaluated for every gas-bearing material, including a sealed one.
The central implementation is `pore_state` in `src/Physics/pressure.jl`.

Mass flux $J_j$ has units kg/(m² s), molar flux $n_j=J_j/M_j$ has units
mol/(m² s), and both use bulk face area. Darcy velocity $v_D$ is superficial.
The model assumes ideal pore gases and local thermal equilibrium in the cells.
It excludes Soret diffusion, thermal transpiration, gravity, barodiffusion
corrections and separate gas temperatures. Uniform partial pressures give zero
transport even when temperature and porosity vary. This is a statement about
the specified pressure potential, not a thermal-diffusion model.

## 9.2 Condensed matrix coefficients

Transport weights use condensed bulk occupancy only:

$$
q_j=\gamma_j\xi_j/\rho_{b,j}(T),\quad q=\sum_{j\in S\cup L}q_j,\quad
w_j=q_j/q.
$$

Gas components have zero matrix weight for every value of their swelling
factor. A condensed component with $\gamma=0$ has zero matrix weight but still
occupies skeletal volume in the pore-state equation. Its matrix coefficients
must be zero. For participating condensed components, the lumped model requires
an explicit `λ`; `κ` defaults to zero. Gas components carry neither coefficient.

For $a_j=\lambda_j$ or $\kappa_j$, `PARALLEL` gives $\sum_jw_ja_j$, `SERIES`
gives $(\sum_jw_j/a_j)^{-1}$, and `WEIGHTED` combines these with weight
$\beta$ on the parallel value. A positive-weight zero coefficient closes the
series channel exactly. No positive coefficient floor opens a closed channel.
A disappearing zero-coefficient constituent can make the exact series rule
discontinuous; setup checks reachable support changes and rejects unsupported
paths. `CARMAN_KOZENY` is a separate permeability closure (§5.7).

A retained matrix-free cell requires explicit `LumpedEndpoint` or
`DustyEndpoint` data. The model blends matrix and endpoint coefficients with a
cubic smoothstep in the declared `blend=(lo,hi)` interval of $q$. Endpoint data
include a maximum supported skeletal occupancy. Zero weighted occupancy alone
does not establish an unobstructed void. Without an endpoint the supported
domain ends at $q=hi$, with outcome `:matrix_exhausted`.

## 9.3 Lumped transfer and Darcy transport

The continuum law is

$$
J_j=-\frac{\Lambda M_j}{RT}\nabla p_j+\rho_j^g v_D,
\qquad v_D=-\frac{\kappa}{\mu}\nabla P.
$$

$\Lambda$ is the effective superficial matrix transfer coefficient [m²/s].
At constant temperature and porosity with $\kappa=0$, the diffusion coefficient
for bulk concentration is $\Lambda/\phi$. The material's `gas_viscosity` supplies
$\mu(T)$ [Pa s]. Either wall transfer or permeability may be exactly zero;
a path can transmit total flux by different mechanisms in different layers.

## 9.4 Dusty-gas friction

The explicit inputs are `tortuosity ≥ 1`, positive `D_ref`, `T_ref`, `P_ref`,
and either an explicit pore-diameter law or the effective capillary closure.
The common molecular diffusivity and species wall diffusivities are

$$
D_m=D_{ref}(T/T_{ref})^{7/4}P_{ref}/P,\quad
K_j=\frac{\phi d_p}{3\tau}\sqrt{\frac{8RT}{\pi M_j}},\quad
 t=\frac{\tau}{\phi D_{ref}P_{ref}}(T_{ref}/T)^{7/4}.
$$

Here $t$ is molecular resistance per pressure. With $g_j=-\nabla p_j/(RT)$,
the wall/molecular flux solves

$$
\left(K_j^{-1}+tP\right)n_j^D-tp_j\sum_kn_k^D=g_j.
$$

Darcy advection is added to this flux before matching paths. For one gas,
molecular friction cancels exactly: $n^D=Kg$. A zero wall coefficient closes
that species' wall channel without an artificial composition floor.

The capillary choice is $d_p=\sqrt{32\tau\kappa_{eff}/\phi}$. It uses the
effective permeability; averaging component square roots defines a different
closure and is not performed. The exact square root is supported where its
state path stays positive or the channel is identically closed. A path crossing
zero requires an explicit diameter or `capillary_regularization=ε>0`, which
uses $\sqrt{32\tau/\phi}\,\kappa/\sqrt{\kappa+\varepsilon}$. This preserves zero
permeability but changes the constitutive law near $\varepsilon$; sensitivity
to that physical regularization must be established for a chosen problem.

## 9.5 Conditioned constitutive evaluation

For one consistent face pressure state, define

$$
u_j=tPK_j,\quad b_j=(1+u_j)^{-1},\quad a_j=K_jb_j,\quad
c_j=u_jb_jx_j,\quad d=\sum_jx_jb_j.
$$

Then $n_t=(\sum_ja_jg_j)/d$ and $n_j^D=a_jg_j+c_jn_t$.
The positive sum for $d$ avoids subtracting nearly equal numbers at strong
coupling. At small $u$, the equivalent $d=1-\sum_jta_jp_j$ is well conditioned
and supplies the finite vacuum continuation. At strong isobaric coupling the
numerator is evaluated as $(\sum_jg_j-\sum_jb_jg_j)/(tP)$, with compensated
summation and a total pressure difference formed from the original endpoint
pressures. This avoids cancelling the common diagonal contribution.

The half-path evaluator uses $\bar P=(P_L+P_R)/2$ throughout. It selects a
collective donor from the sign of $S=\sum_ja_jg_j$ and sets
$p_j^f=\bar P p_{j,up}/P_{up}$. Normalization has no composition floor.
An empty donor on a nonnegative state has zero collective drive. Two-sided
vacuum uses the finite wall-law continuation and its first state derivative.
Darcy density uses its own velocity donor and that donor's actual temperature.
Branch choices inspect primal values; the zero-drive choice is deterministic.

## 9.6 Composite faces and qualified local roots

Each interior face introduces nonnegative interface partial pressures $p_j^*$
and enforces equality of **total** molar species flux through its two half
cells. Each half retains its own coefficients and distance. The interface
thermal state is $T^*=(d_RT_L+d_LT_R)/(d_L+d_R)$; it also supplies carried
enthalpy. For an isolated linear wall channel,

$$
n_j=-\frac{p_{j,R}-p_{j,L}}{R(d_LT_L/K_{j,L}+d_RT_R/K_{j,R})}.
$$

For the full law the mechanisms are combined before continuity is imposed.
In particular, a wall-only layer followed by a Darcy-only layer remains
conducting. The lumped and pure-gas equations have a stable analytic quadratic
reduction. Other cases use a scaled Newton solve, a nonnegative domain, a
merit-function line search and a relative trust bound. Inactive equations are
replaced by explicit choices of their unobservable interface pressures.

Newton steps and implicit derivatives use the same row and column
equilibration. It rescales equations and unknowns without altering pressures,
compositions or fluxes. Conditioning is assessed on the equilibrated Jacobian,
so a small trace-species flux does not by itself imply a singular system.
Acceptance requires residual convergence, a finite nonsingular equilibrated
Jacobian, a resolved correction and matching half-path fluxes. Failure raises
`TransportFailure` with a reason; an unconverged flux is never assembled.
Very strong nearly isobaric coupling can require pressure corrections below
Float64 resolution and is diagnosed rather than accepted. Qualification is
local; no global uniqueness theorem is claimed for arbitrary constitutive data.
First derivatives use the implicit function theorem at the converged root.
They require a regular root and a fixed branch; derivatives at a degenerate
channel or a donor switch require the stated limiting branch. Higher implicit
derivatives are not provided by this first-derivative construction.

## 9.7 Dissipation, positivity and conservation

For positive partial pressures on an isothermal fixed geometry and pore
volume, write $g_j=-c\Delta p_j$ with $c>0$. The diagonal wall contribution is
dissipative species by species. For a left collective donor,

$$
\sum_ja_jp_{j,L}\log(p_{j,R}/p_{j,L})\leq\sum_ja_j\Delta p_j=-S/c\leq0.
$$

For a right donor the reverse logarithmic inequality gives the corresponding
nonpositive collective contribution after multiplying by its signed flux.
The Darcy donor obeys the same argument. Thus
$\sum_jn_j\Delta\log p_j\leq0$. At a positive, flux-matched interface the two
log-pressure jumps telescope. Vacuum is handled by the limiting extensions;
the logarithmic proof itself requires positive pressures. This argument does
not establish nonisothermal entropy stability or time-discrete positivity.

At an absent species endpoint the diagonal term points inward; donor terms
are inward or zero. This gives semi-discrete positivity in the stated domain.
Conservation is separate: a single face vector is added with opposite signs
to its two cells, giving species-wise cancellation independent of the
constitutive model. ALE fluxes and deactivation transfers use their own
conservative assembly.

## 9.8 Reservoirs, films and pressure boundaries

`GasReservoir(P=...,T=...,mole_fractions=(species=...,))` specifies complete
external thermodynamic data. Fractions are named, nonnegative and normalized;
unknown nonzero species are rejected. Pressure, temperature and composition
may depend on time. An explicit vacuum has `P=0` and may omit composition.
`AmbientGasBC(reservoir)` connects this state to the pore half path.

An optional `MassTransferFilm(g)` has common-diffusivity isobaric conductance
$g=\rho D/\delta$ [kg/(m² s)]. With outward mass flux $m=\sum_jJ_j$,

$$
J_j=Y_{j,s}m+g\,\beta(m/g)(Y_{j,s}-Y_{j,\infty}),\qquad
\beta(B)=B/(e^B-1),\quad\beta(0)=1.
$$

At $g=0$ only donor bulk flow remains; at infinite conductance the surface
composition is fixed to the reservoir. A finite film solves $N_g-1$ surface
fractions, with the last from normalization. Pressure, composition and surface
temperature are solved together when coupled. A completely closed pore channel
with zero film conductance has unobservable surface composition, which is
explicitly pinned to the reservoir state.

`NeumannPressureBC(gradient=...,T=...,mole_fractions=...)` specifies the outward
pressure gradient and complete species/energy closure.
`ConvectivePressureBC(reservoir;resistance=...,flow=:darcy_velocity)` specifies
$P_s-P_\infty=R_h f$. `flow` may also be `:mass` or `:molar`; resistance units
follow the selected rate per area. It may carry the same optional mass film.
`ImpermeableBC()` closes all species. `MassFluxBC(fluxes=(...),
inflow_temperature=(...))` prescribes named outward species fluxes, including
negative inflow with required incoming temperature. Mass conditions do not add
with `+`.

## 9.9 Boundary energy and thermal films

All boundary energy is assembled at the direct-reservoir control surface.
Outward species carry $h_j(T_s)$ and inward species carry their specified
$h_j(T_{in})$. Define $H_s=\sum_jJ_jh_j(T_s)$ and the signed donor sum $H_r$.
For an LTE wall heat film with heat $q_s$ positive inward,

$$
q_r=q_s+H_r-H_s,\qquad E_{out}=H_r-q_r=H_s-q_s.
$$

This conversion accounts for incoming heating exactly once. An adiabatic,
prescribed-flux or prescribed-temperature boundary without a thermal film
uses its prescribed heat and donor enthalpy directly. A convective film and
its gas reservoir must share the incoming stream temperature; radiative
surroundings may differ.

For constant conductivity and species heat capacities the Couette film gives
$q_s=h\beta(C/h)(T_\infty-T_s)$ with the **signed**
$C=\sum_jJ_jc_{p,j}$. The configurable $2\lambda$ factor has its classical
value at $\lambda=1/2$. Temperature-dependent properties are evaluated at the
surface as a frozen-property approximation. `PoolingFilm` uses the solved
surface composition when supplied; empirical pooling/property-ratio closures
support outward products or pure freestream suction and diagnose unsupported
mixed counterflow. A custom `FunctionBlowing` receives signed species fluxes
and the surface state.

## 9.10 Admissibility and pressure-aware integration

Physical pore volume must remain positive outside the trial guard:
$\phi\geq2\,10^{-12}$ plus a resolved uncertainty margin from condensed masses,
temperature and geometry. This holds initially, on accepted steps and after
deactivation. A divisor floor does not make a zero-pore physical state valid.
Negative trial masses and distorted trial geometry have finite extensions so
Newton can recover; they are not accepted physical states.

`pressure_tolerance` [Pa] sets an absolute local pressure-error budget shared
among gases. Initial gas-row scales are proportional to
$M_j\phi V p_{tol,j}/(RT)$. A proposed-state error norm also propagates thermal,
condensed-mass and geometry errors through the current pore EOS. The gas
positive-part scales are $\delta_{j}=M_j\phi V p_{tol}/(N_gRT)$ at the
current state. For the cubic positive part, the pressure bias is at most
$0.26p_{tol}$ summed over all gases: the maximum of
$s(1-3s^2+2s^3)$ on $[0,1]$ is below 0.26. The smoothing transition and gas-row
error scales therefore remain commensurate during trace-gas production.
Exactly absent gas still contributes zero pressure. A prescribed reservoir
pressure is used exactly; an unknown hydraulic surface pressure uses the
corresponding total-pressure trial extension. Residual evaluation does not
mutate solver tolerances. Local error control is not a global error
bound; tight-reference comparisons remain necessary as pores evolve.

The 5.8 mm PMMA cone deck uses nominal φ=0.01. Its qualification varies
φ over 0.01–0.03, uses 24/48/72 cells and pressure budgets 1e-4/1e-5 Pa
with reltol 1e-9 and temperature abstol 1e-11 K, in both energy formulations.
Temperature and escaped-product differences contract with mesh refinement;
these bounds qualify this deck and its coefficients, not arbitrary dense melts.
The calorically passive-gas burnout fixture separately covers φ=0.003–0.03
on 15/30/60 PMMA cells with pressure budgets 1e-3/1e-4 Pa. The heterogeneous
burnout fixture uses φ=0.1 and qualifies φ=0.1–0.15 on the same mesh and
pressure-budget ladder. Its moving-mesh redistribution produces a finite
matrix-occupancy error: maximum occupancy is 1.091/1.036/1.016 at
15/30/60 cells. The actual pore-volume check therefore matters even when
nominal porosity is positive. Thermal independence of its calorically
passive gas and a separate species ledger qualify that fixture.

Nominal porosity in a dense melt is material data, separate from the trial
floor. Its consequences must be checked jointly with mesh and pressure
tolerance, including release, temperature, pressure, conservation and solver
work. The separate pressure-work/internal-energy-storage formulation is not
part of these transport equations; both supported energy states retain the
package's sensible-enthalpy convention.

## 9.11 Diagnostics and verification

`get_porosity`, `get_pressure`, `get_partial_pressures`,
`get_intrinsic_gas_densities` and `get_gas_composition` read saved central pore
states. Gas arrays follow `sol.gas.species`. `get_boundary_exchange(sol,tag)`
returns current-area averaged outward species fluxes, surface pressures,
Darcy velocity, carried enthalpy, heat, total energy and local-solve diagnostics.
Dead cells or absent live boundaries are NaN. Velocity means Darcy velocity,
not total mass or molar flow. `gas_species_audit(sol,recorder)` separates
reaction production, escape, ingress and deactivation venting by species.

Independent verification includes high-precision full friction matrices,
original-input vacuum/strong-coupling references, donor inequalities,
interface continuity and implicit derivatives. V-GT-01 tests Gaussian variance
$\sigma^2=\sigma_0^2+2\Lambda t/\phi$; V-GT-02 tests the Darcy pressure-squared
profile; V-GT-03 tests integrated single/two-layer pressure potentials and
mechanism switches; V-GT-04 tests Graham's binary flux relation and its
logarithmic continuum integral. Manufactured cases exercise the coupled
pressure potential and carried energy on fixed and moving meshes.
