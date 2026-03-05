# Global Sensitivity Analysis of Douglas Fir Pyrolysis via Polynomial Chaos Expansion

## Abstract

A variance-based global sensitivity analysis is presented for a one-dimensional
reactive-transport pyrolysis model of Douglas Fir exposed to 50 kW/m² radiant
heat flux in a cone calorimeter. Seventeen uncertain thermophysical input
parameters, grouped into eight statistically independent measurement sources,
are propagated through the model using Sobol quasi-random sampling ($N = 2048$)
and Polynomial Chaos Expansion (PCE) surrogates. Group Sobol indices decompose
the output variance of six scalar quantities of interest into contributions from
each measurement source. The analysis reveals that only two sources—the heats of
pyrolysis reaction (five parameters) and the DSC calibration bias (one
parameter)—account for more than 98% of the output variance across all
quantities. The remaining six groups (HFM calibration, water calibration, virgin
specific heat, char specific heat, and heat of sorption) each contribute less
than 2%. Bootstrap confidence intervals on the Sobol indices are extremely narrow
(half-widths < 0.1 percentage points), and convergence diagnostics confirm index
stability. These results provide a rigorous, GUM-traceable ranking of which
experimental measurements most merit refinement to reduce prediction uncertainty.

---

## 1. Introduction

Pyrolysis models for wood and wood-based materials require numerous
thermophysical input parameters obtained from bench-scale measurements. Each
measurement carries uncertainties from instrument calibration (Type B, systematic)
and regression fitting (Type A, statistical), as categorized by the Guide to the
Expression of Uncertainty in Measurement (GUM; JCGM 100:2008). Understanding
which measurement sources dominate the prediction uncertainty is essential for
directing experimental effort and for providing defensible uncertainty bounds on
model outputs.

This report details a global sensitivity analysis (GSA) of a Douglas Fir cone
calorimeter simulation. Unlike local (one-at-a-time) sensitivity analyses, the
GSA explores the full joint input space simultaneously, capturing nonlinear
effects and parameter interactions. Variance-based Sobol indices are computed via
Polynomial Chaos Expansion, providing an exact analytical decomposition of the
output variance.

## 2. Model Description

The forward model is a one-dimensional reactive-transport pyrolysis simulation
implemented in the ThermaKin framework (Julia). The simulation represents a
12 mm thick Douglas Fir slab exposed to 50 kW/m² external radiant heat flux with
convective cooling ($h = 10$ W/(m²·K)) at the heated surface and an adiabatic
back face. The material undergoes five sequential pyrolysis reactions with
moisture evaporation, producing gaseous products that transport through the
porous char layer via Darcy flow and Fickian diffusion.

**Simulation parameters (fixed):**

| Parameter | Value |
|---|---|
| Sample thickness | 12.0 mm |
| Spatial cells | 240 |
| Initial temperature | 300 K |
| External heat flux | 50 kW/m² |
| Convective coefficient | 10 W/(m²·K) |
| Radiation model | Surface absorption |
| Moisture content (dry basis) | 6% |
| Simulation duration | 600 s |
| ODE tolerances | abstol = $10^{-9}$, reltol = $10^{-7}$ |

**Quantities of interest (QoIs):**

| QoI | Symbol | Units | Description |
|---|---|---|---|
| Peak mass loss rate | $\dot{m}''_\text{peak}$ | kg/(m²·s) | Maximum mass loss rate |
| Time to peak MLR | $t_\text{peak}$ | s | Time at which peak MLR occurs |
| Total mass loss | $\Delta m$ | % | Final mass loss as percentage of initial |
| Peak surface temperature | $T_\text{s,max}$ | K | Maximum temperature at the heated face |
| Peak back face temperature | $T_\text{b,max}$ | K | Maximum temperature at the insulated face |
| Final thickness | $\ell_f$ | m | Residual sample thickness after pyrolysis |

Time to ignition was excluded from the Sobol analysis because it is quantized by
the simulation output time step, producing a discrete distribution unsuitable for
continuous PCE fitting.

## 3. Uncertain Parameters and Covariance Structure

### 3.1 Parameter vector

The 17 uncertain parameters are collected into a vector
$\boldsymbol{\theta} \in \mathbb{R}^{17}$, sampled from eight independent
multivariate normal groups:

$$\boldsymbol{\theta} = (\varepsilon_\text{HFM},\; \varepsilon_\text{STA},\; \varepsilon_\text{water},\; a_1, b_1,\; a_3, b_3,\; a_\text{char}, b_\text{char},\; a_s, b_s, c_s,\; \Delta H_1, \ldots, \Delta H_5)$$

### 3.2 Group definitions

The eight groups correspond to distinct measurement sources—instrument
calibrations and regression fits—such that parameters within a group may be
correlated (sharing a common instrument or fitting procedure) but parameters
across groups are statistically independent.

#### Group 1: HFM calibration bias ($\varepsilon_\text{HFM}$) — 1 parameter

The heat flow meter heat flow accuracy is specified as $\pm 1\%$ by the
manufacturer. Assuming a rectangular distribution and converting to a standard
uncertainty:

$$\varepsilon_\text{HFM} \sim \mathcal{N}(0,\; (0.01/\sqrt{3})^2) = \mathcal{N}(0,\; 0.005774^2)$$

This multiplicative bias scales all HFM-derived specific heat values (virgin
$c_p$ segment 1, $T \leq 323$ K) via the factor $(1 + \varepsilon_\text{HFM})$.

#### Group 2: DSC calibration bias ($\varepsilon_\text{STA}$) — 1 parameter

The DSC enthalpy calibration was performed against seven reference standards
(after exclusion of RbNO₃), yielding:

$$\varepsilon_\text{STA} \sim \mathcal{N}(0,\; 0.02613^2)$$

This 2.61% relative uncertainty (1$\sigma$) is the dominant component of the
combined DSC systematic uncertainty budget (temperature calibration contributes
only 4.6% of the total Type B variance; mass measurement contributes 0.03%).
Applied as $(1 + \varepsilon_\text{STA})$ to all DSC-derived quantities: virgin
$c_p$ segment 3 ($T > 425$ K), char $c_p$, and all heats of reaction.

#### Group 3: Water calibration bias ($\varepsilon_\text{water}$) — 1 parameter

The heat of sorption was measured isothermally on the STA using an independent
water-evaporation calibration (not the standard multi-point DSC calibration),
giving a separate Type B uncertainty:

$$\varepsilon_\text{water} \sim \mathcal{N}(0,\; 0.01732^2)$$

from 22 replicate runs (11 Douglas Fir + 11 red oak), 21 degrees of freedom.
Applied as $(1 + \varepsilon_\text{water})$ to the heat of sorption baseline and
amplitude parameters.

#### Group 4: Virgin $c_p$, low-$T$ (HFM regression) — 2 parameters

Linear regression of HFM-measured specific heat for $T \leq 323$ K:
$c_p(T) = a_1 + b_1 T$ in J/(kg·K).

$$\begin{pmatrix} a_1 \\ b_1 \end{pmatrix} \sim \mathcal{N}\!\left(\begin{pmatrix} -53.05 \\ 4.245 \end{pmatrix},\; \begin{pmatrix} 4.761 \times 10^4 & -153.5 \\ -153.5 & 0.4953 \end{pmatrix}\right)$$

The negative off-diagonal reflects the typical intercept–slope anti-correlation
in linear regression.

#### Group 5: Virgin $c_p$, high-$T$ (DSC regression) — 2 parameters

Linear regression of DSC-measured specific heat for $T > 425$ K:
$c_p(T) = a_3 + b_3 T$ in J/(kg·K).

$$\begin{pmatrix} a_3 \\ b_3 \end{pmatrix} \sim \mathcal{N}\!\left(\begin{pmatrix} -2131.3 \\ 9.119 \end{pmatrix},\; \begin{pmatrix} 3.181 \times 10^3 & -7.200 \\ -7.200 & 1.630 \times 10^{-2} \end{pmatrix}\right)$$

Segment 2 (323–425 K) is linearly interpolated between the endpoints of segments
1 and 3, so it inherits uncertainties from both groups through the boundary
values.

#### Group 6: Char $c_p$ (independent STA) — 2 parameters

Linear regression of independently measured char specific heat:
$c_{p,\text{char}}(T) = a_\text{char} + b_\text{char} T$ in J/(kg·K).

$$\begin{pmatrix} a_\text{char} \\ b_\text{char} \end{pmatrix} \sim \mathcal{N}\!\left(\begin{pmatrix} 347.2 \\ 1.784 \end{pmatrix},\; \begin{pmatrix} 1.215 \times 10^1 & -2.005 \times 10^{-2} \\ -2.005 \times 10^{-2} & 3.343 \times 10^{-5} \end{pmatrix}\right)$$

#### Group 7: Heat of sorption (WLS exponential fit) — 3 parameters

The moisture heat of sorption was fit to a weighted least-squares exponential
model: $h_s(\text{MC}) = a_s + b_s \exp(-c_s \cdot \text{MC})$ in kJ/kg, where
MC is moisture content in percent (dry basis).

$$\begin{pmatrix} a_s \\ b_s \\ c_s \end{pmatrix} \sim \mathcal{N}\!\left(\begin{pmatrix} 2280 \\ 787.0 \\ 0.1144 \end{pmatrix},\; \begin{pmatrix} 174.0 & 22.03 & 0.09442 \\ 22.03 & 941.5 & 0.1843 \\ 0.09442 & 0.1843 & 9.715 \times 10^{-5} \end{pmatrix}\right)$$

The Type B calibration uncertainty for heat of sorption is carried by
$\varepsilon_\text{water}$ (Group 3), applied multiplicatively.

#### Group 8: Heats of reaction ($\Delta H_1$–$\Delta H_5$) — 5 parameters

The five pyrolysis reaction enthalpies (J/kg) were determined by a two-stage
weighted least-squares optimization. The reported covariance is the *conditional*
covariance $\boldsymbol{\Sigma}_{\Delta H | c_{p,\text{virgin}}}$, computed via
the analytical delta method (implicit function theorem + chain rule through the
pseudo-char $c_p$ re-optimization):

$$\boldsymbol{\Sigma}_{\Delta H | c_p} = \boldsymbol{\Sigma}_\text{WLS} + \mathbf{J}_\text{char}\, \boldsymbol{\Sigma}_\text{char}\, \mathbf{J}_\text{char}^\top + \text{cross terms}$$

where $\boldsymbol{\Sigma}_\text{WLS}$ is the Stage 2 regression covariance and
$\mathbf{J}_\text{char}$ propagates Stage 1 pseudo-char $c_p$ estimation
uncertainty.

$$\boldsymbol{\mu}_{\Delta H} = (-19153,\; -5574,\; 266700,\; -219639,\; -12382)^\top \text{ J/kg}$$

$$\boldsymbol{\Sigma}_{\Delta H | c_p} = \begin{pmatrix}
3.423 & 8.046 & 11.37 & 50.16 & 12.33 \\
8.046 & 20.80 & 28.70 & 127.5 & 31.29 \\
11.37 & 28.70 & 41.15 & 178.3 & 45.59 \\
50.16 & 127.5 & 178.3 & 805.9 & 176.3 \\
12.33 & 31.29 & 45.59 & 176.3 & 313.7
\end{pmatrix} \times 10^6 \text{ J}^2/\text{kg}^2$$

**Why the conditional covariance is used:** The full bootstrap covariance of
$\Delta H$ includes the effect of varying virgin $c_p$ within the optimization.
Using it would double-count the virgin $c_p$ uncertainty, since Group 5 already
perturbs $(a_3, b_3)$ directly. Instead, the virgin $c_p \to \Delta H$ coupling
is handled explicitly at runtime via a Jacobian correction:

$$\Delta H_i^{\text{adj}} = \Delta H_i + \sum_{j \in \{a_3, b_3\}} J_{ij}\, \delta_j, \qquad \mathbf{J} = \frac{\partial \Delta H}{\partial (a_3, b_3)}$$

where $\delta_j = \theta_j - \theta_j^{\text{nom}}$ is the deviation of the
virgin $c_p$ regression coefficient from its nominal value, and $\mathbf{J}$ is
the $5 \times 2$ Jacobian computed analytically via the implicit function theorem.
After the Type A correction, the DSC calibration bias is applied multiplicatively:
$\Delta H_i^{\text{final}} = (1 + \varepsilon_\text{STA}) \cdot \Delta H_i^{\text{adj}}$.

### 3.3 Summary of included and omitted uncertainty sources

| Source | Type | Status | Justification |
|---|---|---|---|
| HFM heat flow accuracy | B | **Included** | Dominant HFM systematic |
| DSC multi-point enthalpy calibration | B | **Included** | Dominant DSC systematic |
| Water-evaporation calibration | B | **Included** | Independent HoS calibration |
| Virgin $c_p$ regression (seg 1, seg 3) | A | **Included** | Full covariance |
| Char $c_p$ regression | A | **Included** | Full covariance |
| Heat of sorption WLS fit | A | **Included** | Full covariance |
| Heats of reaction (conditional) | A | **Included** | Conditional covariance + Jacobian |
| STA temperature calibration on $c_p$ | B | Omitted | 4.6% of STA Type B variance |
| STA temperature calibration on $\Delta H$ | B | Omitted | Second-order: $O(\delta T^2 / \Delta T_\text{peak}^2)$ |
| STA mass measurement | B | Omitted | $u_\text{rel} = 0.03\%$; negligible vs $\varepsilon_\text{STA}$ |
| HFM temperature (Type K TC) on $c_p$ | B | Omitted | $\delta c_p \approx 0.004$ J/(g·K) vs SE = 0.218 |
| HFM density systematic | B | Omitted | $u_\text{rel} = 0.06\%$; negligible |
| Kinetic parameters ($A$, $E$) | — | Omitted | Outside scope of thermochemical study |

## 4. Sampling Strategy

### 4.1 Sobol quasi-random sequences

Parameter samples are generated using Sobol low-discrepancy sequences in
$[0,1]^{17}$, which provide more uniform coverage of the input space than
pseudo-random sampling, reducing the integration error for a given sample size.
$N = 2048 = 2^{11}$ samples are used, chosen as a power of 2 to preserve the
balanced properties of the Sobol sequence.

### 4.2 Transformation to correlated samples

For each group $g$ with mean $\boldsymbol{\mu}_g$ and covariance
$\boldsymbol{\Sigma}_g = \mathbf{L}_g \mathbf{L}_g^\top$ (Cholesky
factorization):

1. Draw Sobol points $\mathbf{u}_g \in [0,1]^{k_g}$ for the group's $k_g$ dimensions
2. Transform to standard normals: $z_{g,j} = \Phi^{-1}(u_{g,j})$ for each marginal
3. Induce correlations: $\mathbf{x}_g = \boldsymbol{\mu}_g + \mathbf{L}_g \mathbf{z}_g$

This Nataf-like transformation preserves the low-discrepancy property in each
marginal while reproducing the target covariance structure.

### 4.3 Simulation execution

All 2048 parameter vectors were evaluated by the pyrolysis solver in parallel
using Julia's `Distributed` module. All runs completed successfully (2048/2048),
requiring no rejection or imputation. Scalar QoIs were extracted from each
solution trajectory.

## 5. Polynomial Chaos Expansion Methodology

### 5.1 Transformation to independent standard normal space

Because the PCE basis requires independent inputs, the correlated physical
parameters are first transformed to independent standard normals using the
inverse of the sampling transformation. For each group $g$:

$$\mathbf{z}_g = \mathbf{L}_g^{-1} (\mathbf{x}_g - \boldsymbol{\mu}_g)$$

**Verification:** After transformation, all 17 standard normal variables have
empirical mean $|\mu| < 10^{-4}$ and standard deviation $\sigma = 0.9999$. The
maximum absolute cross-group correlation is $|\rho_\text{max}| = 0.0065$,
confirming effective independence (expected $O(1/\sqrt{N}) \approx 0.022$ for
$N = 2048$).

### 5.2 Hermite polynomial basis

For standard normal inputs, the natural orthogonal polynomial basis is the
probabilist's Hermite polynomials $\text{He}_k(z)$, satisfying:

$$\mathbb{E}[\text{He}_m(Z) \cdot \text{He}_n(Z)] = n!\, \delta_{mn}, \qquad Z \sim \mathcal{N}(0,1)$$

The multivariate basis functions are tensor products of univariate Hermite
polynomials, indexed by a multi-index
$\boldsymbol{\alpha} = (\alpha_1, \ldots, \alpha_{17})$:

$$\Psi_{\boldsymbol{\alpha}}(\mathbf{z}) = \prod_{d=1}^{17} \text{He}_{\alpha_d}(z_d)$$

These are then orthonormalized by dividing by $\sqrt{\prod_d \alpha_d!}$ to obtain
$\mathbb{E}[\Psi_{\boldsymbol{\alpha}}^2] = 1$.

### 5.3 Truncation and basis size

The total-degree truncation retains all basis functions with
$|\boldsymbol{\alpha}| = \sum_d \alpha_d \leq p$, giving:

$$P = \binom{17 + p}{p}$$

| Total degree $p$ | Basis size $P$ | Oversampling ratio $N/P$ |
|---|---|---|
| 2 | 171 | 12.0 |
| 3 | 1140 | 1.8 |
| 4 | 5985 | 0.34 (insufficient) |

### 5.4 Adaptive degree selection

Degrees 2 through 4 are evaluated. For each candidate degree, a PCE surrogate is
fit to all six QoIs, and the average leave-one-out cross-validation score
$\overline{Q^2}$ is computed. The degree maximizing $\overline{Q^2}$ is selected.

**Fitting method.** When the oversampling ratio exceeds the sparse threshold
($N/P > 3$), ordinary least squares (OLS) is used. When it falls below this
threshold (as for degree 3), a sparse LASSO fit is employed (Section 5.5).

**Results of degree selection:**

| Degree | Basis size | Oversampling | Fit method | $\overline{Q^2}$ | $\min Q^2$ | Selected |
|---|---|---|---|---|---|---|
| 2 | 171 | 12.0× | OLS | 0.9999 | 0.9996 | **Yes** |
| 3 | 1140 | 1.8× | LASSO | 0.9992 | 0.9990 | No |
| 4 | 5985 | 0.34× | — | — | — | Skipped |

Degree 2 was selected as it achieved the highest cross-validation score. The
marginal improvement from degree 3 was negligible ($\Delta Q^2 = -0.0007$),
indicating that the model response is well-approximated by linear, quadratic,
and bilinear terms within the GUM uncertainty bounds.

### 5.5 Sparse LASSO regression (for degree $\geq$ 3)

When the oversampling ratio is low, the $\ell_1$-regularized LASSO estimator is
used to select a sparse subset of active basis functions:

$$\hat{\mathbf{c}} = \arg\min_{\mathbf{c}} \left\{ \frac{1}{2N} \|\mathbf{y} - \boldsymbol{\Phi} \mathbf{c}\|_2^2 + \lambda \sum_{j=2}^{P} |c_j| \right\}$$

where the constant term ($j=1$) is unpenalized. The regularization parameter
$\lambda$ is selected by 5-fold cross-validation over a logarithmic grid of 50
values from $\lambda_\text{max}$ to $10^{-4} \lambda_\text{max}$, where
$\lambda_\text{max} = \max_j |\boldsymbol{\Phi}_j^\top \mathbf{y}| / N$ is the
smallest $\lambda$ producing an all-zero solution.

The LASSO is solved by coordinate descent with residual tracking (each coordinate
update costs $O(N)$ rather than $O(NP)$) and warm-starting along the
$\lambda$-path. For degree 3, the LASSO selected 14–30 nonzero coefficients out
of 1139 (1.2–2.6% sparsity), confirming that most higher-order interactions are
negligible.

### 5.6 Leave-one-out cross-validation ($Q^2$)

For the OLS fit, $Q^2$ is computed analytically without refitting, using the hat
matrix $\mathbf{H} = \boldsymbol{\Phi}(\boldsymbol{\Phi}^\top\boldsymbol{\Phi})^{-1}\boldsymbol{\Phi}^\top$:

$$Q^2 = 1 - \frac{\sum_{i=1}^{N} \left[\frac{e_i}{1 - h_{ii}}\right]^2}{\sum_{i=1}^{N} (y_i - \bar{y})^2}$$

where $e_i = y_i - \hat{y}_i$ are the residuals and $h_{ii}$ are the diagonal
elements of $\mathbf{H}$.

### 5.7 Numerical conditioning verification

- **Condition number:** $\kappa(\boldsymbol{\Phi}^\top\boldsymbol{\Phi}) = 3.96$,
  indicating excellent numerical stability.
- **Basis orthonormality:** Empirical verification shows
  $\max_j |\mathbb{E}[\Psi_j^2] - 1| = 0.135$, within acceptable tolerance
  for $N = 2048$ samples (expected sampling fluctuation is $O(1/\sqrt{N}) \approx 0.02$
  for low-order terms; higher-order terms show larger fluctuations due to the
  heavier tails of Hermite polynomials).

## 6. Sobol Index Computation

### 6.1 Variance decomposition from PCE coefficients

For an orthonormal PCE:

$$Y \approx \sum_{\boldsymbol{\alpha}} c_{\boldsymbol{\alpha}}\, \Psi_{\boldsymbol{\alpha}}(\mathbf{Z})$$

the total variance is:

$$D = \text{Var}(Y) = \sum_{\boldsymbol{\alpha} \neq \mathbf{0}} c_{\boldsymbol{\alpha}}^2$$

### 6.2 Group Sobol indices

Let $\mathcal{G}_g$ denote the set of parameter indices belonging to group $g$.
A basis function $\Psi_{\boldsymbol{\alpha}}$ "involves" group $g$ if
$\alpha_d > 0$ for some $d \in \mathcal{G}_g$.

**First-order group index** (fraction of variance due to group $g$ alone):

$$S_g = \frac{1}{D} \sum_{\substack{\boldsymbol{\alpha}:\; \alpha_d > 0 \text{ only for } d \in \mathcal{G}_g}} c_{\boldsymbol{\alpha}}^2$$

**Total group index** (fraction of variance due to group $g$ including all
interactions with other groups):

$$S_{T,g} = \frac{1}{D} \sum_{\substack{\boldsymbol{\alpha}:\; \exists d \in \mathcal{G}_g \text{ with } \alpha_d > 0}} c_{\boldsymbol{\alpha}}^2$$

By construction, $\sum_g S_g \leq 1$ (with equality when there are no
inter-group interactions) and $\sum_g S_{T,g} \geq 1$.

## 7. Results

### 7.1 PCE surrogate quality

| QoI | $R^2$ | $Q^2_\text{LOO}$ | Basis size | Fit method |
|---|---|---|---|---|
| Peak MLR | 1.0000 | 1.0000 | 171 | OLS |
| Time to peak MLR | 0.9999 | 0.9999 | 171 | OLS |
| Total mass loss | 1.0000 | 1.0000 | 171 | OLS |
| Peak surface temp. | 0.9999 | 0.9999 | 171 | OLS |
| Peak back face temp. | 0.9997 | 0.9996 | 171 | OLS |
| Final thickness | 1.0000 | 1.0000 | 171 | OLS |

All QoIs achieve $Q^2 > 0.999$, indicating that the degree-2 PCE is an
effectively exact surrogate over the sampled input domain. The close agreement
$R^2 \approx Q^2$ confirms no overfitting.

### 7.2 Group Sobol indices

The following table presents the first-order ($S_g$) and total ($S_{T,g}$)
Sobol indices for each measurement source group and QoI. Values below 0.001 are
shown as "< 0.1%". The sum $\sum_g S_g$ and the interaction remainder
$1 - \sum_g S_g$ are reported in the final rows.

| Group | Peak MLR ||| Time to peak MLR ||| Total mass loss ||| Peak surface temp. ||| Peak back temp. ||| Final thickness |||
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| | $S_g$ | $S_{T,g}$ | | $S_g$ | $S_{T,g}$ | | $S_g$ | $S_{T,g}$ | | $S_g$ | $S_{T,g}$ | | $S_g$ | $S_{T,g}$ | | $S_g$ | $S_{T,g}$ | |
| HFM calibration | < 0.1% | < 0.1% | | 0.05% | 0.05% | | 0.03% | 0.03% | | < 0.1% | < 0.1% | | < 0.1% | < 0.1% | | 0.03% | 0.03% | |
| **DSC calibration** | **38.30%** | **38.33%** | | **80.48%** | **80.49%** | | **56.32%** | **56.42%** | | 0.09% | 0.11% | | 0.05% | 0.09% | | **56.30%** | **56.42%** | |
| Water calibration | < 0.1% | < 0.1% | | 1.51% | 1.51% | | 0.72% | 0.72% | | < 0.1% | < 0.1% | | < 0.1% | < 0.1% | | 0.72% | 0.72% | |
| Virgin $c_p$, low-$T$ | < 0.1% | < 0.1% | | 0.08% | 0.08% | | 0.04% | 0.04% | | < 0.1% | < 0.1% | | < 0.1% | < 0.1% | | 0.04% | 0.04% | |
| Virgin $c_p$, high-$T$ | < 0.1% | < 0.1% | | 0.29% | 0.29% | | 0.17% | 0.17% | | < 0.1% | < 0.1% | | < 0.1% | < 0.1% | | 0.17% | 0.17% | |
| Char $c_p$ | < 0.1% | < 0.1% | | < 0.1% | < 0.1% | | < 0.1% | < 0.1% | | < 0.1% | < 0.1% | | < 0.1% | < 0.1% | | < 0.1% | < 0.1% | |
| Heat of sorption | < 0.1% | < 0.1% | | 0.14% | 0.14% | | 0.07% | 0.07% | | < 0.1% | < 0.1% | | < 0.1% | < 0.1% | | 0.07% | 0.07% | |
| **Heats of reaction** | **61.67%** | **61.70%** | | **17.44%** | **17.45%** | | **42.53%** | **42.65%** | | **99.89%** | **99.90%** | | **99.91%** | **99.95%** | | **42.55%** | **42.67%** | |
| $\sum S_g$ | 99.97% | | | 99.99% | | | 99.88% | | | 99.98% | | | 99.96% | | | 99.88% | | |
| Remainder (interactions) | 0.03% | | | 0.01% | | | 0.12% | | | 0.02% | | | 0.04% | | | 0.12% | | |

### 7.3 Bootstrap confidence intervals

Nonparametric bootstrap confidence intervals (1000 replicates, resampling
parameter–QoI pairs with replacement) for the dominant Sobol indices:

| QoI | Group | $S_g$ | 95% CI |
|---|---|---|---|
| Peak MLR | Heats of reaction | 0.6167 | [0.6166, 0.6168] |
| Peak MLR | DSC calibration | 0.3830 | [0.3829, 0.3831] |
| Time to peak MLR | DSC calibration | 0.8048 | [0.8045, 0.8050] |
| Time to peak MLR | Heats of reaction | 0.1744 | [0.1741, 0.1746] |
| Total mass loss | DSC calibration | 0.5632 | [0.5631, 0.5632] |
| Total mass loss | Heats of reaction | 0.4253 | [0.4253, 0.4254] |
| Peak surface temp. | Heats of reaction | 0.9989 | [0.9989, 0.9989] |
| Peak back temp. | Heats of reaction | 0.9991 | [0.9990, 0.9992] |
| Final thickness | DSC calibration | 0.5630 | [0.5630, 0.5631] |
| Final thickness | Heats of reaction | 0.4255 | [0.4254, 0.4256] |

The 95% bootstrap confidence intervals have half-widths of $< 0.02$ percentage
points in all cases, indicating that the sample size is far more than sufficient
for precise index estimation.

### 7.4 Convergence diagnostics

Sobol indices were recomputed on nested subsamples (12.5%, 25%, 50%, 75%, 100%
of the full dataset) to assess convergence:

| QoI | $\max\lvert\Delta S_g\rvert$ (75% → 100%) | Status |
|---|---|---|
| Peak MLR | 0.0000 | ✓ Converged |
| Time to peak MLR | 0.0003 | ✓ Converged |
| Total mass loss | 0.0000 | ✓ Converged |
| Peak surface temp. | 0.0000 | ✓ Converged |
| Peak back temp. | 0.0000 | ✓ Converged |
| Final thickness | 0.0000 | ✓ Converged |

All QoIs show $\max|\Delta S_g| < 0.001$ between the last two subsample sizes,
well below common convergence thresholds (0.01–0.02).

## 8. Discussion

### 8.1 Dominance of two measurement sources

The results reveal a striking concentration of sensitivity. Across all six QoIs,
only two of the eight measurement source groups account for more than 1% of the
output variance:

1. **Heats of reaction** ($\Delta H_1$–$\Delta H_5$): Dominant for temperature
   QoIs (99.9% for both peak surface and back temperatures), and a major
   contributor (43–62%) to mass-related QoIs and peak MLR.

2. **DSC calibration bias** ($\varepsilon_\text{STA}$): Dominant for timing
   (80.5% for time to peak MLR) and primary contributor (56%) to total mass loss
   and final thickness. Moderate contributor (38%) to peak MLR.

The remaining six groups collectively account for less than 2% of output variance
for every QoI.

### 8.2 Physical interpretation

**Temperature QoIs.** The near-exclusive sensitivity to heats of reaction is
physically intuitive. Peak temperatures at both faces are governed by the energy
balance at the pyrolysis front: the exothermic and endothermic reaction enthalpies
directly determine how much thermal energy is absorbed or released during
decomposition, setting the maximum achievable temperature. The DSC calibration
multiplier has a negligible effect on temperatures because it scales both $c_p$
and $\Delta H$ proportionally, and the temperature field is governed primarily
by the $\Delta H / c_p$ ratio, in which the calibration largely cancels.

**Timing QoIs.** The dominance of DSC calibration for time to peak MLR arises
because $\varepsilon_\text{STA}$ scales the apparent specific heat uniformly.
A 2.6% (1$\sigma$) uncertainty in $c_p$ directly shifts the thermal lag—the time
required to heat through the decomposition temperature window—without affecting
the shape of the reaction progress curve. The heats of reaction provide a
secondary timing effect (17.4%) through moderating the net endothermicity of the
pyrolysis front.

**Mass-related QoIs.** Total mass loss and final thickness show a roughly 56/43
split between DSC calibration and heats of reaction. These QoIs depend on both
the thermal transport (how quickly heat penetrates to interior layers, controlled
by $c_p$ and hence $\varepsilon_\text{STA}$) and the reaction energetics (which
reactions run to completion within the 600 s exposure, controlled by $\Delta H$).

### 8.3 Negligible inter-group interactions

The sum $\sum_g S_g \geq 0.9988$ for all QoIs, leaving at most 0.12% attributed
to inter-group interactions. Equivalently, $S_{T,g} \approx S_g$ for all groups
and QoIs. This confirms that the model response is essentially additive in the
group-level uncertainty contributions—the effect of each measurement source on
the output is independent of the values of parameters from other sources.

This finding is consistent with the adequacy of degree-2 PCE (which captures
individual-group quadratic effects but only bilinear inter-group terms) and
explains why a first-order Taylor expansion would have yielded nearly identical
rankings.

### 8.4 Implications for experimental design

1. **Highest-impact improvement:** Reducing the conditional covariance of the
   heats of reaction would decrease prediction uncertainty for all QoIs. This can
   be achieved by (a) additional replicate STA experiments to reduce the WLS
   fitting uncertainty, (b) improved baseline correction to reduce the Stage 1
   pseudo-char $c_p$ uncertainty, or (c) direct measurement of individual reaction
   enthalpies via deconvolution DSC.

2. **Second-highest impact:** Tightening the DSC enthalpy calibration (currently
   $u_\text{rel} = 2.61\%$ from 7 standards) would reduce timing and mass loss
   uncertainty. Adding more calibration standards or using certified reference
   materials with smaller certified uncertainties would help.

3. **Low-priority measurements:** The HFM calibration, water calibration, virgin
   $c_p$ regression, char $c_p$ regression, and heat of sorption collectively
   contribute negligibly. Refining these measurements further would have
   essentially no effect on overall prediction uncertainty.

### 8.5 Robustness checks

The following diagnostic checks support the validity of the results:

- **PCE quality:** $R^2 > 0.999$ and $Q^2 > 0.999$ for all QoIs, with
  $R^2 \approx Q^2$ (no overfitting).
- **Independence verification:** Cross-group correlations after transformation
  are below $|\rho| = 0.007$, validated against the expected $O(N^{-1/2})$
  sampling noise.
- **Basis normalization:** Empirical $\mathbb{E}[\Psi_j^2]$ deviates from unity
  by at most 0.135 (within tolerance for finite-sample Hermite polynomial
  estimation).
- **Numerical conditioning:** $\kappa(\boldsymbol{\Phi}^\top\boldsymbol{\Phi}) = 3.96$,
  confirming a well-conditioned regression problem.
- **Bootstrap CIs:** Sub-percent-point half-widths indicate highly precise index
  estimates.
- **Convergence:** Indices are stable to within 0.03% between the 75% and 100%
  subsample sizes.
- **Degree-3 LASSO cross-check:** The degree-3 PCE (with LASSO sparsification)
  achieves $Q^2 = 0.999$, confirming that higher-order effects exist but are too
  small to improve on the degree-2 surrogate for Sobol index estimation.

## 9. Computational Details

| Item | Value |
|---|---|
| Software | Julia + ThermaKin (Pyrolysis.jl) |
| PCE library | PolyChaos.jl |
| Sampling | QuasiMonteCarlo.jl (Sobol sequences) |
| Number of samples | 2048 |
| Simulation time per sample | ~7 min (240 cells, ALE mesh) |
| Total simulation wall time | ~5.8 hours (8 parallel workers) |
| ODE solver | QNDF (KLU sparse direct linear solver) |
| PCE degree selected | 2 (171 basis functions) |
| Bootstrap replicates | 1000 |
| Convergence subsamples | 5 levels (12.5%–100%) |

## 10. Conclusions

A comprehensive, GUM-traceable global sensitivity analysis of a Douglas Fir
pyrolysis model was performed using Polynomial Chaos Expansion with 2048 Sobol
quasi-random samples. The analysis demonstrates that:

1. The heats of pyrolysis reaction (five parameters from a two-stage WLS
   optimization) dominate the prediction uncertainty for all temperature-related
   outputs (99.9%) and contribute substantially (43–62%) to mass-related outputs.

2. The DSC enthalpy calibration bias (a single multiplicative parameter at 2.6%
   relative uncertainty) dominates the timing uncertainty (80.5%) and is the
   primary contributor to mass loss and final thickness uncertainty (56%).

3. The remaining six measurement sources—HFM calibration, water calibration,
   virgin and char specific heat regressions, and heat of sorption—each
   contribute less than 2% to every output, making them low-priority targets for
   experimental refinement.

4. Inter-group interactions are negligible ($< 0.12\%$), confirming that the
   uncertainty propagation is essentially additive within the GUM parameter
   bounds.

5. The degree-2 PCE surrogate achieves $Q^2 > 0.999$ for all QoIs, providing a
   near-exact analytical variance decomposition at negligible computational cost
   relative to Monte Carlo estimation.

These results provide a quantitative basis for prioritizing experimental effort:
improving the heats of reaction estimates and the DSC calibration would address
over 98% of the prediction uncertainty across all quantities of interest.

## References

1. JCGM 100:2008, *Evaluation of measurement data — Guide to the expression of
   uncertainty in measurement* (GUM).

2. Sobol', I. M. (2001). "Global sensitivity indices for nonlinear mathematical
   models and their Monte Carlo estimates." *Mathematics and Computers in
   Simulation*, 55(1–3), 271–280.

3. Sudret, B. (2008). "Global sensitivity analysis using polynomial chaos
   expansions." *Reliability Engineering & System Safety*, 93(7), 964–979.

4. Blatman, G. & Sudret, B. (2011). "Adaptive sparse polynomial chaos expansion
   based on least angle regression." *Journal of Computational Physics*, 230(6),
   2345–2367.

5. Xiu, D. & Karniadakis, G. E. (2002). "The Wiener–Askey polynomial chaos for
   stochastic differential equations." *SIAM Journal on Scientific Computing*,
   24(2), 619–644.

6. Saltelli, A. et al. (2008). *Global Sensitivity Analysis: The Primer*. John
   Wiley & Sons.

## Appendix A: Generated Output Files

| File | Description |
|---|---|
| `sobol_group_indices.csv` | First-order and total group Sobol indices for all QoIs |
| `sobol_bootstrap_ci.csv` | 95% bootstrap CIs for all group × QoI combinations |
| `pce_validation.csv` | PCE fit quality metrics ($R^2$, $Q^2$, oversampling) per QoI |
| `pce_coefficients.csv` | Full PCE coefficient vectors with multi-index labels |
| `sobol_variance_decomposition.png` | Stacked bar chart of group variance contributions |
| `sobol_bar_<qoi>.png` | Individual bar charts of Sobol indices per QoI |
| `pce_validation_<qoi>.png` | PCE prediction vs. simulation scatter plots |
| `sobol_convergence_<qoi>.png` | Convergence of Sobol indices with sample size |
| `srrc_sensitivity.csv` | Standardized Rank Regression Coefficients (screening) |
