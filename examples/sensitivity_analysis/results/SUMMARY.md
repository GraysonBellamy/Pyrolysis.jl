# Sensitivity Analysis Results: Douglas Fir Cone Calorimeter

## What We Did

We propagated uncertainties in 16 thermochemical input parameters through a
ThermaKin3 simulation of Douglas Fir pyrolysis in a cone calorimeter (50 kW/m²,
12 mm thick, 6% moisture content). Kinetic parameters, thermal conductivity, and
boundary conditions were held fixed at their nominal values.

The 16 uncertain parameters come from 7 independent measurement sources:

| # | Measurement source          | Parameters | What was measured                          |
|---|-----------------------------|-----------:|-------------------------------------------|
| 1 | HFM calibration bias        |          1 | Heat flow meter systematic error (0.7%)    |
| 2 | DSC calibration bias        |          1 | DSC/STA systematic error (2.6%)            |
| 3 | Virgin Cp, low-T (HFM)     |          2 | Cp regression coefficients, T < 323 K     |
| 4 | Virgin Cp, high-T (DSC)    |          2 | Cp regression coefficients, T > 425 K     |
| 5 | Char Cp (DSC)               |          2 | Char Cp regression coefficients            |
| 6 | Heat of sorption (WLS fit)  |          3 | Exponential fit to sorption calorimetry    |
| 7 | Heats of reaction (bootstrap)|         5 | Two-stage optimization of 5 reaction enthalpies |

Parameters within each group are correlated (from regression covariances or
bootstrap), but the 7 groups are independent of each other.

We ran 1024 Sobol quasi-random samples (extended to 1600 total), then fitted a
degree-2 Polynomial Chaos Expansion (PCE) surrogate to each output quantity and
computed group Sobol indices to decompose output variance by measurement source.

## How Good Is the Surrogate?

The PCE fits are excellent for all outputs:

| Output                | R²     | Q² (LOO-CV) |
|-----------------------|--------|-------------|
| Peak MLR              | 1.0000 | 1.0000      |
| Time to peak MLR      | 0.9999 | 0.9999      |
| Total mass loss       | 1.0000 | 1.0000      |
| Peak surface temp     | 0.9998 | 0.9997      |
| Peak back temp        | 0.9995 | 0.9993      |
| Final thickness       | 1.0000 | 1.0000      |

R² > 0.999 and Q² > 0.999 for all quantities. The outputs are smooth, nearly
linear functions of the inputs in the sampled range, so a degree-2 polynomial
captures them almost perfectly.

## How Much Do the Outputs Vary?

Very little. The total output uncertainty from these 16 thermochemical parameters
is small:

| Output              | Mean      | 90% CI              | CV    |
|---------------------|-----------|---------------------|-------|
| Peak MLR            | 29.7 g/(m²s) | [28.2, 31.4]    | 3.3%  |
| Time to peak MLR    | 361 s     | [351, 371]          | 1.7%  |
| Total mass loss     | 82.2%     | [82.1, 82.4]        | 0.1%  |
| Peak surface temp   | 926.1 K   | [925.8, 926.5]      | 0.03% |
| Peak back temp      | 926.3 K   | [925.7, 926.9]      | 0.04% |
| Final thickness     | 7.39 mm   | [7.37, 7.42]        | 0.2%  |

Peak temperatures, total mass loss, and final thickness are essentially
deterministic given these input uncertainties. Even the most variable output,
peak MLR, varies by only ~3%.

## Where Does the Uncertainty Come From?

The group Sobol indices tell us what fraction of each output's variance is
attributable to each measurement source:

| Output              | Heats of reaction | DSC calibration | Everything else |
|---------------------|------------------:|----------------:|----------------:|
| Peak MLR            |           **81%** |             18% |             <1% |
| Time to peak MLR    |               20% |         **78%** |              2% |
| Total mass loss     |               47% |         **52%** |             <1% |
| Peak surface temp   |           **99%** |              1% |           <0.1% |
| Peak back temp      |           **99%** |              1% |           <0.1% |
| Final thickness     |               47% |         **52%** |             <1% |

Only two measurement sources matter:

**Heats of reaction (bootstrap optimization)** — This group dominates peak MLR
(81%) and almost entirely controls peak temperatures (99%). The five reaction
enthalpies are highly correlated (r > 0.9) because they share a common source of
uncertainty: the pseudo-char heat capacity estimated in Stage 1 of the two-stage
optimization. When that shared baseline shifts, all five enthalpies shift together.

**DSC calibration bias** — The 2.6% systematic uncertainty in the DSC heat flow
calibration dominates the timing of peak MLR (78%) and contributes roughly half
the variance in total mass loss and final thickness (52%). This is because the
DSC calibration scales all STA-derived heat capacities (virgin Cp above 425 K
and char Cp), which directly affects the thermal inertia of the sample and
therefore how quickly it pyrolyzes.

The other five measurement sources — HFM calibration, virgin Cp low-T regression,
char Cp regression, and heat of sorption — collectively contribute less than 2%
of the variance for any output. Their uncertainties are already small enough to
be irrelevant.

Interactions between groups are negligible (<0.1%), meaning the model responds
nearly additively to these input uncertainties in the sampled range.

## What This Means

1. **The model's predictions are robust to thermochemical input uncertainties.**
   Even with realistic measurement uncertainties propagated through, the outputs
   barely change. If the model disagrees with experiment, the cause is almost
   certainly in parameters that were held fixed (kinetics A and E, thermal
   conductivity, emissivity, convective heat transfer coefficient) rather than
   in these thermochemical properties.

2. **Reducing measurement uncertainty further has limited value.** The DSC
   calibration bias (2.6%) is the single largest contributor to output variance
   for most quantities, but even it only produces ~3% variation in peak MLR.
   Investing in better DSC calibration would marginally tighten the predictions,
   but the uncertainty is already small.

3. **The heats of reaction are the one thermochemical group worth refining if
   peak temperature predictions matter.** They explain 99% of the (already small)
   variance in peak surface and back temperatures. This uncertainty traces back
   to the pseudo-char Cp from Stage 1 of the two-stage optimization — improving
   that estimate (e.g., with direct char Cp measurements at high temperature)
   would tighten all five reaction enthalpies simultaneously.

4. **Heat of sorption, char Cp, and low-temperature virgin Cp can be treated as
   known.** Their measurement uncertainties have no meaningful effect on any
   output at this moisture content (6% dry basis) and heat flux (50 kW/m²).

## Method Details

- **Sampling**: 1024 Sobol quasi-random samples (extended to 1600), with within-group
  correlations induced via Cholesky-transformed marginals.
- **Sensitivity method**: Group Sobol indices computed analytically from a degree-2
  Polynomial Chaos Expansion (153 Hermite basis functions, 10.5x oversampling ratio).
  PCE coefficients fitted by ordinary least squares. Validated by leave-one-out
  cross-validation (Q² > 0.999 for all outputs).
- **Software**: ThermaKin3 (Pyrolysis.jl), PolyChaos.jl for Hermite polynomial basis.

## Files

| File | Description |
|------|-------------|
| `sensitivity_douglas_fir.jl` | Main simulation script (sampling + solver) |
| `sensitivity_sobol.jl` | Post-processing (PCE fitting + Sobol indices) |
| `results/parameters.csv` | All 1600 sampled parameter vectors |
| `results/scalar_qois.csv` | Scalar outputs for all runs |
| `results/sobol_group_indices.csv` | First-order and total group Sobol indices |
| `results/pce_validation.csv` | R² and Q² for each QoI |
| `results/pce_coefficients.csv` | Full PCE coefficient vectors |
| `results/sobol_variance_decomposition.png` | Stacked bar chart (all QoIs) |
| `results/sobol_bar_<qoi>.png` | Individual bar charts per QoI |
| `results/pce_validation_<qoi>.png` | Predicted vs actual scatter plots |
| `results/mlr_uncertainty_envelope.png` | MLR time series with 50%/90% CI |
| `results/surface_temp_uncertainty_envelope.png` | Surface temperature envelope |
| `results/scalar_qoi_histograms.png` | Histograms of scalar outputs |
