# Examples

The repository ships end-to-end studies under
[`examples/`](https://github.com/GraysonBellamy/Pyrolysis.jl/tree/main/examples).
Each script is self-contained and runs in the examples environment
(`julia --project=examples`). The User Guide's
[worked-examples chapter](guide/user_guide/14_worked-examples.md) walks
through the physics and the code of the main studies in detail.

## Cone calorimeter

| Script | What it models |
|:---|:---|
| [`06_cone_calorimeter.jl`](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/examples/06_cone_calorimeter.jl) | PMMA slab under a 25 kW/m² cone heat flux with the incident radiation absorbed at the surface, single-step Arrhenius kinetics, and gas release, overlaid on ThermaKin reference outputs when they are present. |
| [`06_cone_calorimeter_surface_absorption.jl`](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/examples/06_cone_calorimeter_surface_absorption.jl) | The same PMMA case with the surface-absorption radiation model set explicitly; it absorbs the radiation at the exposed face as the default does. |
| [`07_douglas_fir_cone.jl`](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/examples/07_douglas_fir_cone.jl) | Douglas fir with multi-step kinetics, moisture, Darcy–Fick transport, and a mesh that moves with the solid. |
| [`08_pressure_treated_wood_cone.jl`](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/examples/08_pressure_treated_wood_cone.jl) | Pressure-treated wood with moisture, mapped from a ThermaKin calibration. |
| [`macfp/macfp_gasification.jl`](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/examples/macfp/macfp_gasification.jl) | MaCFP gasification benchmark with parallel solves. |

## Transport and physics demonstrations

| Script | What it models |
|:---|:---|
| [`08_darcy_fick_demo.jl`](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/examples/08_darcy_fick_demo.jl) | Pressure-driven Darcy–Fick gas transport through a charring, wood-like material. |
| [`moisture_hos_comparison.jl`](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/examples/moisture_hos_comparison.jl) | Effect of the moisture heat-of-sorption model on drying and pyrolysis. |
| [`convergence_douglas_fir.jl`](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/examples/convergence_douglas_fir.jl) | Grid-convergence study over mesh resolutions (uses `Distributed`). |

## Comparison and validation studies

| Script | What it models |
|:---|:---|
| [`ptw_comparison.jl`](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/examples/ptw_comparison.jl) | Pressure-treated wood results against ThermaKin output (includes a ThermaKin output parser). |
| [`df_sensitivity_comparison.jl`](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/examples/df_sensitivity_comparison.jl) | Inter-code sensitivity comparison for the Douglas fir case. |
| [`literature_comparison/literature_comparison.jl`](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/examples/literature_comparison/literature_comparison.jl) | Runs the cone cases across literature property sets via serialisable run specs. |

## Inverse analysis

[`inverse_analysis/`](https://github.com/GraysonBellamy/Pyrolysis.jl/tree/main/examples/inverse_analysis)
fits white-pine component conductivities and kinetics to cone-calorimeter
data: forward model, objective, optimisation driver, uncertainty
quantification, and plotting/diagnostic helpers. See the
[inverse-analysis walkthrough](guide/user_guide/14_worked-examples.md) and the
[sensitivity-analysis chapter](guide/user_guide/13_sensitivity-analysis.md)
for the underlying pattern.

## Sensitivity analysis

[`sensitivity_analysis/`](https://github.com/GraysonBellamy/Pyrolysis.jl/tree/main/examples/sensitivity_analysis)
contains forward uncertainty propagation via Latin Hypercube sampling
(`sensitivity_douglas_fir.jl`), global Sobol indices via polynomial chaos
expansion (`sensitivity_sobol.jl`), and a PCE degree-selection study
(`compare_degrees.jl`), together with write-ups of the results.

## Running the examples

```julia
# from the repository root
julia --project=examples examples/06_cone_calorimeter.jl
```

The examples environment adds Plots and other post-processing dependencies on
top of the package itself. Scripts that use `Distributed` (the convergence and
sweep studies) document their worker setup at the top of the file. Every
simulation script runs a short check of itself with `PYROLYSIS_SMOKE=1`; the
[examples README](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/examples/README.md)
lists each script's smoke command, its class and the data files it needs.
