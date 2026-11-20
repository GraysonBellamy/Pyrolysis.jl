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
| [`07_douglas_fir_cone.jl`](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/examples/07_douglas_fir_cone.jl) | Douglas fir with multi-step kinetics, the two-way phase change of its moisture on the CAST spruce isotherm (its vapour starting in equilibrium, a humid ambient), lumped wall transfer and Darcy flow, and a mesh that moves with the solid. |
| [`08_pressure_treated_wood_cone.jl`](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/examples/08_pressure_treated_wood_cone.jl) | Pressure-treated wood with multi-step Arrhenius kinetics, transcribed from ThermaKin with its one-way drying step. |
| [`macfp/macfp_gasification.jl`](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/examples/macfp/macfp_gasification.jl) | MaCFP gasification benchmark with parallel solves. |

## Transport and physics demonstrations

| Script | What it models |
|:---|:---|
| [`08_gas_transport.jl`](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/examples/08_gas_transport.jl) | Lumped and dusty-gas transport through a charring slab, with atmospheric reservoir exchange and species diagnostics. |
| [`convergence_douglas_fir.jl`](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/examples/convergence_douglas_fir.jl) | Grid-convergence study over mesh resolutions (uses `Distributed`). |

## Comparison and validation studies

| Script | What it models |
|:---|:---|
| [`ptw_comparison.jl`](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/examples/ptw_comparison.jl) | Pressure-treated wood at 50 kW/m²: Pyrolysis.jl against ThermaKin, from the results `08_pressure_treated_wood_cone.jl` writes. |

## Inverse analysis

[`inverse_analysis/`](https://github.com/GraysonBellamy/Pyrolysis.jl/tree/main/examples/inverse_analysis)
fits white-pine component conductivities and kinetics to cone-calorimeter
data: forward model, objective, optimisation driver, uncertainty
quantification, and plotting/diagnostic helpers. See the
[inverse-analysis walkthrough](guide/user_guide/14_worked-examples.md) and the
[sensitivity-analysis chapter](guide/user_guide/13_sensitivity-analysis.md)
for the underlying pattern.

## Saving results to a file

| Script | What it models |
|:---|:---|
| [`export_results_csv.jl`](https://github.com/GraysonBellamy/Pyrolysis.jl/blob/main/examples/export_results_csv.jl) | Writes a solution to two CSV files that open in Excel: surface and back-face temperature, mass, mass-loss rate and thickness against time, and the temperature and concentration profiles. Its `write_results_csv` function can be copied into any script. |

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
