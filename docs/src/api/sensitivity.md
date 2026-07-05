```@meta
CurrentModule = Pyrolysis
```

# Sensitivity analysis

Forward- and adjoint-mode local sensitivities. These are extension stubs:
concrete methods attach when `SciMLSensitivity` (and `Zygote`) are loaded —
see [Extensions](extensions.md). The User Guide's
[sensitivity chapter](../guide/user_guide/13_sensitivity-analysis.md) covers
the `p_inject` pattern and worked examples; the Technical Reference's
[sensitivity theory chapter](../guide/technical_reference/17_sensitivity-analysis.md)
covers the underlying formulations.

```@docs
forward_sensitivity
adjoint_sensitivity
verify_mixing_derivatives
```
