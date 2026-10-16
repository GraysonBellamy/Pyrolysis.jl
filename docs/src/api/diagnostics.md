```@meta
CurrentModule = Pyrolysis
```

# Conservation diagnostics

Mass- and energy-conservation ledgers, advanced step by step with
`record_substep!`, and the sensible-energy log that the `diagnostics` keyword
of [`solve`](@ref) records in `extras.diagnostics_log`. These names live in
`Pyrolysis.Diagnostics` (also re-exported by `Pyrolysis.Internal`). See the
Technical Reference's
[conservation-diagnostics chapter](../guide/technical_reference/16_conservation-diagnostics.md).

```@autodocs
Modules = [Pyrolysis.Diagnostics]
```
