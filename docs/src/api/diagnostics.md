```@meta
CurrentModule = Pyrolysis
```

# Conservation diagnostics

Runtime mass-, energy-, and per-reaction-ledger tracking. These names live in
`Pyrolysis.Diagnostics` (also re-exported by `Pyrolysis.Internal`); the
common entry point is the `diagnostics` keyword of [`solve`](@ref). See the
Technical Reference's
[conservation-diagnostics chapter](../guide/technical_reference/16_conservation-diagnostics.md).

```@autodocs
Modules = [Pyrolysis.Diagnostics]
```
