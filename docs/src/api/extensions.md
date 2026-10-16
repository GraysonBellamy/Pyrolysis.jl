```@meta
CurrentModule = Pyrolysis
```

# Package extensions

Extensions load automatically when the corresponding optional dependency is
present in the environment:

| Extension | Trigger | Provides |
|:---|:---|:---|
| `PyrolysisPlotsExt` | `using RecipesBase` (or `Plots`) | Plot recipes for `PyrolysisSolution` |
| `PyrolysisSciMLSensitivityExt` | `using SciMLSensitivity, Zygote` | [`forward_sensitivity`](@ref) / [`adjoint_sensitivity`](@ref) |

Calling a sensitivity stub without its trigger packages loaded raises a
`MethodError`. See the User Guide's
[extensions section](../guide/user_guide/11_solver-jacobian-config.md).
