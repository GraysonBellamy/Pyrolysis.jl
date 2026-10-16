```@meta
CurrentModule = Pyrolysis
```

# Package extensions

Extensions load automatically when the corresponding optional dependency is
present in the environment:

| Extension | Trigger | Provides |
|:---|:---|:---|
| `PyrolysisPlotsExt` | `using RecipesBase` (or `Plots`) | Plot recipes for `PyrolysisSolution` |

See the User Guide's
[extensions section](../guide/user_guide/11_solver-jacobian-config.md).
