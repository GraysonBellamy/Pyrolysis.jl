```@meta
CurrentModule = Pyrolysis
```

# Package extensions

Extensions load automatically when the corresponding optional dependency is
present in the environment:

| Extension | Trigger | Provides |
|:---|:---|:---|
| `PyrolysisKLUExt` | `using KLU` | [`klu_factorization`](@ref) linear-solver hook |
| `PyrolysisAMGExt` | `using AlgebraicMultigrid` | [`amg_preconditioner`](@ref) |
| `PyrolysisPardisoExt` | `using Pardiso` | [`pardiso_factorization`](@ref) |
| `PyrolysisPlotsExt` | `using RecipesBase` (or `Plots`) | Plot recipes for `PyrolysisSolution` |
| `PyrolysisSciMLSensitivityExt` | `using SciMLSensitivity, Zygote` | [`forward_sensitivity`](@ref) / [`adjoint_sensitivity`](@ref) |
| `PyrolysisSymbolicsExt` | `using Symbolics` | [`verify_mixing_derivatives`](@ref) (test-only) |

Calling a stub without its trigger dependency loaded raises a `MethodError`
naming the package to load. See the User Guide's
[extensions section](../guide/user_guide/11_solver-jacobian-config.md).

```@docs
klu_factorization
pardiso_factorization
amg_preconditioner
```
