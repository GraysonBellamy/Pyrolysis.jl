# ============================================================================
# Unified Parameter Evaluation Functions
# ============================================================================
#
# These functions provide a single, AD-compatible interface for evaluating
# parameters that may be constants, functions of time, or functions of
# temperature. They replace the previous duplicate implementations in
# boundary_fluxes.jl and radiation_boundary.jl.
#
# Key features:
# - Type-stable: preserves Dual number types for automatic differentiation
# - Unified API: same function names used throughout codebase
# - Handles Real, Function, and AbstractPropertyFunction types
#
# ============================================================================

"""
	evaluate_param(param, t)

Evaluate a time-dependent parameter. Handles:
- Constants (Real): Returns the value, preserving AD type from `t`
- Functions: Calls f(t)

# Notes
- Made generic to support automatic differentiation with Dual numbers
- For constants, preserves the type of `t` (important for AD compatibility)

# Example
```julia
evaluate_param(100.0, t)      # Returns 100.0 * one(t)
evaluate_param(t -> 50*t, t)  # Returns 50*t
```
"""
@inline evaluate_param(param::Real, t) = param * one(t)  # Preserve AD type
@inline evaluate_param(param::Function, t) = param(t)

"""
	evaluate_param_T(param, T)

Evaluate a temperature-dependent parameter. Handles:
- Constants (Real): Returns the value, preserving AD type from `T`
- Functions: Calls f(T)
- AbstractPropertyFunction: Calls the property function with T

# Notes
- Made generic to support automatic differentiation with Dual numbers
- For constants, preserves the type of `T` (important for AD compatibility)

# Example
```julia
evaluate_param_T(0.9, T)                   # Returns 0.9 * one(T)
evaluate_param_T(T -> 0.85 + 0.001*T, T)   # Returns temperature-dependent value
evaluate_param_T(my_property_func, T)      # Calls property function
```
"""
@inline evaluate_param_T(param::Real, T) = param * one(T)  # Preserve AD type
@inline evaluate_param_T(param::Function, T) = param(T)
@inline evaluate_param_T(param::AbstractPropertyFunction, T) = param(T)
