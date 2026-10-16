```@meta
CurrentModule = Pyrolysis
```

# Property functions

Temperature- (and state-) dependent component properties. Plain numbers,
tuples, and callables convert automatically where a property function is
expected; see the User Guide's
[property-function section](../guide/user_guide/03_defining-materials.md).

```@docs
ConstantProperty
LinearProperty
PolynomialProperty
ArrheniusProperty
TableProperty
CallableProperty
SutherlandProperty
StateDependentProperty
```
