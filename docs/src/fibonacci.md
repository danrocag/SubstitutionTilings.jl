# The Fibonacci tiling

In this page, we explain how we can define a substitution tiling, and use this library to generate pictures of it.
We do this in the example of the Fibonacci tiling.


```@example 1
using SubstitutionTilings
using SubstitutionTilings.NumFields
using SubstitutionTilings.CoreDefs
using StructEquality
```

In order to be able to do frequency computation,
we need to work in a group with exact equality, so we can't use floats as coordinates.
We use ou