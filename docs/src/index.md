# EquationsOfState.jl

!!! note

    Starting from `v2.0.0`, `EquationsOfState.jl` has been rewritten.
    The former behaviors of `EquationsOfState.jl` (`v1.1.4` and earlier) will be deprecated.
    Please follow the latest documentation.

## Package Features

- Calculate energy, pressure, and bulk modulus of an `EquationOfState` on a (an)
  volume (array of volumes).
- Fit an `EquationOfState` on a series of volumes using least-squares fitting
  method.
- Fit an `EquationOfState` on a series of volumes linearly.
- Find the corresponding volume of an `EquationOfState` given an (a) energy,
  pressure, and bulk modulus.

See the [Index](@ref main-index) for the complete list of documented functions
and types.

## Manual Outline

```@contents
Pages = [
    "Collections.md",
    "NonlinearFitting.md",
    "Find.md",
]
Depth = 3
```

### [Index](@id main-index)

```@index

```
