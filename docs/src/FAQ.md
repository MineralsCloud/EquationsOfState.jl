# FAQ

## How to make a `Vector` from an `EquationOfState`?

A suggested way is to use the
[`IterTools.fieldvalues` function](https://juliacollections.github.io/IterTools.jl/latest/index.html#IterTools.fieldvalues):

```julia
julia> using IterTools

julia> eos = BirchMurnaghan4th(1, 2.0, 3, 4)
BirchMurnaghan4th{Float64}(1.0, 2.0, 3.0, 4.0, 0.0)

julia> collect(fieldvalues(eos))
5-element Array{Any,1}:
 1.0
 2.0
 3.0
 4.0
 0.0
```

It is lazy and fast. But the `Vector` has element type `Any`.

Or to write a non-lazy version of `fieldvalues` manually:

```julia
julia> fieldvalues(eos::EquationOfState) = [getfield(eos, i) for i in 1:nfields(eos)]
fieldvalues (generic function with 1 method)

julia> Collections.fieldvalues(eos)
5-element Array{Float64,1}:
 1.0
 2.0
 3.0
 4.0
 0.0
```

It is slower than `IterTools.fieldvalues` but has the exact element type `Float64`.
Use it with care.
