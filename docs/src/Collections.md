# Collections

The current `EquationOfState`s contain

```
EquationOfState
├─ AntonSchmidt
├─ BreenanStacey
├─ FiniteStrainEquationOfState
│  ├─ BirchMurnaghan2nd
│  ├─ BirchMurnaghan3rd
│  ├─ BirchMurnaghan4th
│  ├─ PoirierTarantola2nd
│  ├─ PoirierTarantola3rd
│  └─ PoirierTarantola4th
├─ Murnaghan
├─ PolynomialEquationOfState
└─ Vinet
```

## Guide

We will use `BirchMurnaghan3rd` as an example.

`BirchMurnaghan3rd` can be constructed from scratch:

```jldoctest
julia> BirchMurnaghan3rd(1, 2, 3)
4-element BirchMurnaghan3rd{Int64}:
 1
 2
 3
 0

julia> BirchMurnaghan3rd(1, 2, 3, 4)
4-element BirchMurnaghan3rd{Int64}:
 1
 2
 3
 4

julia> BirchMurnaghan3rd(1, 2, 3, 4.0)
4-element BirchMurnaghan3rd{Float64}:
 1.0
 2.0
 3.0
 4.0
```

It can also be constructed from an existing `BirchMurnaghan3rd`:

```jldoctest
julia> BirchMurnaghan3rd(BirchMurnaghan3rd(1, 2, 3, 4.0), b0=10, e0=5)
4-element Birch{Float64}:
  1.0
 10.0
  3.0
  5.0

julia> BirchMurnaghan3rd(BirchMurnaghan3rd(1, 2, 3, 4.0), Dict(:b0=>10, :e0=>5))
4-element Birch{Float64}:
  1.0
 10.0
  3.0
  5.0

julia> BirchMurnaghan3rd(BirchMurnaghan3rd(1, 2, 3, 4.0), (:b0, 10))
4-element Birch{Float64}:
  1.0
 10.0
  3.0
  4.0
```

Users can access `BirchMurnaghan3rd`'s element by either "dot notation" or indexing:

```jldoctest
julia> b = BirchMurnaghan3rd(1, 2, 3, 4.0)
4-element BirchMurnaghan3rd{Float64}:
 1.0
 2.0
 3.0
 4.0

julia> b.v0
1.0

julia> b[1]
1.0
```
