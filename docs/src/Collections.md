# Collections

```@meta
CurrentModule = EquationsOfState.Collections
```

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

## Types

```@docs
EquationOfState
FiniteStrainEquationOfState
Murnaghan
BirchMurnaghan2nd
BirchMurnaghan3rd
BirchMurnaghan4th
PoirierTarantola2nd
PoirierTarantola3rd
PoirierTarantola4th
Vinet
```

## Usage

### Construct an `EquationOfState`
We will use `BirchMurnaghan3rd` as an example.

`BirchMurnaghan3rd` can be constructed from scratch:

```julia
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

```julia
julia> BirchMurnaghan3rd(BirchMurnaghan3rd(1, 2, 3, 4.0), b0=10, e0=5)
4-element BirchMurnaghan3rd{Float64}:
  1.0
 10.0
  3.0
  5.0

julia> BirchMurnaghan3rd(BirchMurnaghan3rd(1, 2, 3, 4.0), Dict(:b0=>10, :e0=>5))
4-element BirchMurnaghan3rd{Float64}:
  1.0
 10.0
  3.0
  5.0

julia> BirchMurnaghan3rd(BirchMurnaghan3rd(1, 2, 3, 4.0), (:b0, 10))
4-element BirchMurnaghan3rd{Float64}:
  1.0
 10.0
  3.0
  4.0
```

Users can access `BirchMurnaghan3rd`'s element by either "dot notation" or indexing:

```julia
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

### Calculate energies on an `EquationOfState`

The $E(V)$ relation of equations of state are listed as below:

1. `Murnaghan`:
   ```math
   E(V) = E_{0}+K_{0} V_{0}\left[\frac{1}{K_{0}^{\prime}\left(K_{0}^{\prime}-1\right)}\left(\frac{V}{V_{0}}\right)^{1-K_{0}^{\prime}}+\frac{1}{K_{0}^{\prime}} \frac{V}{V_{0}}-\frac{1}{K_{0}^{\prime}-1}\right].
   ```
   
2. `BirchMurnaghan2nd`:
   
   ```math
   E(V) = E_{0} + \frac{9}{8} B_{0} V_{0} \left(\left( V / V_0 \right)^{-2 / 3}-1\right)^{2}.
   ```
   
3. `BirchMurnaghan3rd`:

   ```math
   E(V) = E_{0}+\frac{9}{16} V_{0} B_{0} \frac{\left(x^{2 / 3}-1\right)^{2}}{x^{7 / 3}}\left\{x^{1 / 3}\left(B_{0}^{\prime}-4\right)-x\left(B_{0}^{\prime}-6\right)\right\}.
   ```

   where ``x = V / V_0``, and ``f = \frac{ 1 }{ 2 } \bigg[ \bigg( \frac{ V_0 }{ V } \bigg)^{2/3} - 1 \bigg]``.

4. `BirchMurnaghan4th`:

   ```math
   E(V) = E_{0}+\frac{3}{8} V_{0} B_{0} f^{2}\left[\left(9 H-63 B_{0}^{\prime}+143\right) f^{2}+12\left(B_{0}^{\prime}-4\right) f+12\right].
   ```

   where ``H = B_0 B_0'' + (B_0')^2``.

5. `PoirierTarantola2nd`:

   ```math
   E(V) = E_{0}+\frac{1}{2} B_{0} V_{0} \ln ^{2} x.
   ```

6. `PoirierTarantola3rd`:

   ```math
   E(V) = E_{0}+\frac{1}{6} B_{0} V_{0} \ln ^{2} x\left[\left(B_{0}^{\prime}+2\right) \ln x+3\right].
   ```

   
## Public interfaces

```@docs
apply(::EnergyForm, eos::EquationOfState)
apply(::EnergyForm, eos::Murnaghan, v::Real)
apply(::EnergyForm, eos::BirchMurnaghan2nd, v::Real)
apply(::EnergyForm, eos::BirchMurnaghan3rd, v::Real)
apply(::EnergyForm, eos::BirchMurnaghan4th, v::Real)

apply(::PressureForm, eos::EquationOfState)

apply(::BulkModulusForm, eos::EquationOfState)
```


