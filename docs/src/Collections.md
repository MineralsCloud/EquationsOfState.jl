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

### Calculate energy on an `EquationOfState`

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

7. `PoirierTarantola4th`:

   ```math
   E(V) = E_{0}+\frac{1}{24} B_{0} V_{0} \ln ^{2} x\left\{\left(H+3 B_{0}^{\prime}+3\right) \ln ^{2} x\right. \left.+4\left(B_{0}^{\prime}+2\right) \ln x+12\right\}.
   ```
   where ``H = B_0 B_0'' + (B_0')^2``.

8. `Vinet`:

   ```math
   E(V) = E_{0}+\frac{9}{16} V_{0} B_{0} \frac{\left(x^{2 / 3}-1\right)^{2}}{x^{7 / 3}}\left\{x^{1 / 3}\left(B_{0}^{\prime}-4\right)-x\left(B_{0}^{\prime}-6\right)\right\}.
   ```

9. `AntonSchmidt`:

   ```math
   E(V)=\frac{\beta V_{0}}{n+1}\left(\frac{V}{V_{0}}\right)^{n+1}\left[\ln \left(\frac{V}{V_{0}}\right)-\frac{1}{n+1}\right]+E_{\infty}.
   ```

### Calculate pressure on an `EquationOfState`

The $P(V)$ relation of equations of state are listed as below:

1. `Murnaghan`:

   ```math
   P(V) = \frac{B_{0}}{B_{0}^{\prime}}\left[\left(\frac{V_{0}}{V}\right)^{B_{0}^{\prime}}-1\right].
   ```

2. `BirchMurnaghan2nd`:

   ```math
   P(V) = \frac{3}{2} B_{0}\left(x^{-7 / 3}-x^{-5 / 3}\right).
   ```

3. `BirchMurnaghan3rd`:

   ```math
   P(V) = \frac{3}{8} B_{0} \frac{x^{2 / 3}-1}{x^{10 / 3}}\left\{3 B_{0}^{\prime} x-16 x-3 x^{1 / 3}\left(B_{0}^{\prime}-4\right)\right\}.
   ```

4. `BirchMurnaghan4th`:

   ```math
   P(V) = \frac{1}{2} B_{0}(2 f+1)^{5 / 2}\left\{\left(9 H-63 B_{0}^{\prime}+143\right) f^{2}\right.\left.+9\left(B_{0}^{\prime}-4\right) f+6\right\}.
   ```

5. `PoirierTarantola2nd`:

   ```math
   P(V) = -\frac{B_{0}}{x} \ln x.
   ```

6. `PoirierTarantola3rd`:

   ```math
   P(V) = -\frac{B_{0} \ln x}{2 x}\left[\left(B_{0}^{\prime}+2\right) \ln x+2\right].
   ```

7. `PoirierTarantola4th`:

   ```math
   P(V) = -\frac{B_{0} \ln x}{6 x}\left\{\left(H+3 B_{0}^{\prime}+3\right) \ln ^{2} x+3\left(B_{0}^{\prime}+6\right) \ln x+6\right\}.
   ```

8. `Vinet`:

   ```math
   P(V) = 3 B_{0} \frac{1-\eta}{\eta^{2}} \exp \left\{-\frac{3}{2}\left(B_{0}^{\prime}-1\right)(\eta-1)\right\}.
   ```

9. `AntonSchmidt`:

   ```math
   P(V) = -\beta\left(\frac{V}{V_{0}}\right)^{n} \ln \left(\frac{V}{V_{0}}\right).
   ```

### Calculate bulk modulus on an `EquationOfState`

The $B(V)$ relation of equations of state are listed as below:

1. `BirchMurnaghan2nd`:

   ```math
   B(V) = B_{0}(7 f+1)(2 f+1)^{5 / 2}.
   ```

2. `BirchMurnaghan3rd`:

   ```math
   B(V) = \frac{B_{0}}{8 \chi^{10 / 3}}\left\{x^{5 / 3}\left(15 B_{0}^{\prime}-80\right)-x\left(42 B_{0}^{\prime}-196\right)\right.\left.+27 x^{1 / 3}\left(B_{0}^{\prime}-4\right)\right\}.
   ```

3. `BirchMurnaghan4th`:

   ```math
   B(V) = \frac{1}{6} B_{0}(2 f+1)^{5 / 2}\left\{\left(99 H-693 B_{0}^{\prime}+1573\right) f^{3}\right.\left.+\left(27 H-108 B_{0}^{\prime}+105\right) f^{2}+6\left(3 B_{0}^{\prime}-5\right) f+6\right\}.
   ```

4. `PoirierTarantola2nd`:

   ```math
   B(V) = \frac{B_{0}}{x}(1-\ln x).
   ```

5. `PoirierTarantola3rd`:

   ```math
   B(V) = -\frac{B_{0}}{2 x}\left[\left(B_{0}^{\prime}+2\right) \ln x(\ln x-1)-2\right].
   ```

6. `PoirierTarantola4th`:

   ```math
   B(V) = -\frac{B_{0}}{6 x}\left\{\left(H+3 B_{0}^{\prime}+3\right) \ln ^{3} x-3\left(H+2 B_{0}^{\prime}+1\right) \ln ^{2} x\right.\left.-6\left(B_{0}^{\prime}+1\right) \ln x-6\right\}.
   ```

7. `Vinet`:

   ```math
   B(V) = -\frac{B_{0}}{2 \eta^{2}}\left[3 \eta(\eta-1)\left(B_{0}^{\prime}-1\right)+2(\eta-2)\right]\times \exp \left\{-\frac{3}{2}\left(B_{0}^{\prime}-1\right)(\eta-1)\right\}.
   ```

8. `AntonSchmidt`:

   ```math
   B(V) = \beta\left(\frac{V}{V_{0}}\right)^{n}\left[1+n \ln \frac{V}{V_{0}}\right].
   ```


## Public interfaces

```@docs
apply(::EnergyForm, eos::EquationOfState)
apply(::EnergyForm, eos::Murnaghan, v::Real)
apply(::EnergyForm, eos::BirchMurnaghan2nd, v::Real)
apply(::EnergyForm, eos::BirchMurnaghan3rd, v::Real)
apply(::EnergyForm, eos::BirchMurnaghan4th, v::Real)
apply(::EnergyForm, eos::PoirierTarantola2nd, v::Real)
apply(::EnergyForm, eos::PoirierTarantola3rd, v::Real)
apply(::EnergyForm, eos::PoirierTarantola4th, v::Real)
apply(::EnergyForm, eos::Vinet, v::Real)
apply(::EnergyForm, eos::AntonSchmidt, v::Real)
apply(::PressureForm, eos::EquationOfState)
apply(::PressureForm, eos::Murnaghan, v::Real)
apply(::PressureForm, eos::BirchMurnaghan2nd, v::Real)
apply(::PressureForm, eos::BirchMurnaghan3rd, v::Real)
apply(::PressureForm, eos::BirchMurnaghan4th, v::Real)
apply(::PressureForm, eos::PoirierTarantola2nd, v::Real)
apply(::PressureForm, eos::PoirierTarantola3rd, v::Real)
apply(::PressureForm, eos::PoirierTarantola4th, v::Real)
apply(::PressureForm, eos::Vinet, v::Real)
apply(::PressureForm, eos::AntonSchmidt, v::Real)
apply(::PressureForm, eos::BreenanStacey, v::Real)
apply(::BulkModulusForm, eos::EquationOfState)
apply(::BulkModulusForm, eos::BirchMurnaghan2nd, v::Real)
apply(::BulkModulusForm, eos::BirchMurnaghan3rd, v::Real)
apply(::BulkModulusForm, eos::BirchMurnaghan4th, v::Real)
apply(::BulkModulusForm, eos::PoirierTarantola2nd, v::Real)
apply(::BulkModulusForm, eos::PoirierTarantola3rd, v::Real)
apply(::BulkModulusForm, eos::PoirierTarantola4th, v::Real)
apply(::BulkModulusForm, eos::Vinet, v::Real)
apply(::BulkModulusForm, eos::AntonSchmidt, v::Real)
```
