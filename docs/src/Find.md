# Find

```@meta
CurrentModule = EquationsOfState.Find
```

This module contains a function `findvolume`, which is used to find an approximate
volume at a given pressure, energy, or bulk modulus based on an equation of state.
A result is not always guaranteed, especially when the equation of state is not a
monotonic function of volume. However, according to experience, ``P(V)`` relation
is usually a monotonic function. So we suggest using `PressureForm` to find the
corresponding volume.

## Usage

```julia
julia> using EquationsOfState, EquationsOfState.Collections, EquationsOfState.Find, Unitful, UnitfulAtomic

julia> pressures = collect(0:20:200) .* u"GPa";

julia> eos = BirchMurnaghan3rd(167u"angstrom^3", 2600u"kbar", 4.0);

julia> volumes = map(
           p -> findvolume(PressureForm(), eos, p, (eps() * u"bohr^3", eos.v0 * 1.3)),
           pressures
       )
[ Info: Using method "Roots.Bisection"...
[ Info: Using method "Roots.Bisection"...
[ Info: Using method "Roots.Bisection"...
[ Info: Using method "Roots.Bisection"...
[ Info: Using method "Roots.Bisection"...
[ Info: Using method "Roots.Bisection"...
[ Info: Using method "Roots.Bisection"...
[ Info: Using method "Roots.Bisection"...
[ Info: Using method "Roots.Bisection"...
[ Info: Using method "Roots.Bisection"...
[ Info: Using method "Roots.Bisection"...
11-element Array{Quantity{Float64,𝐋^3,Unitful.FreeUnits{(Å^3,),𝐋^3,nothing}},1}:
              167.0 Å^3
 156.14036210727835 Å^3
 147.99803635986564 Å^3
 141.51093713795865 Å^3
 136.13864615965332 Å^3
 131.56784031939347 Å^3
 127.60046278645824 Å^3
 124.10332447387113 Å^3
 120.98257680606459 Å^3
 118.16962836248427 Å^3
 115.61284838696814 Å^3
```

Here we let the algorithm choose the bisection root-finding method to find the
`volumes` corresponding to `pressures`.

A figure is plotted below to verify our results, and it fits very well.

![findvolume](assets/findvolume.png)

## Public interfaces

```@docs
findvolume(form::EquationForm, eos::EquationOfState, y, x0, method)
```
