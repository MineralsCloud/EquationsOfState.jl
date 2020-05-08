"""
This module provides `EquationOfState` types and calculate
energy, pressure, or bulk modulus of an `EquationOfState` on
a (an) volume (array of volumes).
"""
module Collections

using IterTools: FieldValues, fieldvalues
using Unitful: AbstractQuantity, dimension, upreferred, @u_str

import Unitful

export Energy,
    Pressure,
    BulkModulus,
    Murnaghan,
    BirchMurnaghan2nd,
    BirchMurnaghan3rd,
    BirchMurnaghan4th,
    PoirierTarantola2nd,
    PoirierTarantola3rd,
    PoirierTarantola4th,
    Vinet,
    AntonSchmidt,
    BreenanStacey,
    Shanker

# ============================================================================ #
#                                     Types                                    #
# ============================================================================ #
abstract type PhysicalProperty end
"""
    Energy()
    (::EquationOfState)(::Energy)(v)
    (::EquationOfState)(::Energy)

Return the energy of an `EquationOfState` on volume `v`. If `eos` has units,
`v` must also has.

Return a [function-like object](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects-1) that takes a volume as a variable, suitable for mapping onto an array.

# Examples
```jldoctest
julia> f = Vinet(1, 2, 3)(Energy());

julia> map(f, 1:1:10)
10-element Array{Float64,1}:
 0.0
 0.367905230584308
 0.7652477289745814
 1.0516459435179233
 1.2560420090256408
 1.405149833626178
 1.5165867441792136
 1.6017034530570884
 1.6679539823686644
 1.7203642945516917
```

However, these methods are preserved for special cases
(see [#52](https://github.com/MineralsCloud/EquationsOfState.jl/issues/52#issuecomment-555856194)).
In most cases, the Julia [`do` block syntax](http://docs.julialang.org/en/v1/base/base/#do)
is preferred:
```jldoctest
julia> map(1:1:10) do v
           Vinet(1, 2, 3)(Energy())(v)
       end
10-element Array{Float64,1}:
 0.0
 0.367905230584308
 0.7652477289745814
 1.0516459435179235
 1.2560420090256412
 1.405149833626178
 1.5165867441792138
 1.6017034530570884
 1.6679539823686644
 1.7203642945516917
```
"""
struct Energy <: PhysicalProperty end
"""
    Pressure()
    (::EquationOfState)(::Pressure)(v)
    (::EquationOfState)(::Pressure)

Return the pressure of an `EquationOfState` on volume `v`. If `eos` has units,
`v` must also has.

# Examples
```jldoctest
julia> f = Vinet(1, 2, 3)(Pressure());

julia> map(f, 1:1:10)
10-element Array{Float64,1}:
  0.0
 -0.45046308428750254
 -0.3384840350043251
 -0.24010297221667418
 -0.17314062272722755
 -0.12795492664586872
 -0.09677154467733216
 -0.07468060255179591
 -0.05864401631176751
 -0.04674768462396211
```
"""
struct Pressure <: PhysicalProperty end
"""
    BulkModulus()
    (::EquationOfState)(::BulkModulus)(v)
    (::EquationOfState)(::BulkModulus)

Return the bulk modulus of an `EquationOfState` on volume `v`. If `eos` has units,
`v` must also has.

# Examples
```jldoctest
julia> f = BirchMurnaghan3rd(1, 2, 3)(BulkModulus());

julia> map(f, 1:1:10)
10-element Array{Float64,1}:
 2.0
 0.9216086833346415
 0.444903691617472
 0.2540009203153288
 0.16193296566524193
 0.11130584492987289
 0.08076305569984538
 0.06103515625
 0.047609811583958425
 0.03808959181078831
```
"""
struct BulkModulus <: PhysicalProperty end

"""
    EquationOfState{T}

An abstraction of equations of state, where `T` specifies the elements' common type.
"""
abstract type EquationOfState{T} end

"""
    FiniteStrainEOS{T} <: EquationOfState{T}

An abstraction of finite strain equations of state, where `T` specifies the elements' common type.
"""
abstract type FiniteStrainEOS{T} <: EquationOfState{T} end

"""
    Murnaghan(v0, b0, b′0, e0)

Create a Murnaghan equation of state. The elements' type will be handled automatically.

This equation of state can have units. The units are specified in [`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)'s
`@u_str` style.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

# Examples
```jldoctest
julia> Murnaghan(1, 2, 3.0)
Murnaghan{Float64}(1.0, 2.0, 3.0, 0.0)

julia> Murnaghan(Int8(1), 2//1, 3.0, 4)
Murnaghan{Float64}(1.0, 2.0, 3.0, 4.0)

julia> Murnaghan(1u"nm^3", 2u"GPa", 3, 3.0u"eV")
Murnaghan{Quantity{Float64,D,U} where U where D}(1.0 nm³, 2.0 GPa, 3.0, 3.0 eV)
```
"""
struct Murnaghan{T} <: EquationOfState{T}
    v0::T
    b0::T
    b′0::T
    e0::T
end
function Murnaghan(v0, b0, b′0, e0)
    T = Base.promote_typeof(v0, b0, b′0, e0)
    return Murnaghan{T}(convert.(T, [v0, b0, b′0, e0])...)
end
Murnaghan(v0::Real, b0::Real, b′0::Real) = Murnaghan(v0, b0, b′0, 0)
Murnaghan(v0::AbstractQuantity, b0::AbstractQuantity, b′0) =
    Murnaghan(v0, b0, b′0, 0 * upreferred(Unitful.J))

"""
    BirchMurnaghan2nd(v0, b0, e0)

Create a Birch–Murnaghan 2nd order equation of state. The elements' type will be handled automatically.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

See also: [`BirchMurnaghan3rd`](@ref), [`BirchMurnaghan4th`](@ref)

# Examples
```jldoctest
julia> BirchMurnaghan2nd(1, 2.0)
BirchMurnaghan2nd{Float64}(1.0, 2.0, 0.0)

julia> BirchMurnaghan2nd(Int8(1), 2//1, 0.0)
BirchMurnaghan2nd{Float64}(1.0, 2.0, 0.0)

julia> BirchMurnaghan2nd(1u"nm^3", 2u"GPa", 3.0u"eV")
BirchMurnaghan2nd{Quantity{Float64,D,U} where U where D}(1.0 nm³, 2.0 GPa, 3.0 eV)
```
"""
struct BirchMurnaghan2nd{T} <: FiniteStrainEOS{T}
    v0::T
    b0::T
    e0::T
end
function BirchMurnaghan2nd(v0, b0, e0)
    T = Base.promote_typeof(v0, b0, e0)
    return BirchMurnaghan2nd{T}(convert.(T, [v0, b0, e0])...)
end
BirchMurnaghan2nd(v0::Real, b0::Real) = BirchMurnaghan2nd(v0, b0, 0)
BirchMurnaghan2nd(v0::AbstractQuantity, b0::AbstractQuantity) =
    BirchMurnaghan2nd(v0, b0, 0 * upreferred(Unitful.J))

"""
    BirchMurnaghan3rd(v0, b0, b′0, e0)

Create a Birch–Murnaghan 3rd order equation of state. The elements' type will be handled automatically.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

See also: [`BirchMurnaghan2nd`](@ref), [`BirchMurnaghan4th`](@ref)

# Examples
```jldoctest
julia> BirchMurnaghan3rd(1, 2.0, 3)
BirchMurnaghan3rd{Float64}(1.0, 2.0, 3.0, 0.0)

julia> BirchMurnaghan3rd(Int8(1), 2//1, 4, 0.0)
BirchMurnaghan3rd{Float64}(1.0, 2.0, 4.0, 0.0)

julia> BirchMurnaghan3rd(1u"nm^3", 2u"GPa", 4.0, 3u"eV")
BirchMurnaghan3rd{Quantity{Float64,D,U} where U where D}(1.0 nm³, 2.0 GPa, 4.0, 3.0 eV)
```
"""
struct BirchMurnaghan3rd{T} <: FiniteStrainEOS{T}
    v0::T
    b0::T
    b′0::T
    e0::T
end
function BirchMurnaghan3rd(v0, b0, b′0, e0)
    T = Base.promote_typeof(v0, b0, b′0, e0)
    return BirchMurnaghan3rd{T}(convert.(T, [v0, b0, b′0, e0])...)
end
BirchMurnaghan3rd(v0::Real, b0::Real, b′0::Real) = BirchMurnaghan3rd(v0, b0, b′0, 0)
BirchMurnaghan3rd(v0::AbstractQuantity, b0::AbstractQuantity, b′0) =
    BirchMurnaghan3rd(v0, b0, b′0, 0 * upreferred(Unitful.J))

"""
    BirchMurnaghan4th(v0, b0, b′0, b′′0, e0)

Create a Birch–Murnaghan 4th order equation of state. The elements' type will be handled automatically.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `b′′0`: the second-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

See also: [`BirchMurnaghan2nd`](@ref), [`BirchMurnaghan4th`](@ref)

# Examples
```jldoctest
julia> BirchMurnaghan4th(1, 2.0, 3, 4)
BirchMurnaghan4th{Float64}(1.0, 2.0, 3.0, 4.0, 0.0)

julia> BirchMurnaghan4th(Int8(1), 2//1, 4, 5.0, Float16(6))
BirchMurnaghan4th{Float64}(1.0, 2.0, 4.0, 5.0, 6.0)

julia> BirchMurnaghan4th(1u"nm^3", 2u"GPa", 3.0, 4u"1/GPa", 5u"eV")
BirchMurnaghan4th{Quantity{Float64,D,U} where U where D}(1.0 nm³, 2.0 GPa, 3.0, 4.0 GPa⁻¹, 5.0 eV)
```
"""
struct BirchMurnaghan4th{T} <: FiniteStrainEOS{T}
    v0::T
    b0::T
    b′0::T
    b′′0::T
    e0::T
end
function BirchMurnaghan4th(v0, b0, b′0, b′′0, e0)
    T = Base.promote_typeof(v0, b0, b′0, b′′0, e0)
    return BirchMurnaghan4th{T}(convert.(T, [v0, b0, b′0, b′′0, e0])...)
end
BirchMurnaghan4th(v0::Real, b0::Real, b′0::Real, b′′0::Real) =
    BirchMurnaghan4th(v0, b0, b′0, b′′0, 0)
BirchMurnaghan4th(v0::AbstractQuantity, b0::AbstractQuantity, b′0, b′′0::AbstractQuantity) =
    BirchMurnaghan4th(v0, b0, b′0, b′′0, 0 * upreferred(Unitful.J))

"""
    PoirierTarantola2nd(v0, b0, e0)

Create a Poirier–Tarantola order equation of state. The elements' type will be handled automatically.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

See also: [`PoirierTarantola3rd`](@ref), [`PoirierTarantola4th`](@ref)

# Examples
```jldoctest
julia> PoirierTarantola2nd(1, 2.0)
PoirierTarantola2nd{Float64}(1.0, 2.0, 0.0)

julia> PoirierTarantola2nd(Int8(1), 2//1, 3.0)
PoirierTarantola2nd{Float64}(1.0, 2.0, 3.0)

julia> PoirierTarantola2nd(1u"nm^3", 2u"GPa", 3.0u"eV")
PoirierTarantola2nd{Quantity{Float64,D,U} where U where D}(1.0 nm³, 2.0 GPa, 3.0 eV)
```
"""
struct PoirierTarantola2nd{T} <: FiniteStrainEOS{T}
    v0::T
    b0::T
    e0::T
end
function PoirierTarantola2nd(v0, b0, e0)
    T = Base.promote_typeof(v0, b0, e0)
    return PoirierTarantola2nd{T}(convert.(T, [v0, b0, e0])...)
end
PoirierTarantola2nd(v0::Real, b0::Real) = PoirierTarantola2nd(v0, b0, 0)
PoirierTarantola2nd(v0::AbstractQuantity, b0::AbstractQuantity) =
    PoirierTarantola2nd(v0, b0, 0 * upreferred(Unitful.J))

"""
    PoirierTarantola3rd(v0, b0, b′0, e0)

Create a Poirier–Tarantola 3rd order equation of state. The elements' type will be handled automatically.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

See also: [`PoirierTarantola2nd`](@ref), [`PoirierTarantola4th`](@ref)

# Examples
```jldoctest
julia> PoirierTarantola3rd(1, 2.0, 3)
PoirierTarantola3rd{Float64}(1.0, 2.0, 3.0, 0.0)

julia> PoirierTarantola3rd(Int8(1), 2//1, 3.0, Float16(4))
PoirierTarantola3rd{Float64}(1.0, 2.0, 3.0, 4.0)

julia> PoirierTarantola3rd(1u"nm^3", 2u"GPa", 3, 4.0u"eV")
PoirierTarantola3rd{Quantity{Float64,D,U} where U where D}(1.0 nm³, 2.0 GPa, 3.0, 4.0 eV)
```
"""
struct PoirierTarantola3rd{T} <: FiniteStrainEOS{T}
    v0::T
    b0::T
    b′0::T
    e0::T
end
function PoirierTarantola3rd(v0, b0, b′0, e0)
    T = Base.promote_typeof(v0, b0, b′0, e0)
    return PoirierTarantola3rd{T}(convert.(T, [v0, b0, b′0, e0])...)
end
PoirierTarantola3rd(v0::Real, b0::Real, b′0::Real) = PoirierTarantola3rd(v0, b0, b′0, 0)
PoirierTarantola3rd(v0::AbstractQuantity, b0::AbstractQuantity, b′0) =
    PoirierTarantola3rd(v0, b0, b′0, 0 * upreferred(Unitful.J))

"""
    PoirierTarantola4th(v0, b0, b′0, b′′0, e0)

Create a Birch–Murnaghan 4th order equation of state. The elements' type will be handled automatically.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `b′′0`: the second-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

See also: [`PoirierTarantola2nd`](@ref), [`PoirierTarantola3rd`](@ref)

# Examples
```jldoctest
julia> PoirierTarantola4th(1, 2.0, 3, 4)
PoirierTarantola4th{Float64}(1.0, 2.0, 3.0, 4.0, 0.0)

julia> PoirierTarantola4th(Int8(1), 2//1, 3.0, Float16(4), 5)
PoirierTarantola4th{Float64}(1.0, 2.0, 3.0, 4.0, 5.0)

julia> PoirierTarantola4th(1u"nm^3", 2u"GPa", 3, 4u"1/GPa", 5.0u"eV")
PoirierTarantola4th{Quantity{Float64,D,U} where U where D}(1.0 nm³, 2.0 GPa, 3.0, 4.0 GPa⁻¹, 5.0 eV)
```
"""
struct PoirierTarantola4th{T} <: FiniteStrainEOS{T}
    v0::T
    b0::T
    b′0::T
    b′′0::T
    e0::T
end
function PoirierTarantola4th(v0, b0, b′0, b′′0, e0)
    T = Base.promote_typeof(v0, b0, b′0, b′′0, e0)
    return PoirierTarantola4th{T}(convert.(T, [v0, b0, b′0, b′′0, e0])...)
end
PoirierTarantola4th(v0::Real, b0::Real, b′0::Real, b′′0::Real) =
    PoirierTarantola4th(v0, b0, b′0, b′′0, 0)
PoirierTarantola4th(
    v0::AbstractQuantity,
    b0::AbstractQuantity,
    b′0,
    b′′0::AbstractQuantity,
) = PoirierTarantola4th(v0, b0, b′0, b′′0, 0 * upreferred(Unitful.J))

"""
    Vinet(v0, b0, b′0, e0)

Create a Vinet equation of state. The elements' type will be handled automatically.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

# Examples
```jldoctest
julia> Vinet(1, 2.0, 3)
Vinet{Float64}(1.0, 2.0, 3.0, 0.0)

julia> Vinet(Int8(1), 2//1, 3.0, Float16(4))
Vinet{Float64}(1.0, 2.0, 3.0, 4.0)

julia> Vinet(1u"nm^3", 2u"GPa", 3, 4.0u"eV")
Vinet{Quantity{Float64,D,U} where U where D}(1.0 nm³, 2.0 GPa, 3.0, 4.0 eV)
```
"""
struct Vinet{T} <: EquationOfState{T}
    v0::T
    b0::T
    b′0::T
    e0::T
end
function Vinet(v0, b0, b′0, e0)
    T = Base.promote_typeof(v0, b0, b′0, e0)
    return Vinet{T}(convert.(T, [v0, b0, b′0, e0])...)
end
Vinet(v0::Real, b0::Real, b′0::Real) = Vinet(v0, b0, b′0, 0)
Vinet(v0::AbstractQuantity, b0::AbstractQuantity, b′0) =
    Vinet(v0, b0, b′0, 0 * upreferred(Unitful.J))

struct AntonSchmidt{T} <: EquationOfState{T}
    v0::T
    β::T
    n::T
    e∞::T
end
function AntonSchmidt(v0, β, n, e∞)
    T = Base.promote_typeof(v0, β, n, e∞)
    return AntonSchmidt{T}(convert.(T, [v0, β, n, e∞])...)
end
AntonSchmidt(v0::Real, β::Real, n::Real) = AntonSchmidt(v0, β, n, 0)

struct BreenanStacey{T} <: EquationOfState{T}
    v0::T
    b0::T
    γ0::T
    e0::T
end
function BreenanStacey(v0, b0, γ0, e0)
    T = Base.promote_typeof(v0, b0, γ0, e0)
    return BreenanStacey{T}(convert.(T, [v0, b0, γ0, e0])...)
end
BreenanStacey(v0::Real, b0::Real, γ0::Real) = BreenanStacey(v0, b0, γ0, 0)

struct Shanker{T} <: EquationOfState{T}
    v0::T
    b0::T
    b′0::T
    e0::T
end
function Shanker(v0, b0, b′0, e0)
    T = Base.promote_typeof(v0, b0, b′0, e0)
    return Shanker{T}(convert.(T, [v0, b0, b′0, e0])...)
end
Shanker(v0::Real, b0::Real, b′0::Real) = Shanker(v0, b0, b′0, 0)
Shanker(v0::AbstractQuantity, b0::AbstractQuantity, b′0) =
    Shanker(v0, b0, b′0, 0 * upreferred(Unitful.J))
# =================================== Types ================================== #

# Energy evaluation
function _evaluate(eos::Murnaghan, ::Energy, v)
    v0, b0, b′0, e0 = fieldvalues(eos)
    x, y = b′0 - 1, (v0 / v)^b′0
    return e0 + b0 / b′0 * v * (y / x + 1) - v0 * b0 / x
end
function _evaluate(eos::BirchMurnaghan2nd, ::Energy, v)
    v0, b0, e0 = fieldvalues(eos)
    f = (cbrt(v0 / v)^2 - 1) / 2
    return e0 + 9 / 2 * b0 * v0 * f^2
end
function _evaluate(eos::BirchMurnaghan3rd, ::Energy, v)
    v0, b0, b′0, e0 = fieldvalues(eos)
    eta = cbrt(v0 / v)
    xi = eta^2 - 1
    return e0 + 9 / 16 * b0 * v0 * xi^2 * (6 + b′0 * xi - 4 * eta^2)
end
function _evaluate(eos::BirchMurnaghan4th, ::Energy, v)
    v0, b0, b′0, b′′0, e0 = fieldvalues(eos)
    f, h = (cbrt(v0 / v)^2 - 1) / 2, b0 * b′′0 + b′0^2
    return e0 + 3 / 8 * v0 * b0 * f^2 * ((9h - 63b′0 + 143) * f^2 + 12 * (b′0 - 4) * f + 12)
end
function _evaluate(eos::PoirierTarantola2nd, ::Energy, v)
    v0, b0, e0 = fieldvalues(eos)
    return e0 + b0 / 2 * v0 * cbrt(log(v / v0))^2
end
function _evaluate(eos::PoirierTarantola3rd, ::Energy, v)
    v0, b0, b′0, e0 = fieldvalues(eos)
    x = cbrt(v / v0)
    xi = -3 * log(x)
    return e0 + b0 / 6 * v0 * xi^2 * ((b′0 - 2) * xi + 3)
end
function _evaluate(eos::PoirierTarantola4th, ::Energy, v)
    v0, b0, b′0, b′′0, e0 = fieldvalues(eos)
    x = cbrt(v / v0)
    xi = log(x)
    h = b0 * b′′0 + b′0^2
    return e0 + b0 / 24v0 * xi^2 * ((h + 3b′0 + 3) * xi^2 + 4 * (b′0 + 2) * xi + 12)
end
function _evaluate(eos::Vinet, ::Energy, v)
    v0, b0, b′0, e0 = fieldvalues(eos)
    x, xi = cbrt(v / v0), 3 / 2 * (b′0 - 1)
    return e0 + 9b0 * v0 / xi^2 * (1 + (xi * (1 - x) - 1) * exp(xi * (1 - x)))
end
function _evaluate(eos::AntonSchmidt, ::Energy, v)
    v0, β, n, e∞ = fieldvalues(eos)
    x, η = v / v0, n + 1
    return e∞ + β * v0 / η * x^η * (log(x) - 1 / η)
end
# Pressure evaluation
function _evaluate(eos::Murnaghan, ::Pressure, v)
    v0, b0, b′0 = fieldvalues(eos)
    return b0 / b′0 * ((v0 / v)^b′0 - 1)
end
function _evaluate(eos::BirchMurnaghan2nd, ::Pressure, v)
    v0, b0 = fieldvalues(eos)
    f = (cbrt(v0 / v)^2 - 1) / 2
    return 3b0 * f * (1 + 2f)^(5 / 2)
end
function _evaluate(eos::BirchMurnaghan3rd, ::Pressure, v)
    v0, b0, b′0 = fieldvalues(eos)
    eta = cbrt(v0 / v)
    return 3 / 2 * b0 * (eta^7 - eta^5) * (1 + 3 / 4 * (b′0 - 4) * (eta^2 - 1))
end
function _evaluate(eos::BirchMurnaghan4th, ::Pressure, v)
    v0, b0, b′0, b′′0 = fieldvalues(eos)
    f, h = (cbrt(v0 / v)^2 - 1) / 2, b0 * b′′0 + b′0^2
    return b0 / 2 * (2f + 1)^(5 / 2) * ((9h - 63b′0 + 143) * f^2 + 9 * (b′0 - 4) * f + 6)
end
function _evaluate(eos::PoirierTarantola2nd, ::Pressure, v)
    v0, b0 = fieldvalues(eos)
    x = cbrt(v / v0)
    return -b0 / x * log(x)
end
function _evaluate(eos::PoirierTarantola3rd, ::Pressure, v)
    v0, b0, b′0 = fieldvalues(eos)
    x = v / v0
    xi = log(x)
    return -b0 * xi / 2x * ((b′0 - 2) * xi - 2)
end
function _evaluate(eos::PoirierTarantola4th, ::Pressure, v)
    v0, b0, b′0, b′′0 = fieldvalues(eos)
    x = cbrt(v / v0)
    xi = log(x)
    h = b0 * b′′0 + b′0^2
    return -b0 * xi / 6 / x * ((h + 3b′0 + 3) * xi^2 + 3 * (b′0 + 6) * xi + 6)
end
function _evaluate(eos::Vinet, ::Pressure, v)
    v0, b0, b′0 = fieldvalues(eos)
    x = cbrt(v / v0)
    xi = 3 / 2 * (b′0 - 1)
    return 3b0 / x^2 * (1 - x) * exp(xi * (1 - x))
end
function _evaluate(eos::AntonSchmidt, ::Pressure, v)
    v0, β, n = fieldvalues(eos)
    x = v / v0
    return -β * x^n * log(x)
end
function _evaluate(eos::BreenanStacey, ::Pressure, v)
    v0, b0, γ0 = fieldvalues(eos)
    x = v0 / v
    return b0 / 2 / γ0 * x^(4 / 3) * (exp(2γ0 * (1 - x)) - 1)
end
function _evaluate(eos::Shanker, ::Pressure, v)
    v0, b0, b′0 = fieldvalues(eos)
    x = v / v0
    y = 1 - x
    t = b′0 - 8 / 3
    return b0 / (x^(4 / 3) * t) *
           ((1 - 1 / t + 2 / t^2) * exp(t * y - 1) + y * (1 + y - 2 / t) * exp(t * y))
end
# Bulk modulus evaluation
function _evaluate(eos::BirchMurnaghan2nd, ::BulkModulus, v)
    v0, b0 = fieldvalues(eos)
    f = (cbrt(v0 / v)^2 - 1) / 2
    return b0 * (7f + 1) * (2f + 1)^(5 / 2)
end
function _evaluate(eos::BirchMurnaghan3rd, ::BulkModulus, v)
    v0, b0, b′0 = fieldvalues(eos)
    f = (cbrt(v0 / v)^2 - 1) / 2
    return b0 / 2 * (2f + 1)^(5 / 2) * ((27 * f^2 + 6f) * (b′0 - 4) - 4f + 2)
end
function _evaluate(eos::BirchMurnaghan4th, ::BulkModulus, v)
    v0, b0, b′0, b′′0 = fieldvalues(eos)
    f, h = (cbrt(v0 / v)^2 - 1) / 2, b0 * b′′0 + b′0^2
    return b0 / 6 *
           (2f + 1)^(5 / 2) *
           ((99h - 693b′0 + 1573) * f^3 + (27h - 108b′0 + 105) * f^2 + 6f * (3b′0 - 5) + 6)
end
function _evaluate(eos::PoirierTarantola2nd, ::BulkModulus, v)
    v0, b0 = fieldvalues(eos)
    x = cbrt(v / v0)
    return b0 / x * (1 - log(x))
end
function _evaluate(eos::PoirierTarantola3rd, ::BulkModulus, v)
    v0, b0, b′0 = fieldvalues(eos)
    x = v / v0
    xi = log(x)
    return -b0 / 2x * (((b′0 - 2) * xi + 2 - 2b′0) * xi + 2)
end
function _evaluate(eos::PoirierTarantola4th, ::BulkModulus, v)
    v0, b0, b′0, b′′0 = fieldvalues(eos)
    x = cbrt(v / v0)
    xi = log(x)
    h = b0 * b′′0 + b′0^2
    return -b0 / (6x) *
           ((h + 3b′0 + 3) * xi^3 - 3 * xi^2 * (h + 2b′0 + 1) - 6xi * (b′0 + 1) - 6)
end
function _evaluate(eos::Vinet, ::BulkModulus, v)
    v0, b0, b′0 = fieldvalues(eos)
    x, xi = cbrt(v / v0), 3 / 2 * (b′0 - 1)
    return -b0 / (2 * x^2) * (3x * (x - 1) * (b′0 - 1) + 2 * (x - 2)) * exp(-xi * (x - 1))
end
function _evaluate(eos::AntonSchmidt, ::BulkModulus, v)
    v0, β, n = fieldvalues(eos)
    x = v / v0
    return β * x^n * (1 + n * log(x))
end
function _evaluate(eos::Shanker, ::BulkModulus, v)
    v0, b0, b′0 = fieldvalues(eos)
    x = v / v0
    y = 1 - x
    t = b′0 - 8 / 3
    return b0 / cbrt(x) * (1 + y + y^2) * exp(t * y) + 4 / 3 * eos(Pressure())(v)
end

# Miscellaneous
if VERSION >= v"1.3"
    (eos::EquationOfState)(property::PhysicalProperty) = v -> _evaluate(eos, property, v)
else
    for T in (
        :Murnaghan,
        :BirchMurnaghan2nd,
        :BirchMurnaghan3rd,
        :BirchMurnaghan4th,
        :PoirierTarantola2nd,
        :PoirierTarantola3rd,
        :PoirierTarantola4th,
        :Vinet,
        :AntonSchmidt,
        :BreenanStacey,
        :Shanker,
    )
        eval(quote
            (eos::$T)(property::PhysicalProperty) = v -> _evaluate(eos, property, v)
        end)
    end  # Julia 1.0-1.2 does not support adding methods to abstract types.
end

Base.:(==)(x::T, y::T) where {T<:EquationOfState} = all(fieldvalues(x) .== fieldvalues(y))

Base.eltype(::FieldValues{<:EquationOfState{T}}) where {T} = T
Base.eltype(::Type{<:EquationOfState{T}}) where {T} = T

function Base.getproperty(eos::EquationOfState, name::Symbol)
    if name ∈ (:bp0, :bd0)
        return getfield(eos, :b′0)
    elseif name ∈ (:bpp0, :bdd0)
        return getfield(eos, :b′′0)
    else
        return getfield(eos, name)
    end
end

Unitful.upreferred(::typeof(dimension(u"J"))) = u"eV"
Unitful.upreferred(::typeof(dimension(u"m^3"))) = u"angstrom^3"
Unitful.upreferred(::typeof(dimension(u"Pa"))) = u"eV/angstrom^3"

end
