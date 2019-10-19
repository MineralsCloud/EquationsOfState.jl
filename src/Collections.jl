"""
This module provides `EquationOfState` types and `apply` methods to calculate
energy, pressure, or bulk modulus of an `EquationOfState` on
a (an) volume (array of volumes).
"""
module Collections

using Unitful
using Unitful: AbstractQuantity
import Base: ==

using EquationsOfState: EnergyForm, PressureForm, BulkModulusForm

export apply,
       EquationOfState,
       FiniteStrainEquationOfState,
       Murnaghan,
       BirchMurnaghan2nd,
       BirchMurnaghan3rd,
       BirchMurnaghan4th,
       PoirierTarantola2nd,
       PoirierTarantola3rd,
       PoirierTarantola4th,
       Vinet,
       AntonSchmidt,
       BreenanStacey

# ============================================================================ #
#                                     Types                                    #
# ============================================================================ #
"""
    EquationOfState{T}

An abstraction of equations of state, where `T` specifies the elements' common type.
"""
abstract type EquationOfState{T} end

"""
    FiniteStrainEquationOfState{T} <: EquationOfState{T}

An abstraction of finite strain equations of state, where `T` specifies the elements' common type.
"""
abstract type FiniteStrainEquationOfState{T} <: EquationOfState{T} end

"""
    Murnaghan(v0, b0, bp0, e0)

Create a Murnaghan equation of state. The elements' type will be handled automatically.

This equation of state can have units. The units are specified in [`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl)'s
`@u_str` style.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `bp0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

# Examples
```jldoctest
julia> Murnaghan(1, 2, 3.0)
Murnaghan{Float64}(1.0, 2.0, 3.0, 0.0)

julia> Murnaghan(Int8(1), 2//1, 3.0, 4)
Murnaghan{Float64}(1.0, 2.0, 3.0, 4.0)

julia> Murnaghan(1u"nm^3", 2u"GPa", 3, 3.0u"eV")
Murnaghan{Quantity{Float64,D,U} where U where D}(1.0 nm^3, 2.0 GPa, 3.0, 3.0 eV)
```
"""
struct Murnaghan{T} <: EquationOfState{T}
    v0::T
    b0::T
    bp0::T
    e0::T
end
function Murnaghan(v0, b0, bp0, e0)
    T = Base.promote_typeof(v0, b0, bp0, e0)
    return Murnaghan{T}(convert.(T, [v0, b0, bp0, e0])...)
end
Murnaghan(v0::Real, b0::Real, bp0::Real) = Murnaghan(v0, b0, bp0, 0)
Murnaghan(v0::AbstractQuantity, b0::AbstractQuantity, bp0) =
    Murnaghan(v0, b0, bp0, 0 * upreferred(Unitful.J))

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
BirchMurnaghan2nd{Quantity{Float64,D,U} where U where D}(1.0 nm^3, 2.0 GPa, 3.0 eV)
```
"""
struct BirchMurnaghan2nd{T} <: FiniteStrainEquationOfState{T}
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
    BirchMurnaghan3rd(v0, b0, bp0, e0)

Create a Birch–Murnaghan 3rd order equation of state. The elements' type will be handled automatically.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `bp0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

See also: [`BirchMurnaghan2nd`](@ref), [`BirchMurnaghan4th`](@ref)

# Examples
```jldoctest
julia> BirchMurnaghan3rd(1, 2.0, 3)
BirchMurnaghan3rd{Float64}(1.0, 2.0, 3.0, 0.0)

julia> BirchMurnaghan3rd(Int8(1), 2//1, 4, 0.0)
BirchMurnaghan3rd{Float64}(1.0, 2.0, 4.0, 0.0)

julia> BirchMurnaghan3rd(1u"nm^3", 2u"GPa", 4.0, 3u"eV")
BirchMurnaghan3rd{Quantity{Float64,D,U} where U where D}(1.0 nm^3, 2.0 GPa, 4.0, 3.0 eV)
```
"""
struct BirchMurnaghan3rd{T} <: FiniteStrainEquationOfState{T}
    v0::T
    b0::T
    bp0::T
    e0::T
end
function BirchMurnaghan3rd(v0, b0, bp0, e0)
    T = Base.promote_typeof(v0, b0, bp0, e0)
    return BirchMurnaghan3rd{T}(convert.(T, [v0, b0, bp0, e0])...)
end
BirchMurnaghan3rd(v0::Real, b0::Real, bp0::Real) = BirchMurnaghan3rd(v0, b0, bp0, 0)
BirchMurnaghan3rd(v0::AbstractQuantity, b0::AbstractQuantity, bp0) =
    BirchMurnaghan3rd(v0, b0, bp0, 0 * upreferred(Unitful.J))

"""
    BirchMurnaghan4th(v0, b0, bp0, bpp0, e0)

Create a Birch–Murnaghan 4th order equation of state. The elements' type will be handled automatically.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `bp0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `bpp0`: the second-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

See also: [`BirchMurnaghan2nd`](@ref), [`BirchMurnaghan4th`](@ref)

# Examples
```jldoctest
julia> BirchMurnaghan4th(1, 2.0, 3, 4)
BirchMurnaghan4th{Float64}(1.0, 2.0, 3.0, 4.0, 0.0)

julia> BirchMurnaghan4th(Int8(1), 2//1, 4, 5.0, Float16(6))
BirchMurnaghan4th{Float64}(1.0, 2.0, 4.0, 5.0, 6.0)

julia> BirchMurnaghan4th(1u"nm^3", 2u"GPa", 3.0, 4u"1/GPa", 5u"eV")
BirchMurnaghan4th{Quantity{Float64,D,U} where U where D}(1.0 nm^3, 2.0 GPa, 3.0, 4.0 GPa^-1, 5.0 eV)
```
"""
struct BirchMurnaghan4th{T} <: FiniteStrainEquationOfState{T}
    v0::T
    b0::T
    bp0::T
    bpp0::T
    e0::T
end
function BirchMurnaghan4th(v0, b0, bp0, bpp0, e0)
    T = Base.promote_typeof(v0, b0, bp0, bpp0, e0)
    return BirchMurnaghan4th{T}(convert.(T, [v0, b0, bp0, bpp0, e0])...)
end
BirchMurnaghan4th(v0::Real, b0::Real, bp0::Real, bpp0::Real) =
    BirchMurnaghan4th(v0, b0, bp0, bpp0, 0)
BirchMurnaghan4th(
    v0::AbstractQuantity,
    b0::AbstractQuantity,
    bp0,
    bpp0::AbstractQuantity,
) = BirchMurnaghan4th(v0, b0, bp0, bpp0, 0 * upreferred(Unitful.J))

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
PoirierTarantola2nd{Quantity{Float64,D,U} where U where D}(1.0 nm^3, 2.0 GPa, 3.0 eV)
```
"""
struct PoirierTarantola2nd{T} <: FiniteStrainEquationOfState{T}
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
    PoirierTarantola3rd(v0, b0, bp0, e0)

Create a Poirier–Tarantola 3rd order equation of state. The elements' type will be handled automatically.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `bp0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

See also: [`PoirierTarantola2nd`](@ref), [`PoirierTarantola4th`](@ref)

# Examples
```jldoctest
julia> PoirierTarantola3rd(1, 2.0, 3)
PoirierTarantola3rd{Float64}(1.0, 2.0, 3.0, 0.0)

julia> PoirierTarantola3rd(Int8(1), 2//1, 3.0, Float16(4))
PoirierTarantola3rd{Float64}(1.0, 2.0, 3.0, 4.0)

julia> PoirierTarantola3rd(1u"nm^3", 2u"GPa", 3, 4.0u"eV")
PoirierTarantola3rd{Quantity{Float64,D,U} where U where D}(1.0 nm^3, 2.0 GPa, 3.0, 4.0 eV)
```
"""
struct PoirierTarantola3rd{T} <: FiniteStrainEquationOfState{T}
    v0::T
    b0::T
    bp0::T
    e0::T
end
function PoirierTarantola3rd(v0, b0, bp0, e0)
    T = Base.promote_typeof(v0, b0, bp0, e0)
    return PoirierTarantola3rd{T}(convert.(T, [v0, b0, bp0, e0])...)
end
PoirierTarantola3rd(v0::Real, b0::Real, bp0::Real) = PoirierTarantola3rd(v0, b0, bp0, 0)
PoirierTarantola3rd(v0::AbstractQuantity, b0::AbstractQuantity, bp0) =
    PoirierTarantola3rd(v0, b0, bp0, 0 * upreferred(Unitful.J))

"""
    PoirierTarantola4th(v0, b0, bp0, bpp0, e0)

Create a Birch–Murnaghan 4th order equation of state. The elements' type will be handled automatically.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `bp0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `bpp0`: the second-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

See also: [`PoirierTarantola2nd`](@ref), [`PoirierTarantola3rd`](@ref)

# Examples
```jldoctest
julia> PoirierTarantola4th(1, 2.0, 3, 4)
PoirierTarantola4th{Float64}(1.0, 2.0, 3.0, 4.0, 0.0)

julia> PoirierTarantola4th(Int8(1), 2//1, 3.0, Float16(4), 5)
PoirierTarantola4th{Float64}(1.0, 2.0, 3.0, 4.0, 5.0)

julia> PoirierTarantola4th(1u"nm^3", 2u"GPa", 3, 4u"1/GPa", 5.0u"eV")
PoirierTarantola4th{Quantity{Float64,D,U} where U where D}(1.0 nm^3, 2.0 GPa, 3.0, 4.0 GPa^-1, 5.0 eV)
```
"""
struct PoirierTarantola4th{T} <: FiniteStrainEquationOfState{T}
    v0::T
    b0::T
    bp0::T
    bpp0::T
    e0::T
end
function PoirierTarantola4th(v0, b0, bp0, bpp0, e0)
    T = Base.promote_typeof(v0, b0, bp0, bpp0, e0)
    return PoirierTarantola4th{T}(convert.(T, [v0, b0, bp0, bpp0, e0])...)
end
PoirierTarantola4th(v0::Real, b0::Real, bp0::Real, bpp0::Real) =
    PoirierTarantola4th(v0, b0, bp0, bpp0, 0)
PoirierTarantola4th(
    v0::AbstractQuantity,
    b0::AbstractQuantity,
    bp0,
    bpp0::AbstractQuantity,
) = PoirierTarantola4th(v0, b0, bp0, bpp0, 0 * upreferred(Unitful.J))

"""
    Vinet(v0, b0, bp0, e0)

Create a Vinet equation of state. The elements' type will be handled automatically.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `bp0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

# Examples
```jldoctest
julia> Vinet(1, 2.0, 3)
Vinet{Float64}(1.0, 2.0, 3.0, 0.0)

julia> Vinet(Int8(1), 2//1, 3.0, Float16(4))
Vinet{Float64}(1.0, 2.0, 3.0, 4.0)

julia> Vinet(1u"nm^3", 2u"GPa", 3, 4.0u"eV")
Vinet{Quantity{Float64,D,U} where U where D}(1.0 nm^3, 2.0 GPa, 3.0, 4.0 eV)
```
"""
struct Vinet{T} <: EquationOfState{T}
    v0::T
    b0::T
    bp0::T
    e0::T
end
function Vinet(v0, b0, bp0, e0)
    T = Base.promote_typeof(v0, b0, bp0, e0)
    return Vinet{T}(convert.(T, [v0, b0, bp0, e0])...)
end
Vinet(v0::Real, b0::Real, bp0::Real) = Vinet(v0, b0, bp0, 0)
Vinet(v0::AbstractQuantity, b0::AbstractQuantity, bp0) =
    Vinet(v0, b0, bp0, 0 * upreferred(Unitful.J))

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
# =================================== Types ================================== #


# ============================================================================ #
#                               Energy evaluation                              #
# ============================================================================ #
"""
    apply(EnergyForm(), eos::EquationOfState)
    apply(PressureForm(), eos::EquationOfState)
    apply(BulkModulusForm(), eos::EquationOfState)

Return a function that takes a volume as a variable, suitable for mapping onto an array.

# Examples
```jldoctest
julia> using EquationsOfState, EquationsOfState.Collections

julia> f = apply(EnergyForm(), Vinet(1, 2, 3));

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

julia> g = apply(PressureForm(), Vinet(1, 2, 3));

julia> map(g, 1:1:10)
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

julia> h = apply(BulkModulusForm(), BirchMurnaghan3rd(1, 2, 3));

julia> map(h, 1:1:10)
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
apply(form::EnergyForm, eos::EquationOfState) = v -> apply(form, eos, v)
"""
    apply(EnergyForm(), eos::EquationOfState, v)

Return the energy of an `EquationOfState` on volume `v`. If `eos` has units,
`v` must also has.
"""
function apply(::EnergyForm, eos::Murnaghan, v)
    v0, b0, bp0, e0 = fieldvalues(eos)

    x = bp0 - 1
    y = (v0 / v)^bp0
    return e0 + b0 / bp0 * v * (y / x + 1) - v0 * b0 / x
end
function apply(::EnergyForm, eos::BirchMurnaghan2nd, v)
    v0, b0, e0 = fieldvalues(eos)

    f = (cbrt(v0 / v)^2 - 1) / 2
    return e0 + 9 / 2 * b0 * v0 * f^2
end
function apply(::EnergyForm, eos::BirchMurnaghan3rd, v)
    v0, b0, bp0, e0 = fieldvalues(eos)

    eta = cbrt(v0 / v)
    xi = eta^2 - 1
    return e0 + 9 / 16 * b0 * v0 * xi^2 * (6 + bp0 * xi - 4 * eta^2)
end
function apply(::EnergyForm, eos::BirchMurnaghan4th, v)
    v0, b0, bp0, bpp0, e0 = fieldvalues(eos)

    f = (cbrt(v0 / v)^2 - 1) / 2
    h = b0 * bpp0 + bp0^2
    return e0 + 3 / 8 * v0 * b0 * f^2 * ((9h - 63bp0 + 143) * f^2 + 12 * (bp0 - 4) * f + 12)
end
function apply(::EnergyForm, eos::PoirierTarantola2nd, v)
    v0, b0, e0 = fieldvalues(eos)

    return e0 + b0 / 2 * v0 * log(v / v0)^(2 / 3)
end
function apply(::EnergyForm, eos::PoirierTarantola3rd, v)
    v0, b0, bp0, e0 = fieldvalues(eos)

    x = cbrt(v / v0)
    xi = -3 * log(x)
    return e0 + b0 / 6 * v0 * xi^2 * ((bp0 - 2) * xi + 3)
end
function apply(::EnergyForm, eos::PoirierTarantola4th, v)
    v0, b0, bp0, bpp0, e0 = fieldvalues(eos)

    x = cbrt(v / v0)
    xi = log(x)
    h = b0 * bpp0 + bp0^2
    return e0 + b0 / 24v0 * xi^2 * ((h + 3bp0 + 3) * xi^2 + 4 * (bp0 + 2) * xi + 12)
end
function apply(::EnergyForm, eos::Vinet, v)
    v0, b0, bp0, e0 = fieldvalues(eos)

    x = cbrt(v / v0)
    xi = 3 / 2 * (bp0 - 1)
    return e0 + 9b0 * v0 / xi^2 * (1 + (xi * (1 - x) - 1) * exp(xi * (1 - x)))
end
function apply(::EnergyForm, eos::AntonSchmidt, v)
    v0, β, n, e∞ = fieldvalues(eos)

    x = v / v0
    η = n + 1
    return e∞ + β * v0 / η * x^η * (log(x) - 1 / η)
end
# ============================= Energy evaluation ============================ #


# ============================================================================ #
#                              Pressure evaluation                             #
# ============================================================================ #
apply(::PressureForm, eos::EquationOfState) = v -> apply(PressureForm(), eos, v)
"""
    apply(PressureForm(), eos::EquationOfState, v)

Return the pressure of an `EquationOfState` on volume `v`. If `eos` has units,
`v` must also has.
"""
function apply(::PressureForm, eos::Murnaghan, v)
    v0, b0, bp0 = fieldvalues(eos)

    return b0 / bp0 * ((v0 / v)^bp0 - 1)
end
function apply(::PressureForm, eos::BirchMurnaghan2nd, v)
    v0, b0 = fieldvalues(eos)

    f = ((v0 / v)^(2 / 3) - 1) / 2
    return 3b0 * f * (1 + 2f)^(5 / 2)
end
function apply(::PressureForm, eos::BirchMurnaghan3rd, v)
    v0, b0, bp0 = fieldvalues(eos)

    eta = (v0 / v)^(1 / 3)
    return 3 / 2 * b0 * (eta^7 - eta^5) * (1 + 3 / 4 * (bp0 - 4) * (eta^2 - 1))
end
function apply(::PressureForm, eos::BirchMurnaghan4th, v)
    v0, b0, bp0, bpp0 = fieldvalues(eos)

    f = ((v0 / v)^(2 / 3) - 1) / 2
    h = b0 * bpp0 + bp0^2
    return b0 / 2 * (2f + 1)^(5 / 2) * ((9h - 63bp0 + 143) * f^2 + 9 * (bp0 - 4) * f + 6)
end
function apply(::PressureForm, eos::PoirierTarantola2nd, v)
    v0, b0 = fieldvalues(eos)

    x = (v / v0)^(1 / 3)
    return -b0 / x * log(x)
end
function apply(::PressureForm, eos::PoirierTarantola3rd, v)
    v0, b0, bp0 = fieldvalues(eos)

    x = v / v0
    xi = log(x)
    return -b0 * xi / 2x * ((bp0 - 2) * xi - 2)
end
function apply(::PressureForm, eos::PoirierTarantola4th, v)
    v0, b0, bp0, bpp0 = fieldvalues(eos)

    x = (v / v0)^(1 / 3)
    xi = log(x)
    h = b0 * bpp0 + bp0^2
    return -b0 * xi / 6 / x * ((h + 3bp0 + 3) * xi^2 + 3 * (bp0 + 6) * xi + 6)
end
function apply(::PressureForm, eos::Vinet, v)
    v0, b0, bp0 = fieldvalues(eos)

    x = (v / v0)^(1 / 3)
    xi = 3 / 2 * (bp0 - 1)
    return 3b0 / x^2 * (1 - x) * exp(xi * (1 - x))
end
function apply(::PressureForm, eos::AntonSchmidt, v)
    v0, β, n = fieldvalues(eos)

    x = v / v0
    return -β * x^n * log(x)
end
function apply(::PressureForm, eos::BreenanStacey, v)
    v0, b0, γ0 = fieldvalues(eos)

    x = v0 / v
    return b0 / 2 / γ0 * x^(4 / 3) * (exp(2γ0 * (1 - x)) - 1)
end
# ============================ Pressure evaluation =========================== #


# ============================================================================ #
#                            Bulk modulus evaluation                           #
# ============================================================================ #
apply(::BulkModulusForm, eos::EquationOfState) = v -> apply(BulkModulusForm(), eos, v)
"""
    apply(BulkModulusForm(), eos::EquationOfState, v)

Return the bulk modulus of an `EquationOfState` on volume `v`. If `eos` has units,
`v` must also has.
"""
function apply(::BulkModulusForm, eos::BirchMurnaghan2nd, v)
    v0, b0 = fieldvalues(eos)

    f = ((v0 / v)^(2 / 3) - 1) / 2
    return b0 * (7f + 1) * (2f + 1)^(5 / 2)
end
function apply(::BulkModulusForm, eos::BirchMurnaghan3rd, v)
    v0, b0, bp0 = fieldvalues(eos)

    f = ((v0 / v)^(2 / 3) - 1) / 2
    return b0 / 2 * (2f + 1)^(5 / 2) * ((27 * f^2 + 6f) * (bp0 - 4) - 4f + 2)
end
function apply(::BulkModulusForm, eos::BirchMurnaghan4th, v)
    v0, b0, bp0, bpp0 = fieldvalues(eos)

    f = ((v0 / v)^(2 / 3) - 1) / 2
    h = b0 * bpp0 + bp0^2
    return b0 / 6 * (2f + 1)^(5 / 2) *
           ((99h - 693bp0 + 1573) * f^3 + (27h - 108bp0 + 105) * f^2 + 6f * (3bp0 - 5) + 6)
end
function apply(::BulkModulusForm, eos::PoirierTarantola2nd, v)
    v0, b0 = fieldvalues(eos)

    x = (v / v0)^(1 / 3)
    return b0 / x * (1 - log(x))
end
function apply(::BulkModulusForm, eos::PoirierTarantola3rd, v)
    v0, b0, bp0 = fieldvalues(eos)

    x = v / v0
    xi = log(x)
    return -b0 / 2x * (((bp0 - 2) * xi + 2 - 2bp0) * xi + 2)
end
function apply(::BulkModulusForm, eos::PoirierTarantola4th, v)
    v0, b0, bp0, bpp0 = fieldvalues(eos)

    x = (v / v0)^(1 / 3)
    xi = log(x)
    h = b0 * bpp0 + bp0^2
    return -b0 / (6x) *
           ((h + 3bp0 + 3) * xi^3 - 3 * xi^2 * (h + 2bp0 + 1) - 6xi * (bp0 + 1) - 6)
end
function apply(::BulkModulusForm, eos::Vinet, v)
    v0, b0, bp0 = fieldvalues(eos)

    x = (v / v0)^(1 / 3)
    xi = 3 / 2 * (bp0 - 1)
    return -b0 / (2 * x^2) * (3x * (x - 1) * (bp0 - 1) + 2 * (x - 2)) * exp(-xi * (x - 1))
end
function apply(::BulkModulusForm, eos::AntonSchmidt, v)
    v0, β, n = fieldvalues(eos)

    x = v / v0
    return β * x^n * (1 + n * log(x))
end
# ========================== Bulk modulus evaluation ========================= #


# ============================================================================ #
#                                 Miscellaneous                                #
# ============================================================================ #
# This is a helper function and should not be exported.
fieldvalues(eos::EquationOfState) = [getfield(eos, i) for i in 1:nfields(eos)]

function ==(x::T, y::T) where {T<:EquationOfState}
    return all(getfield(x, i) == getfield(y, i) for i in 1:fieldcount(T))
end

Unitful.upreferred(::typeof(dimension(u"J"))) = u"eV"
Unitful.upreferred(::typeof(dimension(u"m^3"))) = u"angstrom^3"
Unitful.upreferred(::typeof(dimension(u"Pa"))) = u"eV/angstrom^3"
# =============================== Miscellaneous ============================== #

end
