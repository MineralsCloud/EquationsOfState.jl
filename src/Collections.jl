"""
This module provides `EOSParameters` types and calculate
energy, pressure, or bulk modulus of an `EOSParameters` on
a (an) volume (array of volumes).
"""
module Collections

using AutoHashEquals: @auto_hash_equals
using IterTools: FieldValues, fieldvalues
using LinearAlgebra: dot
using StaticArrays: SVector
using Unitful: AbstractQuantity, @u_str

import Unitful

export Murnaghan,
    BirchMurnaghan2nd,
    BirchMurnaghan3rd,
    BirchMurnaghan4th,
    BirchMurnaghan5th,
    PoirierTarantola2nd,
    PoirierTarantola3rd,
    PoirierTarantola4th,
    PoirierTarantola5th,
    Vinet,
    AntonSchmidt,
    BreenanStacey,
    Shanker,
    PolynomialEOS,
    EnergyEquation,
    PressureEquation,
    BulkModulusEquation

"""
    EOSParameters{T}

An abstraction of equations of state, where `T` specifies the elements' common type.
"""
abstract type EOSParameters{T} end

"""
    FiniteStrainEOS{T} <: EOSParameters{T}

An abstraction of finite strain equations of state, where `T` specifies the elements' common type.
"""
abstract type FiniteStrainEOS{T} <: EOSParameters{T} end

"""
    Murnaghan(v0, b0, b′0, e0)

Create a Murnaghan equation of state.

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
@auto_hash_equals struct Murnaghan{T} <: EOSParameters{T}
    v0::T
    b0::T
    b′0::T
    e0::T
end
function Murnaghan(v0, b0, b′0, e0)
    T = Base.promote_typeof(v0, b0, b′0, e0)
    return Murnaghan{T}((convert(T, x) for x in (v0, b0, b′0, e0))...)  # Cannot use `T.(args...)`! For `AbstractQuantity` they will fail!
end
Murnaghan(v0::Real, b0::Real, b′0::Real) = Murnaghan(v0, b0, b′0, 0)
Murnaghan(v0::AbstractQuantity, b0::AbstractQuantity, b′0) =
    Murnaghan(v0, b0, b′0, 0 * u"eV")

"""
    BirchMurnaghan2nd(v0, b0, e0)

Create a Birch–Murnaghan 2nd order equation of state.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

See also: [`BirchMurnaghan3rd`](@ref), [`BirchMurnaghan4th`](@ref), [`BirchMurnaghan5th`](@ref)

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
@auto_hash_equals struct BirchMurnaghan2nd{T} <: FiniteStrainEOS{T}
    v0::T
    b0::T
    e0::T
end
function BirchMurnaghan2nd(v0, b0, e0)
    T = Base.promote_typeof(v0, b0, e0)
    return BirchMurnaghan2nd{T}((convert(T, x) for x in (v0, b0, e0))...)
end
BirchMurnaghan2nd(v0::Real, b0::Real) = BirchMurnaghan2nd(v0, b0, 0)
BirchMurnaghan2nd(v0::AbstractQuantity, b0::AbstractQuantity) =
    BirchMurnaghan2nd(v0, b0, 0 * u"eV")

"""
    BirchMurnaghan3rd(v0, b0, b′0, e0)

Create a Birch–Murnaghan 3rd order equation of state.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

See also: [`BirchMurnaghan2nd`](@ref), [`BirchMurnaghan4th`](@ref), [`BirchMurnaghan5th`](@ref)

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
@auto_hash_equals struct BirchMurnaghan3rd{T} <: FiniteStrainEOS{T}
    v0::T
    b0::T
    b′0::T
    e0::T
end
function BirchMurnaghan3rd(v0, b0, b′0, e0)
    T = Base.promote_typeof(v0, b0, b′0, e0)
    return BirchMurnaghan3rd{T}((convert(T, x) for x in (v0, b0, b′0, e0))...)
end
BirchMurnaghan3rd(v0::Real, b0::Real, b′0::Real) = BirchMurnaghan3rd(v0, b0, b′0, 0)
BirchMurnaghan3rd(v0::AbstractQuantity, b0::AbstractQuantity, b′0) =
    BirchMurnaghan3rd(v0, b0, b′0, 0 * u"eV")

"""
    BirchMurnaghan4th(v0, b0, b′0, b′′0, e0)

Create a Birch–Murnaghan 4th order equation of state.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `b′′0`: the second-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

See also: [`BirchMurnaghan2nd`](@ref), [`BirchMurnaghan4th`](@ref), [`BirchMurnaghan5th`](@ref)

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
@auto_hash_equals struct BirchMurnaghan4th{T} <: FiniteStrainEOS{T}
    v0::T
    b0::T
    b′0::T
    b′′0::T
    e0::T
end
function BirchMurnaghan4th(v0, b0, b′0, b′′0, e0)
    T = Base.promote_typeof(v0, b0, b′0, b′′0, e0)
    return BirchMurnaghan4th{T}((convert(T, x) for x in (v0, b0, b′0, b′′0, e0))...)
end
BirchMurnaghan4th(v0::Real, b0::Real, b′0::Real, b′′0::Real) =
    BirchMurnaghan4th(v0, b0, b′0, b′′0, 0)
BirchMurnaghan4th(v0::AbstractQuantity, b0::AbstractQuantity, b′0, b′′0::AbstractQuantity) =
    BirchMurnaghan4th(v0, b0, b′0, b′′0, 0 * u"eV")

"""
    BirchMurnaghan5th(v0, b0, b′0, b′′0, b′′′0, e0)

Create a Birch–Murnaghan 5th order equation of state.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `b′′0`: the second-order pressure-derivative bulk modulus of solid at zero pressure.
- `b′′′0`: the third-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

See also: [`BirchMurnaghan2nd`](@ref), [`BirchMurnaghan3rd`](@ref), [`BirchMurnaghan4th`](@ref)

# Examples
```jldoctest
julia> BirchMurnaghan5th(1, 2.0, 3, 4, 5 // 1)
BirchMurnaghan5th{Float64}(1.0, 2.0, 3.0, 4.0, 5.0, 0.0)

julia> BirchMurnaghan5th(Int8(1), 2//1, 3.0, Float16(4), 5)
BirchMurnaghan5th{Float64}(1.0, 2.0, 3.0, 4.0, 5.0, 0.0)

julia> BirchMurnaghan5th(1u"nm^3", 2u"GPa", 3, 4u"1/GPa", 5u"1/GPa^2", 6.0u"eV")
BirchMurnaghan5th{Quantity{Float64,D,U} where U where D}(1.0 nm³, 2.0 GPa, 3.0, 4.0 GPa⁻¹, 5.0 GPa⁻², 6.0 eV)
```
"""
@auto_hash_equals struct BirchMurnaghan5th{T} <: FiniteStrainEOS{T}
    v0::T
    b0::T
    b′0::T
    b′′0::T
    b′′′0::T
    e0::T
end
function BirchMurnaghan5th(v0, b0, b′0, b′′0, b′′′0, e0)
    T = Base.promote_typeof(v0, b0, b′0, b′′0, b′′′0, e0)
    return BirchMurnaghan5th{T}((convert(T, x) for x in (v0, b0, b′0, b′′0, b′′′0, e0))...)
end
BirchMurnaghan5th(v0::Real, b0::Real, b′0::Real, b′′0::Real, b′′′0::Real) =
    BirchMurnaghan5th(v0, b0, b′0, b′′0, b′′′0, 0)
BirchMurnaghan5th(
    v0::AbstractQuantity,
    b0::AbstractQuantity,
    b′0,
    b′′0::AbstractQuantity,
    b′′′0::AbstractQuantity,
) = BirchMurnaghan5th(v0, b0, b′0, b′′0, b′′′0, 0 * u"eV")

"""
    PoirierTarantola2nd(v0, b0, e0)

Create a Poirier–Tarantola order equation of state.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

See also: [`PoirierTarantola3rd`](@ref), [`PoirierTarantola4th`](@ref), [`PoirierTarantola5th`](@ref)

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
@auto_hash_equals struct PoirierTarantola2nd{T} <: FiniteStrainEOS{T}
    v0::T
    b0::T
    e0::T
end
function PoirierTarantola2nd(v0, b0, e0)
    T = Base.promote_typeof(v0, b0, e0)
    return PoirierTarantola2nd{T}((convert(T, x) for x in (v0, b0, e0))...)
end
PoirierTarantola2nd(v0::Real, b0::Real) = PoirierTarantola2nd(v0, b0, 0)
PoirierTarantola2nd(v0::AbstractQuantity, b0::AbstractQuantity) =
    PoirierTarantola2nd(v0, b0, 0 * u"eV")

"""
    PoirierTarantola3rd(v0, b0, b′0, e0)

Create a Poirier–Tarantola 3rd order equation of state.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

See also: [`PoirierTarantola2nd`](@ref), [`PoirierTarantola4th`](@ref), [`PoirierTarantola5th`](@ref)

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
@auto_hash_equals struct PoirierTarantola3rd{T} <: FiniteStrainEOS{T}
    v0::T
    b0::T
    b′0::T
    e0::T
end
function PoirierTarantola3rd(v0, b0, b′0, e0)
    T = Base.promote_typeof(v0, b0, b′0, e0)
    return PoirierTarantola3rd{T}((convert(T, x) for x in (v0, b0, b′0, e0))...)
end
PoirierTarantola3rd(v0::Real, b0::Real, b′0::Real) = PoirierTarantola3rd(v0, b0, b′0, 0)
PoirierTarantola3rd(v0::AbstractQuantity, b0::AbstractQuantity, b′0) =
    PoirierTarantola3rd(v0, b0, b′0, 0 * u"eV")

"""
    PoirierTarantola4th(v0, b0, b′0, b′′0, e0)

Create a Poirier–Tarantola 4th order equation of state.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `b′′0`: the second-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

See also: [`PoirierTarantola2nd`](@ref), [`PoirierTarantola3rd`](@ref), [`PoirierTarantola5th`](@ref)

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
@auto_hash_equals struct PoirierTarantola4th{T} <: FiniteStrainEOS{T}
    v0::T
    b0::T
    b′0::T
    b′′0::T
    e0::T
end
function PoirierTarantola4th(v0, b0, b′0, b′′0, e0)
    T = Base.promote_typeof(v0, b0, b′0, b′′0, e0)
    return PoirierTarantola4th{T}((convert(T, x) for x in (v0, b0, b′0, b′′0, e0))...)
end
PoirierTarantola4th(v0::Real, b0::Real, b′0::Real, b′′0::Real) =
    PoirierTarantola4th(v0, b0, b′0, b′′0, 0)
PoirierTarantola4th(
    v0::AbstractQuantity,
    b0::AbstractQuantity,
    b′0,
    b′′0::AbstractQuantity,
) = PoirierTarantola4th(v0, b0, b′0, b′′0, 0 * u"eV")

"""
    PoirierTarantola5th(v0, b0, b′0, b′′0, b′′′0, e0)

Create a Poirier–Tarantola 5th order equation of state.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `b′0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `b′′0`: the second-order pressure-derivative bulk modulus of solid at zero pressure.
- `b′′′0`: the third-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0`: the energy of solid at zero pressure. Its default value is `0u"eV"` (`0`), if other parameters have (no) units.

See also: [`PoirierTarantola2nd`](@ref), [`PoirierTarantola3rd`](@ref), [`PoirierTarantola4th`](@ref)

# Examples
```jldoctest
julia> PoirierTarantola5th(1, 2.0, 3, 4, 5 // 1)
PoirierTarantola5th{Float64}(1.0, 2.0, 3.0, 4.0, 5.0, 0.0)

julia> PoirierTarantola5th(Int8(1), 2//1, 3.0, Float16(4), 5)
PoirierTarantola5th{Float64}(1.0, 2.0, 3.0, 4.0, 5.0, 0.0)

julia> PoirierTarantola5th(1u"nm^3", 2u"GPa", 3, 4u"1/GPa", 5u"1/GPa^2", 6.0u"eV")
PoirierTarantola5th{Quantity{Float64,D,U} where U where D}(1.0 nm³, 2.0 GPa, 3.0, 4.0 GPa⁻¹, 5.0 GPa⁻², 6.0 eV)
```
"""
@auto_hash_equals struct PoirierTarantola5th{T} <: FiniteStrainEOS{T}
    v0::T
    b0::T
    b′0::T
    b′′0::T
    b′′′0::T
    e0::T
end
function PoirierTarantola5th(v0, b0, b′0, b′′0, b′′′0, e0)
    T = Base.promote_typeof(v0, b0, b′0, b′′0, b′′′0, e0)
    return PoirierTarantola5th{T}((
        convert(T, x) for x in (v0, b0, b′0, b′′0, b′′′0, e0)
    )...)
end
PoirierTarantola5th(v0::Real, b0::Real, b′0::Real, b′′0::Real, b′′′0::Real) =
    PoirierTarantola5th(v0, b0, b′0, b′′0, b′′′0, 0)
PoirierTarantola5th(
    v0::AbstractQuantity,
    b0::AbstractQuantity,
    b′0,
    b′′0::AbstractQuantity,
    b′′′0::AbstractQuantity,
) = PoirierTarantola5th(v0, b0, b′0, b′′0, b′′′0, 0 * u"eV")

"""
    Vinet(v0, b0, b′0, e0)

Create a Vinet equation of state.

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
@auto_hash_equals struct Vinet{T} <: EOSParameters{T}
    v0::T
    b0::T
    b′0::T
    e0::T
end
function Vinet(v0, b0, b′0, e0)
    T = Base.promote_typeof(v0, b0, b′0, e0)
    return Vinet{T}((convert(T, x) for x in (v0, b0, b′0, e0))...)
end
Vinet(v0::Real, b0::Real, b′0::Real) = Vinet(v0, b0, b′0, 0)
Vinet(v0::AbstractQuantity, b0::AbstractQuantity, b′0) = Vinet(v0, b0, b′0, 0 * u"eV")

@auto_hash_equals struct AntonSchmidt{T} <: EOSParameters{T}
    v0::T
    β::T
    n::T
    e∞::T
end
function AntonSchmidt(v0, β, n, e∞)
    T = Base.promote_typeof(v0, β, n, e∞)
    return AntonSchmidt{T}((convert(T, x) for x in (v0, β, n, e∞))...)
end
AntonSchmidt(v0::Real, β::Real, n::Real) = AntonSchmidt(v0, β, n, 0)

@auto_hash_equals struct BreenanStacey{T} <: EOSParameters{T}
    v0::T
    b0::T
    γ0::T
    e0::T
end
function BreenanStacey(v0, b0, γ0, e0)
    T = Base.promote_typeof(v0, b0, γ0, e0)
    return BreenanStacey{T}((convert(T, x) for x in (v0, b0, γ0, e0))...)
end
BreenanStacey(v0::Real, b0::Real, γ0::Real) = BreenanStacey(v0, b0, γ0, 0)

@auto_hash_equals struct Shanker{T} <: EOSParameters{T}
    v0::T
    b0::T
    b′0::T
    e0::T
end
function Shanker(v0, b0, b′0, e0)
    T = Base.promote_typeof(v0, b0, b′0, e0)
    return Shanker{T}((convert(T, x) for x in (v0, b0, b′0, e0))...)
end
Shanker(v0::Real, b0::Real, b′0::Real) = Shanker(v0, b0, b′0, 0)
Shanker(v0::AbstractQuantity, b0::AbstractQuantity, b′0) = Shanker(v0, b0, b′0, 0 * u"eV")

@auto_hash_equals struct PolynomialEOS{N,T} <: EOSParameters{T}
    v0::T
    p0::SVector{N,T}
    e0::T
end
function PolynomialEOS(v0, p0::AbstractVector, e0)
    T = Base.promote_typeof(v0, p0..., e0)
    return PolynomialEOS{length(p0),T}(
        convert(T, v0),
        [convert(T, x) for x in p0],
        convert(T, e0),
    )
end
PolynomialEOS(v0::Real, p0::AbstractVector{<:Real}) = PolynomialEOS(v0, p0, 0)
PolynomialEOS(v0::AbstractQuantity, p0::AbstractVector{<:AbstractQuantity}) =
    PolynomialEOS(v0, p0, 0 * u"eV")

abstract type EquationOfState{T<:EOSParameters} end
"""
    EnergyEquation()
    (::EOSParameters)(::EnergyEquation)(v)
    (::EOSParameters)(::EnergyEquation)

Return the energy of an `EOSParameters` on volume `v`. If `eos` has units,
`v` must also has.

Return a [function-like object](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects-1) that takes a volume as a variable, suitable for mapping onto an array.

# Examples
```jldoctest
julia> f = Vinet(1, 2, 3)(EnergyEquation());

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
           Vinet(1, 2, 3)(EnergyEquation())(v)
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
struct EnergyEquation{T} <: EquationOfState{T}
    params::T
end
"""
    PressureEquation()
    (::EOSParameters)(::PressureEquation)(v)
    (::EOSParameters)(::PressureEquation)

Return the pressure of an `EOSParameters` on volume `v`. If `eos` has units,
`v` must also has.

# Examples
```jldoctest
julia> f = Vinet(1, 2, 3)(PressureEquation());

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
struct PressureEquation{T} <: EquationOfState{T}
    params::T
end
"""
    BulkModulusEquation()
    (::EOSParameters)(::BulkModulusEquation)(v)
    (::EOSParameters)(::BulkModulusEquation)

Return the bulk modulus of an `EOSParameters` on volume `v`. If `eos` has units,
`v` must also has.

# Examples
```jldoctest
julia> f = BirchMurnaghan3rd(1, 2, 3)(BulkModulusEquation());

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
struct BulkModulusEquation{T} <: EquationOfState{T}
    params::T
end

# EnergyEquation evaluation
function (f::EnergyEquation{<:Murnaghan})(v)
    v0, b0, b′0, e0 = fieldvalues(f.params)
    x, y = b′0 - 1, (v0 / v)^b′0
    return e0 + b0 / b′0 * v * (y / x + 1) - v0 * b0 / x
end
function (f::EnergyEquation{<:BirchMurnaghan2nd})(v)
    v0, b0, e0 = fieldvalues(f.params)
    f = (cbrt(v0 / v)^2 - 1) / 2
    return e0 + 9 / 2 * b0 * v0 * f^2
end
function (f::EnergyEquation{<:BirchMurnaghan3rd})(v)
    v0, b0, b′0, e0 = fieldvalues(f.params)
    x = cbrt(v0 / v)
    y = x^2 - 1
    return e0 + 9 / 16 * b0 * v0 * y^2 * (6 - 4 * x^2 + b′0 * y)
end
function (f::EnergyEquation{<:BirchMurnaghan4th})(v)
    v0, b0, b′0, b′′0, e0 = fieldvalues(f.params)
    f, h = (cbrt(v0 / v)^2 - 1) / 2, b0 * b′′0 + b′0^2
    return e0 + 3 / 8 * v0 * b0 * f^2 * ((9h - 63b′0 + 143) * f^2 + 12f * (b′0 - 4) + 12)
end
function (f::EnergyEquation{<:BirchMurnaghan5th})(v)
    v0, b0, b′0, b′′0, b′′′0, e0 = fieldvalues(f.params)
    f = (cbrt(v0 / v)^2 - 1) / 2
    c2 = 9 / 2 * b0 * v0
    c3 = c2 * (b′0 - 4)
    c4 = 3 / 8 * b0 * v0 * (9 * (b′′0 * b0 + b′0^2) - 63b′0 + 143)
    c5 =
        (
            432 * c2 * c3 * c4 + 576 * c2^2 * c4 - 243 * c3^3 - 648 * c2 * c3^2 -
            1350 * c2^2 * c3 - 2520 * c2^3
        ) / (180 * c2^2) + b′′′0 * c2^3 / (45 * v0^2)
    return e0 + f^2 * (f * (f * (f * c5 + c4) + c3) + c2)
end
function (f::EnergyEquation{<:PoirierTarantola2nd})(v)
    v0, b0, e0 = fieldvalues(f.params)
    return e0 + b0 / 2 * v0 * cbrt(log(v / v0))^2
end
function (f::EnergyEquation{<:PoirierTarantola3rd})(v)
    v0, b0, b′0, e0 = fieldvalues(f.params)
    x = log(v0 / v)
    return e0 + b0 * v0 / 6 * x^2 * ((b′0 - 2) * x + 3)
end
function (f::EnergyEquation{<:PoirierTarantola4th})(v)
    v0, b0, b′0, b′′0, e0 = fieldvalues(f.params)
    x = cbrt(v / v0)
    xi = log(x)
    h = b0 * b′′0 + b′0^2
    return e0 + b0 / 24v0 * xi^2 * ((h + 3b′0 + 3) * xi^2 + 4 * (b′0 + 2) * xi + 12)
end
function (f::EnergyEquation{<:PoirierTarantola5th})(v)
    v0, b0, b′0, b′′0, b′′′0, e0 = fieldvalues(f.params)
    f = log(v / v0) / 3
    c = 9 / 2 * b0 * v0
    d = c * (2 - b′0)
    ee = (9v0 * (d^2 - c * d + c^2) + 2 * c^3 * b′′0) / (12c * v0)
    g =
        (
            (432 * (d - c) * c * ee - 243 * d^3 + 486c * d * (d - c) + 324 * c^3) * v0^2 -
            4 * c^5 * b′′′0
        ) / (180 * (c * v0)^2)
    return e0 + f^2 * (f * (f * (f * g + ee) + d) + c)
end
function (f::EnergyEquation{<:Vinet})(v)
    v0, b0, b′0, e0 = fieldvalues(f.params)
    x, y = 1 - cbrt(v / v0), 3 / 2 * (b′0 - 1)
    return e0 + 9b0 * v0 / y^2 * (1 + (x * y - 1) * exp(x * y))
end
function (f::EnergyEquation{<:AntonSchmidt})(v)
    v0, β, n, e∞ = fieldvalues(f.params)
    x, η = v / v0, n + 1
    return e∞ + β * v0 / η * x^η * (log(x) - 1 / η)
end
function (f::EnergyEquation{<:PolynomialEOS{N}})(v) where {N}
    v0, p0, e0 = fieldvalues(f.params)
    return e0 + dot(p0, (v - v0)^n for n = 1:N)  # It cannot be `v0 - v`!
end # function _evaluate
# PressureEquation evaluation
function (f::PressureEquation{<:Murnaghan})(v)
    v0, b0, b′0 = fieldvalues(f.params)
    return b0 / b′0 * ((v0 / v)^b′0 - 1)
end
function (f::PressureEquation{<:BirchMurnaghan2nd})(v)
    v0, b0 = fieldvalues(f.params)
    f = (cbrt(v0 / v)^2 - 1) / 2
    return 3b0 * f * (1 + 2f)^(5 / 2)
end
function (f::PressureEquation{<:BirchMurnaghan3rd})(v)
    v0, b0, b′0 = fieldvalues(f.params)
    x = cbrt(v0 / v)
    return 3 / 2 * b0 * (x^7 - x^5) * (1 + 3 / 4 * (b′0 - 4) * (x^2 - 1))
end
function (f::PressureEquation{<:BirchMurnaghan4th})(v)
    v0, b0, b′0, b′′0 = fieldvalues(f.params)
    f, h = (cbrt(v0 / v)^2 - 1) / 2, b0 * b′′0 + b′0^2
    return b0 / 2 * (2f + 1)^(5 / 2) * ((9h - 63b′0 + 143) * f^2 + 9f * (b′0 - 4) + 6)
end
function (f::PressureEquation{<:PoirierTarantola2nd})(v)
    v0, b0 = fieldvalues(f.params)
    x = cbrt(v / v0)
    return -b0 / x * log(x)
end
function (f::PressureEquation{<:PoirierTarantola3rd})(v)
    v0, b0, b′0 = fieldvalues(f.params)
    x = v0 / v
    ξ = log(x)
    return b0 * x * ξ * (1 + (b′0 - 2) / 2 * ξ)
end
function (f::PressureEquation{<:PoirierTarantola4th})(v)
    v0, b0, b′0, b′′0 = fieldvalues(f.params)
    x = cbrt(v / v0)
    xi = log(x)
    h = b0 * b′′0 + b′0^2
    return -b0 * xi / 6 / x * ((h + 3b′0 + 3) * xi^2 + 3xi * (b′0 + 6) + 6)
end
function (f::PressureEquation{<:Vinet})(v)
    v0, b0, b′0 = fieldvalues(f.params)
    x, y = cbrt(v / v0), 3 // 2 * (b′0 - 1)
    return 3b0 / x^2 * (1 - x) * exp(y * (1 - x))
end
function (f::PressureEquation{<:AntonSchmidt})(v)
    v0, β, n = fieldvalues(f.params)
    x = v / v0
    return -β * x^n * log(x)
end
function (f::PressureEquation{<:BreenanStacey})(v)
    v0, b0, γ0 = fieldvalues(f.params)
    x = v0 / v
    return b0 / 2 / γ0 * x^(4 / 3) * (exp(2γ0 * (1 - x)) - 1)
end
function (f::PressureEquation{<:Shanker})(v)
    v0, b0, b′0 = fieldvalues(f.params)
    x = v / v0
    y = 1 - x
    t = b′0 - 8 / 3
    return b0 / (x^(4 / 3) * t) *
           ((1 - 1 / t + 2 / t^2) * exp(t * y - 1) + y * (1 + y - 2 / t) * exp(t * y))
end
# Bulk modulus evaluation
function (f::BulkModulusEquation{<:BirchMurnaghan2nd})(v)
    v0, b0 = fieldvalues(f.params)
    f = (cbrt(v0 / v)^2 - 1) / 2
    return b0 * (7f + 1) * (2f + 1)^(5 / 2)
end
function (f::BulkModulusEquation{<:BirchMurnaghan3rd})(v)
    v0, b0, b′0 = fieldvalues(f.params)
    f = (cbrt(v0 / v)^2 - 1) / 2
    return b0 / 2 * (2f + 1)^(5 / 2) * ((27 * f^2 + 6f) * (b′0 - 4) - 4f + 2)
end
function (f::BulkModulusEquation{<:BirchMurnaghan4th})(v)
    v0, b0, b′0, b′′0 = fieldvalues(f.params)
    f, h = (cbrt(v0 / v)^2 - 1) / 2, b0 * b′′0 + b′0^2
    return b0 / 6 *
           (2f + 1)^(5 / 2) *
           ((99h - 693b′0 + 1573) * f^3 + (27h - 108b′0 + 105) * f^2 + 6f * (3b′0 - 5) + 6)
end
function (f::BulkModulusEquation{<:PoirierTarantola2nd})(v)
    v0, b0 = fieldvalues(f.params)
    x = cbrt(v / v0)
    return b0 / x * (1 - log(x))
end
function (f::BulkModulusEquation{<:PoirierTarantola3rd})(v)
    v0, b0, b′0 = fieldvalues(f.params)
    x = v / v0
    xi = log(x)
    return -b0 / 2x * (((b′0 - 2) * xi + 2 - 2b′0) * xi + 2)
end
function (f::BulkModulusEquation{<:PoirierTarantola4th})(v)
    v0, b0, b′0, b′′0 = fieldvalues(f.params)
    x = cbrt(v / v0)
    xi = log(x)
    h = b0 * b′′0 + b′0^2
    return -b0 / (6x) *
           ((h + 3b′0 + 3) * xi^3 - 3 * xi^2 * (h + 2b′0 + 1) - 6xi * (b′0 + 1) - 6)
end
function (f::BulkModulusEquation{<:Vinet})(v)
    v0, b0, b′0 = fieldvalues(f.params)
    x, xi = cbrt(v / v0), 3 / 2 * (b′0 - 1)
    return -b0 / (2 * x^2) * (3x * (x - 1) * (b′0 - 1) + 2 * (x - 2)) * exp(-xi * (x - 1))
end
function (f::BulkModulusEquation{<:AntonSchmidt})(v)
    v0, β, n = fieldvalues(f.params)
    x = v / v0
    return β * x^n * (1 + n * log(x))
end
function (f::BulkModulusEquation{<:Shanker})(v)
    v0, b0, b′0 = fieldvalues(f.params)
    x = v / v0
    y = 1 - x
    t = b′0 - 8 / 3
    return b0 / cbrt(x) * (1 + y + y^2) * exp(t * y) + 4 / 3 * PressureEquation(f.params)(v)
end

Base.eltype(::FieldValues{<:EOSParameters{T}}) where {T} = T
# See https://docs.julialang.org/en/latest/manual/methods/#Extracting-the-type-parameter-from-a-super-type-1
Base.eltype(::Type{<:EOSParameters{T}}) where {T} = T  # Use triangular dispatch

function Base.show(io::IO, eos::EOSParameters)  # Ref: https://github.com/mauro3/Parameters.jl/blob/3c1d72b/src/Parameters.jl#L542-L549
    if get(io, :compact, false)
        Base.show_default(IOContext(io, :limit => true), eos)
    else
        # just dumping seems to give ok output, in particular for big data-sets:
        T = typeof(eos)
        println(io, T)
        for f in fieldnames(T)
            println(io, " ", f, " = ", getfield(eos, f))
        end
    end
end # function Base.show

end
