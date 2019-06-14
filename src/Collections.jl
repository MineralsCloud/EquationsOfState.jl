"""
# module Collections



# Examples

```jldoctest
julia>
```
"""
module Collections

using InteractiveUtils
using Parameters
using StaticArrays: FieldVector, Size

import StaticArrays: similar_type

export eval_energy,
    eval_pressure,
    eval_bulk_modulus,
    EquationOfState,
    FiniteStrainEquationOfState,
    Birch,
    Murnaghan,
    BirchMurnaghan2nd, BirchMurnaghan3rd, BirchMurnaghan4th,
    PoirierTarantola2nd, PoirierTarantola3rd, PoirierTarantola4th,
    Vinet,
    AntonSchmidt,
    BreenanStacey,
    similar_type

# ============================================================================ #
#                                     Types                                    #
# ============================================================================ #
abstract type EquationOfState{T,N} <: FieldVector{N,T} end

abstract type FiniteStrainEquationOfState{T,N} <: EquationOfState{T,N} end

@with_kw struct Birch{T <: Real} <: FiniteStrainEquationOfState{T,4}
    v0::T
    b0::T
    bp0::T
    e0::T = 0
end
function Birch(v0::Real, b0::Real, bp0::Real, e0::Real)
    T = Base.promote_typeof(v0, b0, bp0, e0)
    Birch{T}(v0, b0, bp0, e0)
end
Birch(v0, b0, bp0) = Birch(v0, b0, bp0, 0)

@with_kw struct Murnaghan{T <: Real} <: EquationOfState{T,4}
    v0::T
    b0::T
    bp0::T
    e0::T = 0
end
function Murnaghan(v0::Real, b0::Real, bp0::Real, e0::Real)
    T = Base.promote_typeof(v0, b0, bp0, e0)
    Murnaghan{T}(v0, b0, bp0, e0)
end
Murnaghan(v0, b0, bp0) = Murnaghan(v0, b0, bp0, 0)

@with_kw struct BirchMurnaghan2nd{T <: Real} <: FiniteStrainEquationOfState{T,3}
    v0::T
    b0::T
    e0::T = 0
end
function BirchMurnaghan2nd(v0::Real, b0::Real, e0::Real)
    T = Base.promote_typeof(v0, b0, e0)
    BirchMurnaghan2nd{T}(v0, b0, e0)
end
BirchMurnaghan2nd(v0, b0) = BirchMurnaghan2nd(v0, b0, 0)

@with_kw struct BirchMurnaghan3rd{T <: Real} <: FiniteStrainEquationOfState{T,4}
    v0::T
    b0::T
    bp0::T
    e0::T = 0
end
function BirchMurnaghan3rd(v0::Real, b0::Real, bp0::Real, e0::Real)
    T = Base.promote_typeof(v0, b0, bp0, e0)
    BirchMurnaghan3rd{T}(v0, b0, bp0, e0)
end
BirchMurnaghan3rd(v0, b0, bp0) = BirchMurnaghan3rd(v0, b0, bp0, 0)

@with_kw struct BirchMurnaghan4th{T <: Real} <: FiniteStrainEquationOfState{T,5}
    v0::T
    b0::T
    bp0::T
    bpp0::T
    e0::T = 0
end
function BirchMurnaghan4th(v0::Real, b0::Real, bp0::Real, bpp0::Real, e0::Real)
    T = Base.promote_typeof(v0, b0, bp0, bpp0, e0)
    BirchMurnaghan4th{T}(v0, b0, bp0, bpp0, e0)
end
BirchMurnaghan4th(v0, b0, bp0, bpp0) = BirchMurnaghan4th(v0, b0, bp0, bpp0, 0)

@with_kw struct PoirierTarantola2nd{T <: Real} <: FiniteStrainEquationOfState{T,3}
    v0::T
    b0::T
    e0::T = 0
end
function PoirierTarantola2nd(v0::Real, b0::Real, e0::Real)
    T = Base.promote_typeof(v0, b0, e0)
    PoirierTarantola2nd{T}(v0, b0, e0)
end
PoirierTarantola2nd(v0, b0) = PoirierTarantola2nd(v0, b0, 0)

@with_kw struct PoirierTarantola3rd{T <: Real} <: FiniteStrainEquationOfState{T,4}
    v0::T
    b0::T
    bp0::T
    e0::T = 0
end
function PoirierTarantola3rd(v0::Real, b0::Real, bp0::Real, e0::Real)
    T = Base.promote_typeof(v0, b0, bp0, e0)
    PoirierTarantola3rd{T}(v0, b0, bp0, e0)
end
PoirierTarantola3rd(v0, b0, bp0) = PoirierTarantola3rd(v0, b0, bp0, 0)

@with_kw struct PoirierTarantola4th{T <: Real} <: FiniteStrainEquationOfState{T,5}
    v0::T
    b0::T
    bp0::T
    bpp0::T
    e0::T = 0
end
function PoirierTarantola4th(v0::Real, b0::Real, bp0::Real, bpp0::Real, e0::Real)
    T = Base.promote_typeof(v0, b0, bp0, bpp0, e0)
    PoirierTarantola4th{T}(v0, b0, bp0, bpp0, e0)
end
PoirierTarantola4th(v0, b0, bp0, bpp0) = PoirierTarantola4th(v0, b0, bp0, bpp0, 0)

@with_kw struct Vinet{T <: Real} <: EquationOfState{T,4}
    v0::T
    b0::T
    bp0::T
    e0::T = 0
end
function Vinet(v0::Real, b0::Real, bp0::Real, e0::Real)
    T = Base.promote_typeof(v0, b0, bp0, e0)
    Vinet{T}(v0, b0, bp0, e0)
end
Vinet(v0, b0, bp0) = Vinet(v0, b0, bp0, 0)

@with_kw struct AntonSchmidt{T <: Real} <: EquationOfState{T,4}
    v0::T
    β::T
    n::T
    e∞::T = 0
end
function AntonSchmidt(v0::Real, β::Real, n::Real, e∞::Real)
    T = Base.promote_typeof(v0, β, n, e∞)
    AntonSchmidt{T}(v0, β, n, e∞)
end
AntonSchmidt(v0, β, n) = AntonSchmidt(v0, β, n, 0)

@with_kw struct BreenanStacey{T <: Real} <: EquationOfState{T,4}
    v0::T
    b0::T
    γ0::T
    e0::T = 0
end
function BreenanStacey(v0::Real, b0::Real, γ0::Real, e0::Real)
    T = Base.promote_typeof(v0, b0, γ0, e0)
    BreenanStacey{T}(v0, b0, γ0, e0)
end
BreenanStacey(v0, b0, γ0) = BreenanStacey(v0, b0, γ0, 0)
# =================================== Types ================================== #


# ============================================================================ #
#                               Energy evaluation                              #
# ============================================================================ #
eval_energy(eos::EquationOfState, vec::AbstractVector{<: Real}) = map(x->eval_energy(eos, x), vec)
function eval_energy(eos::Birch, v::Real)
    @unpack v0, b0, bp0, e0 = eos

    x = (v0 / v)^(2 / 3) - 1
    xi = 9 / 16 * b0 * v0 * x^2
    return e0 + 2 * xi + (bp0 - 4) * xi * x
end
function eval_energy(eos::Murnaghan, v::Real)
    @unpack v0, b0, bp0, e0 = eos

    x = bp0 - 1
    y = (v0 / v)^bp0
    return e0 + b0 / bp0 * v * (y / x + 1) - v0 * b0 / x
end
function eval_energy(eos::BirchMurnaghan2nd, v::Real)
    @unpack v0, b0, e0 = eos

    f = ((v0 / v)^(2 / 3) - 1) / 2
    return e0 + 9 / 2 * b0 * v0 * f^2
end
function eval_energy(eos::BirchMurnaghan3rd, v::Real)
    @unpack v0, b0, bp0, e0 = eos

    eta = (v0 / v)^(1 / 3)
    xi = eta^2 - 1
    return e0 + 9 / 16 * b0 * v0 * xi^2 * (6 + bp0 * xi - 4eta^2)
end
function eval_energy(eos::BirchMurnaghan4th, v::Real)
    @unpack v0, b0, bp0, bpp0, e0 = eos

    f = ((v0 / v)^(2 / 3) - 1) / 2
    h = b0 * bpp0 + bp0^2
    return e0 + 3 / 8 * v0 * b0 * f^2 * ((9h - 63bp0 + 143) * f^2 + 12(bp0 - 4) * f + 12)
end
function eval_energy(eos::PoirierTarantola2nd, v::Real)
    @unpack v0, b0, e0 = eos

    return e0 + b0 / 2 * v0 * log(v / v0)^(2 / 3)
end
function eval_energy(eos::PoirierTarantola3rd, v::Real)
    @unpack v0, b0, bp0, e0 = eos

    x = (v / v0)^(1 / 3)
    xi = -3log(x)
    return e0 + b0 / 6 * v0 * xi^2 * ((bp0 - 2) * xi + 3)
end
function eval_energy(eos::PoirierTarantola4th, v::Real)
    @unpack v0, b0, bp0, bpp0, e0 = eos

    x = (v / v0)^(1 / 3)
    xi = log(x)
    h = b0 * bpp0 + bp0^2
    return e0 + b0 / 24v0 * xi^2 * ((h + 3bp0 + 3) * xi^2 + 4(bp0 + 2) * xi + 12)
end
function eval_energy(eos::Vinet, v::Real)
    @unpack v0, b0, bp0, e0 = eos

    x = (v / v0)^(1 / 3)
    xi = 3 / 2 * (bp0 - 1)
    return e0 + 9b0 * v0 / xi^2 * (1 + (xi * (1 - x) - 1) * exp(xi * (1 - x)))
end
function eval_energy(eos::AntonSchmidt, v::Real)
    @unpack v0, β, n, e∞ = eos

    x = v / v0
    η = n + 1
    return e∞ + β * v0 / η * x^η * (log(x) - 1 / η)
end
# ============================= Energy evaluation ============================ #


# ============================================================================ #
#                              Pressure evaluation                             #
# ============================================================================ #
eval_pressure(eos::EquationOfState, vec::AbstractVector{<: Real}) = map(x->eval_pressure(eos, x), vec)
function eval_pressure(eos::Birch, v::Real)
    @unpack v0, b0, bp0 = eos

    x = v0 / v
    xi = x^(2 / 3) - 1
    return 3 / 8 * b0 * x^(5 / 3) * xi * (4 + 3(bp0 - 4) * xi)
end
function eval_pressure(eos::Murnaghan, v::Real)
    @unpack v0, b0, bp0 = eos

    return b0 / bp0 * ((v0 / v)^bp0 - 1)
end
function eval_pressure(eos::BirchMurnaghan2nd, v::Real)
    @unpack v0, b0 = eos

    f = ((v0 / v)^(2 / 3) - 1) / 2
    return 3b0 * f * (1 + 2f)^(5 / 2)
end
function eval_pressure(eos::BirchMurnaghan3rd, v::Real)
    @unpack v0, b0, bp0 = eos

    eta = (v0 / v)^(1 / 3)
    return 3 / 2 * b0 * (eta^7 - eta^5) * (1 + 3 / 4 * (bp0 - 4) * (eta^2 - 1))
end
function eval_pressure(eos::BirchMurnaghan4th, v::Real)
    @unpack v0, b0, bp0, bpp0 = eos

    f = ((v0 / v)^(2 / 3) - 1) / 2
    h = b0 * bpp0 + bp0^2
    return b0 / 2 * (2f + 1)^(5 / 2) * ((9h - 63bp0 + 143) * f^2 + 9(bp0 - 4) * f + 6)
end
function eval_pressure(eos::PoirierTarantola2nd, v::Real)
    @unpack v0, b0 = eos

    x = (v / v0)^(1 / 3)
    return -b0 / x * log(x)
end
function eval_pressure(eos::PoirierTarantola3rd, v::Real)
    @unpack v0, b0, bp0 = eos

    x = (v / v0)^(1 / 3)
    xi = log(x)
    return -b0 * xi / (2x) * ((bp0 + 2) * xi + 2)
end
function eval_pressure(eos::PoirierTarantola4th, v::Real)
    @unpack v0, b0, bp0, bpp0 = eos

    x = (v / v0)^(1 / 3)
    xi = log(x)
    h = b0 * bpp0 + bp0^2
    return -b0 * xi / 6 / x * ((h + 3bp0 + 3) * xi^2 + 3(bp0 + 6) * xi + 6)
end
function eval_pressure(eos::Vinet, v::Real)
    @unpack v0, b0, bp0 = eos

    x = (v / v0)^(1 / 3)
    xi = 3 / 2 * (bp0 - 1)
    return 3b0 / x^2 * (1 - x) * exp(xi * (1 - x))
end
function eval_pressure(eos::AntonSchmidt, v::Real)
    @unpack v0, β, n = eos

    x = v / v0
    return -β * x^n * log(x)
end
function eval_pressure(eos::BreenanStacey, v::Real)
    @unpack v0, b0, γ0 = eos

    x = v0 / v
    return b0 / 2 / γ0 * x^(4 / 3) * (exp(2γ0 * (1 - x)) - 1)
end
# ============================ Pressure evaluation =========================== #


# ============================================================================ #
#                            Bulk modulus evaluation                           #
# ============================================================================ #
eval_bulk_modulus(eos::EquationOfState, vec::AbstractVector{<: Real}) = map(x->eval_bulk_modulus(eos, x), vec)
function eval_bulk_modulus(eos::BirchMurnaghan2nd, v::Real)
    @unpack v0, b0 = eos

    f = ((v0 / v)^(2 / 3) - 1) / 2
    return b0 * (7f + 1) * (2f + 1)^(5 / 2)
end
function eval_bulk_modulus(eos::BirchMurnaghan3rd, v::Real)
    @unpack v0, b0, bp0 = eos

    f = ((v0 / v)^(2 / 3) - 1) / 2
    return b0 / 2 * (2f + 1)^(5 / 2) * ((27f^2 + 6f) * (bp0 - 4) - 4f + 2)
end
function eval_bulk_modulus(eos::BirchMurnaghan4th, v::Real)
    @unpack v0, b0, bp0, bpp0 = eos

    f = ((v0 / v)^(2 / 3) - 1) / 2
    h = b0 * bpp0 + bp0^2
    return b0 / 6 * (2f + 1)^(5 / 2) * ((99h - 693bp0 + 1573) * f^3 + (27h - 108bp0 + 105) * f^2 + 6f * (3bp0 - 5) + 6)
end
function eval_bulk_modulus(eos::PoirierTarantola2nd, v::Real)
    @unpack v0, b0 = eos

    x = (v / v0)^(1 / 3)
    return b0 / x * (1 - log(x))
end
function eval_bulk_modulus(eos::PoirierTarantola3rd, v::Real)
    @unpack v0, b0, bp0 = eos

    x = (v / v0)^(1 / 3)
    xi = log(x)
    return -b0 / (2x) * ((bp0 + 2) * xi * (xi - 1) - 2)
end
function eval_bulk_modulus(eos::PoirierTarantola4th, v::Real)
    @unpack v0, b0, bp0, bpp0 = eos

    x = (v / v0)^(1 / 3)
    xi = log(x)
    h = b0 * bpp0 + bp0^2
    return -b0 / (6x) * ((h + 3bp0 + 3) * xi^3 - 3xi^2 * (h + 2bp0 + 1) - 6xi * (bp0 + 1) - 6)
end
function eval_bulk_modulus(eos::Vinet, v::Real)
    @unpack v0, b0, bp0 = eos

    x = (v / v0)^(1 / 3)
    xi = 3 / 2 * (bp0 - 1)
    return -b0 / (2x^2) * (3x * (x - 1) * (bp0 - 1) + 2(x - 2)) * exp(-xi * (x - 1))
end
function eval_bulk_modulus(eos::AntonSchmidt, v::Real)
    @unpack v0, β, n = eos

    x = v / v0
    return β * x^n * (1 + n * log(x))
end
# ========================== Bulk modulus evaluation ========================= #


# ============================================================================ #
#                                 Miscellaneous                                #
# ============================================================================ #
function allsubtypes(t::Type, types = Type[])::Vector{Type}
    for s in subtypes(t)
        types = allsubtypes(s, push!(types, s))
    end
    types
end

nonabstract(t::Type)::Vector{Type} = filter(!isabstracttype, allsubtypes(t))

for E in nonabstract(EquationOfState)
    eval(quote
        similar_type(::Type{A}, ::Type{T}, size::Size{(fieldcount($E),)}) where {A <: $E,T} = $E{T}
    end)
end
# =============================== Miscellaneous ============================== #

end
