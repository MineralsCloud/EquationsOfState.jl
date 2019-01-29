"""
# module Collections



# Examples

```jldoctest
julia>
```
"""
module Collections

using GSL: sf_gamma_inc
using StaticArrays
using Unitful

export eval_energy,
    eval_pressure,
    EquationOfState,
    FiniteStrainEquationOfState,
    NonFittingParameter,
    get_parameters,
    Birch,
    Murnaghan,
    BirchMurnaghan2nd, BirchMurnaghan3rd, BirchMurnaghan4th,
    Vinet,
    PoirierTarantola2nd, PoirierTarantola3rd, PoirierTarantola4th,
    Holzapfel,
    AntonSchmidt,
    BreenanStacey

struct NonFittingParameter{T <: Real}
    data::T
end

abstract type EquationOfState{N, T <: Real} <: FieldVector{N, T} end

abstract type FiniteStrainEquationOfState{N, T} <: EquationOfState{N, T} end

struct Birch{T} <: FiniteStrainEquationOfState{3, T}
    v0::T
    b0::T
    bp0::T
end
Birch(v0, b0, bp0) = Birch(promote(v0, b0, bp0))

struct Murnaghan{T} <: EquationOfState{3, T}
    v0::T
    b0::T
    bp0::T
end
Murnaghan(v0, b0, bp0) = Murnaghan(promote(v0, b0, bp0))

struct BirchMurnaghan2nd{T} <: FiniteStrainEquationOfState{2, T}
    v0::T
    b0::T
end
BirchMurnaghan2nd(v0, b0) = BirchMurnaghan2nd(promote(v0, b0))

struct BirchMurnaghan3rd{T} <: FiniteStrainEquationOfState{3, T}
    v0::T
    b0::T
    bp0::T
end
BirchMurnaghan3rd(v0, b0, bp0) = BirchMurnaghan3rd(promote(v0, b0, bp0))

struct BirchMurnaghan4th{T} <: FiniteStrainEquationOfState{4, T}
    v0::T
    b0::T
    bp0::T
    bpp0::T
end
BirchMurnaghan4th(v0, b0, bp0, bpp0) = BirchMurnaghan4th(promote(v0, b0, bp0, bpp0))

struct Vinet{T} <: EquationOfState{3, T}
    v0::T
    b0::T
    bp0::T
end
Vinet(v0, b0, bp0) = Vinet(promote(v0, b0, bp0))

struct PoirierTarantola2nd{T} <: FiniteStrainEquationOfState{2, T}
    v0::T
    b0::T
end
PoirierTarantola2nd(v0, b0) = PoirierTarantola2nd(promote(v0, b0))

struct PoirierTarantola3rd{T} <: FiniteStrainEquationOfState{3, T}
    v0::T
    b0::T
    bp0::T
end
PoirierTarantola3rd(v0, b0, bp0) = PoirierTarantola3rd(promote(v0, b0, bp0))

struct PoirierTarantola4th{T} <: FiniteStrainEquationOfState{4, T}
    v0::T
    b0::T
    bp0::T
    bpp0::T
end
PoirierTarantola4th(v0, b0, bp0, bpp0) = PoirierTarantola4th(promote(v0, b0, bp0, bpp0))

struct Holzapfel{T} <: EquationOfState{4, T}
    v0::T
    b0::T
    bp0::T
    z::NonFittingParameter
end
Holzapfel(v0, b0, bp0, z::NonFittingParameter) = Holzapfel(promote(v0, b0, bp0)..., z)
Holzapfel(v0, b0, bp0, z::Number) = Holzapfel(v0, b0, bp0, NonFittingParameter(z))

struct AntonSchmidt{T} <: EquationOfState{3, T}
    v0::T
    β::T
    n::T
end
AntonSchmidt(v0, β, n) = AntonSchmidt(promote(v0, β, n))

struct BreenanStacey{T} <: EquationOfState{3, T}
    v0::T
    b0::T
    γ0::T
end
BreenanStacey(v0, b0, γ0) = BreenanStacey(promote(v0, b0, γ0))

function get_parameters(eos::T) where {T <: EquationOfState}
    map(f -> getfield(eos, f), fieldnames(T))
end

function eval_energy(eos::Birch)::Function
    v0, b0, bp0 = get_parameters(eos)

    function (v::T, e0 = zero(T)) where {T <: Real}
        x = (v0 / v)^(2 / 3) - 1
        xi = 9 / 16 * b0 * v0 * x^2
        return e0 + 2 * xi + (bp0 - 4) * xi * x
    end
end

function eval_pressure(eos::Birch)::Function
    v0, b0, bp0 = get_parameters(eos)

    function (v::Real)
        x = v0 / v
        xi = x^(2 / 3) - 1
        return 3 / 8 * b0 * x^(5 / 3) * xi * (4 + 3 * (bp0 - 4) * xi)
    end
end

function eval_energy(eos::Murnaghan)::Function
    v0, b0, bp0 = get_parameters(eos)

    function (v::T, e0 = zero(T)) where {T <: Real}
        x = bp0 - 1
        y = (v0 / v)^bp0
        return e0 + b0 / bp0 * v * (y / x + 1) - v0 * b0 / x
    end
end

function eval_pressure(eos::Murnaghan)::Function
    v0, b0, bp0 = get_parameters(eos)

    function (v::Real)
        return b0 / bp0 * ((v0 / v)^bp0 - 1)
    end
end

function eval_energy(eos::BirchMurnaghan2nd)::Function
    v0, b0 = get_parameters(eos)

    function (v::T, e0 = zero(T)) where {T <: Real}
        f = ((v0 / v)^(2 / 3) - 1) / 2
        return e0 + 9 / 2 * b0 * v0 * f^2
    end
end

function eval_pressure(eos::BirchMurnaghan2nd)::Function
    v0, b0 = get_parameters(eos)

    function (v::Real)
        f = ((v0 / v)^(2 / 3) - 1) / 2
        return 3 * b0 * f * (1 + 2 * f)^(5 / 2)
    end
end

function eval_energy(eos::BirchMurnaghan3rd)::Function
    v0, b0, bp0 = get_parameters(eos)

    function (v::T, e0 = zero(T)) where {T <: Real}
        eta = (v0 / v)^(1 / 3)
        xi = eta^2 - 1
        return e0 + 9 / 16 * b0 * v0 * xi^2 * (6 + bp0 * xi - 4 * eta^2)
    end
end

function eval_pressure(eos::BirchMurnaghan3rd)::Function
    v0, b0, bp0 = get_parameters(eos)

    function (v::Real)
        eta = (v0 / v)^(1 / 3)
        return 3 / 2 * b0 * (eta^7 - eta^5) * (1 + 3 / 4 * (bp0 - 4) * (eta^2 - 1))
    end
end

function eval_energy(eos::BirchMurnaghan4th)::Function
    v0, b0, bp0, bpp0 = get_parameters(eos)

    function (v::T, e0 = zero(T)) where {T <: Real}
        f = ((v0 / v)^(2 / 3) - 1) / 2
        h = b0 * bpp0 + bp0^2
        return e0 + 3 / 8 * v0 * b0 * f^2 * ((9 * h - 63 * bp0 + 143) * f^2 + 12 * (bp0 - 4) * f + 12)
    end
end

function eval_pressure(eos::BirchMurnaghan4th)::Function
    v0, b0, bp0, bpp0 = get_parameters(eos)

    function (v::Real)
        f = ((v0 / v)^(2 / 3) - 1) / 2
        h = b0 * bpp0 + bp0^2
        return 1 / 2 * b0 * (2 * f + 1)^(5 / 2) * ((9 * h - 63 * bp0 + 143) * f^2 + 9 * (bp0 - 4) * f + 6)
    end
end

function eval_energy(eos::Vinet)::Function
    v0, b0, bp0 = get_parameters(eos)

    function (v::T, e0 = zero(T)) where {T <: Real}
        x = (v / v0)^(1 / 3)
        xi = 3 / 2 * (bp0 - 1)
        return e0 + 9 * b0 * v0 / xi^2 * (1 + (xi * (1 - x) - 1) * exp(xi * (1 - x)))
    end
end

function eval_pressure(eos::Vinet)::Function
    v0, b0, bp0 = get_parameters(eos)

    function (v::Real)
        x = (v / v0)^(1 / 3)
        xi = 3 / 2 * (bp0 - 1)
        return 3 * b0 / x^2 * (1 - x) * exp(xi * (1 - x))
    end
end

function eval_energy(eos::PoirierTarantola2nd)::Function
    v0, b0 = get_parameters(eos)

    function (v::T, e0 = zero(T)) where {T <: Real}
        return e0 + 1 / 2 * b0 * v0 * log(v / v0)^(2 / 3)
    end
end

function eval_pressure(eos::PoirierTarantola2nd)::Function
    v0, b0 = get_parameters(eos)

    function (v::Real)
        x = (v / v0)^(1 / 3)
        return -b0 / x * log(x)
    end
end

function eval_energy(eos::PoirierTarantola3rd)::Function
    v0, b0, bp0 = get_parameters(eos)

    function (v::T, e0 = zero(T)) where {T <: Real}
        x = (v / v0)^(1 / 3)
        xi = log(x)
        return e0 + 1 / 6 * b0 * v0 * xi^2 * ((bp0 + 2) * xi + 3)
    end
end

function eval_pressure(eos::PoirierTarantola3rd)::Function
    v0, b0, bp0 = get_parameters(eos)

    function (v::Real)
        x = (v / v0)^(1 / 3)
        xi = log(x)
        return -b0 * xi / (2 * x) * ((bp0 + 2) * xi + 2)
    end
end

function eval_energy(eos::PoirierTarantola4th)::Function
    v0, b0, bp0, bpp0 = get_parameters(eos)

    function (v::T, e0 = zero(T)) where {T <: Real}
        x = (v / v0)^(1 / 3)
        xi = log(x)
        h = b0 * bpp0 + bp0^2
        return e0 + 1 / 24 * b0 * v0 * xi^2 * ((h + 3 * bp0 + 3) * xi^2 + 4 * (bp0 + 2) * xi + 12)
    end
end

function eval_pressure(eos::PoirierTarantola4th)::Function
    v0, b0, bp0, bpp0 = get_parameters(eos)

    function (v::Real)
        x = (v / v0)^(1 / 3)
        xi = log(x)
        h = b0 * bpp0 + bp0^2
        return -b0 * xi / 6 / x * ((h + 3 * bp0 + 3) * xi^2 + 3 * (bp0 + 6) * xi + 6)
    end
end

function eval_energy(eos::Holzapfel)::Function
    v0, b0, bp0, z = get_parameters(eos)

    function (v::T, e0 = zero(T)) where {T <: Real}
        η = (v / v0)^(1 / 3)
        pfg0 = 3.8283120002509214 * (z / v0)^(5 / 3)
        c0 = -log(3 * b0 / pfg0)
        c2 = 3 / 2 * (bp0 - 3) - c0
        term1 = (sf_gamma_inc(-2, c0 * η) - sf_gamma_inc(-2, c0)) * c0^2 * exp(c0)
        term2 = (sf_gamma_inc(-1, c0 * η) - sf_gamma_inc(-1, c0)) * c0 * (c2 - 1) * exp(c0)
        term3 = (sf_gamma_inc(0, c0 * η) - sf_gamma_inc(0, c0)) * 2 * c2 * exp(c0)
        term4 = c2 / c0 * (exp(c0 * (1 - η)) - 1)
        return e0 + 9 * b0 * v0 * (term1 + term2 - term3 + term4)
    end
end

function eval_pressure(eos::Holzapfel)::Function
    v0, b0, bp0, z = get_parameters(eos)

    function (v::Real)
        η = (v / v0)^(1 / 3)
        pfg0 = 3.8283120002509214 * (z / v0)^(5 / 3)
        c0 = -log(3 * b0 / pfg0)
        c2 = 3 / 2 * (bp0 - 3) - c0
        return p0 + 3 * b0 * (1 - η) / η^5 * exp(c0 * (1 - η)) * (1 + c2 * η * (1 - η))
    end
end

function eval_energy(eos::AntonSchmidt)::Function
    v0, β, n = get_parameters(eos)

    function (v::T, e∞::T=0) where {T <: Real}
        x = v / v0
        η = n + 1
        return e∞ + β * v0 / η * x^η * (log(x) - 1 / η)
    end
end

function eval_pressure(eos::AntonSchmidt)::Function
    v0, β, n = get_parameters(eos)

    function (v::Real)
        x = v / v0
        return -β * x^n * log(x)
    end
end

function eval_pressure(eos::BreenanStacey)::Function
    v0, b0, γ0 = get_parameters(eos)

    function (v::Real)
        x = v0 / v
        return b0 / 2 / γ0 * x^(4 / 3) * (exp(2 * γ0 * (1 - x)) - 1)
    end
end

end