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
using UnitfulHartree

export eval_energy,
    eval_pressure,
    EquationOfState,
    FiniteStrainEquationOfState,
    NonFittingParameter,
    collect_parameters,
    Birch,
    Murnaghan,
    BirchMurnaghan2nd, BirchMurnaghan3rd, BirchMurnaghan4th,
    Vinet,
    PoirierTarantola2nd, PoirierTarantola3rd, PoirierTarantola4th,
    Holzapfel,
    AntonSchmidt

const VolumeType = typeof(1.0u"a₀^3")
const BulkModulusType = typeof(1.0u"GPa")
const BulkModulus2ndDerivativeType = typeof(1.0u"GPa^-1")
const EnergyType = typeof(1.0u"eV")

abstract type EquationOfState{N, Float64} <: FieldVector{N, Float64} end

abstract type FiniteStrainEquationOfState{N, Float64} <: EquationOfState{N, Float64} end

primitive type NonFittingParameter <: AbstractFloat 64 end

struct Birch <: FiniteStrainEquationOfState{3, Float64}
    v0::VolumeType
    b0::BulkModulusType
    bp0::Unitful.DimensionlessQuantity{Float64}
end

struct Murnaghan <: EquationOfState{3, Float64}
    v0::VolumeType
    b0::BulkModulusType
    bp0::Unitful.DimensionlessQuantity{Float64}
end

struct BirchMurnaghan2nd <: FiniteStrainEquationOfState{2, Float64}
    v0::VolumeType
    b0::BulkModulusType
end

struct BirchMurnaghan3rd <: FiniteStrainEquationOfState{3, Float64}
    v0::VolumeType
    b0::BulkModulusType
    bp0::Unitful.DimensionlessQuantity{Float64}
end

struct BirchMurnaghan4th <: FiniteStrainEquationOfState{4, Float64}
    v0::VolumeType
    b0::BulkModulusType
    bp0::Unitful.DimensionlessQuantity{Float64}
    bpp0::BulkModulus2ndDerivativeType
end

struct Vinet <: EquationOfState{3, Float64}
    v0::VolumeType
    b0::BulkModulusType
    bp0::Unitful.DimensionlessQuantity{Float64}
end

struct PoirierTarantola2nd <: FiniteStrainEquationOfState{2, Float64}
    v0::VolumeType
    b0::BulkModulusType
end

struct PoirierTarantola3rd <: FiniteStrainEquationOfState{3, Float64}
    v0::VolumeType
    b0::BulkModulusType
    bp0::Unitful.DimensionlessQuantity{Float64}
end

struct PoirierTarantola4th <: FiniteStrainEquationOfState{4, Float64}
    v0::VolumeType
    b0::BulkModulusType
    bp0::Unitful.DimensionlessQuantity{Float64}
    bpp0::BulkModulus2ndDerivativeType
end

struct Holzapfel <: EquationOfState{4, Float64}
    v0::VolumeType
    b0::BulkModulusType
    bp0::Unitful.DimensionlessQuantity{Float64}
    z::NonFittingParameter
end

struct AntonSchmidt <: EquationOfState{3, Float64}
    v0::VolumeType
    β::Float64
    n::Float64
end

struct BreenanStacey <: EquationOfState{3, Float64}
    v0::VolumeType
    b0::BulkModulusType
    γ0::Float64
end

function collect_parameters(eos::T) where {T <: EquationOfState}
    map(f -> getfield(eos, f), fieldnames(T)) |> collect
end

function eval_energy(eos::Birch)::Function
    v0, b0, bp0 = collect_parameters(eos)

    function (v::VolumeType, e0::EnergyType=0.0u"eV")
        x = (v0 / v)^(2 / 3) - 1
        xi = 9 / 16 * b0 * v0 * x^2
        return e0 + 2 * xi + (bp0 - 4) * xi * x
    end
end

function eval_pressure(eos::Birch)::Function
    v0, b0, bp0 = collect_parameters(eos)

    function (v::VolumeType)
        x = v0 / v
        xi = x^(2 / 3) - 1
        return 3 / 8 * b0 * x^(5 / 3) * xi * (4 + 3 * (bp0 - 4) * xi)
    end
end

function eval_energy(eos::Murnaghan)::Function
    v0, b0, bp0 = collect_parameters(eos)

    function (v::VolumeType, e0::EnergyType=0.0u"eV")
        x = bp0 - 1
        y = (v0 / v)^bp0
        return e0 + b0 / bp0 * v * (y / x + 1) - v0 * b0 / x
    end
end

function eval_pressure(eos::Murnaghan)::Function
    v0, b0, bp0 = collect_parameters(eos)

    function (v::VolumeType)
        return b0 / bp0 * ((v0 / v)^bp0 - 1)
    end
end

function eval_energy(eos::BirchMurnaghan2nd)::Function
    v0, b0 = collect_parameters(eos)

    function (v::VolumeType, e0::EnergyType=0.0u"eV")
        f = ((v0 / v)^(2 / 3) - 1) / 2
        return e0 + 9 / 2 * b0 * v0 * f^2
    end
end

function eval_pressure(eos::BirchMurnaghan2nd)::Function
    v0, b0 = collect_parameters(eos)

    function (v::VolumeType)
        f = ((v0 / v)^(2 / 3) - 1) / 2
        return 3 * b0 * f * (1 + 2 * f)^(5 / 2)
    end
end

function eval_energy(eos::BirchMurnaghan3rd)::Function
    v0, b0, bp0 = collect_parameters(eos)

    function (v::VolumeType, e0::EnergyType=0.0u"eV")
        eta = (v0 / v)^(1 / 3)
        xi = eta^2 - 1
        return e0 + 9 / 16 * b0 * v0 * xi^2 * (6 + bp0 * xi - 4 * eta^2)
    end
end

function eval_pressure(eos::BirchMurnaghan3rd)::Function
    v0, b0, bp0 = collect_parameters(eos)

    function (v::VolumeType)
        eta = (v0 / v)^(1 / 3)
        return 3 / 2 * b0 * (eta^7 - eta^5) * (1 + 3 / 4 * (bp0 - 4) * (eta^2 - 1))
    end
end

function eval_energy(eos::BirchMurnaghan4th)::Function
    v0, b0, bp0, bpp0 = collect_parameters(eos)

    function (v::VolumeType, e0::EnergyType=0.0u"eV")
        f = ((v0 / v)^(2 / 3) - 1) / 2
        h = b0 * bpp0 + bp0^2
        return e0 + 3 / 8 * v0 * b0 * f^2 * ((9 * h - 63 * bp0 + 143) * f^2 + 12 * (bp0 - 4) * f + 12)
    end
end

function eval_pressure(eos::BirchMurnaghan4th)::Function
    v0, b0, bp0, bpp0 = collect_parameters(eos)

    function (v::VolumeType)
        f = ((v0 / v)^(2 / 3) - 1) / 2
        h = b0 * bpp0 + bp0^2
        return 1 / 2 * b0 * (2 * f + 1)^(5 / 2) * ((9 * h - 63 * bp0 + 143) * f^2 + 9 * (bp0 - 4) * f + 6)
    end
end

function eval_energy(eos::Vinet)::Function
    v0, b0, bp0 = collect_parameters(eos)

    function (v::VolumeType, e0::EnergyType=0.0u"eV")
        x = (v / v0)^(1 / 3)
        xi = 3 / 2 * (bp0 - 1)
        return e0 + 9 * b0 * v0 / xi^2 * (1 + (xi * (1 - x) - 1) * exp(xi * (1 - x)))
    end
end

function eval_pressure(eos::Vinet)::Function
    v0, b0, bp0 = collect_parameters(eos)

    function (v::VolumeType)
        x = (v / v0)^(1 / 3)
        xi = 3 / 2 * (bp0 - 1)
        return 3 * b0 / x^2 * (1 - x) * exp(xi * (1 - x))
    end
end

function eval_energy(eos::PoirierTarantola2nd)::Function
    v0, b0 = collect_parameters(eos)

    function (v::VolumeType, e0::EnergyType=0.0u"eV")
        return e0 + 1 / 2 * b0 * v0 * log(v / v0)^(2 / 3)
    end
end

function eval_pressure(eos::PoirierTarantola2nd)::Function
    v0, b0 = collect_parameters(eos)

    function (v::VolumeType)
        x = (v / v0)^(1 / 3)
        return -b0 / x * log(x)
    end
end

function eval_energy(eos::PoirierTarantola3rd)::Function
    v0, b0, bp0 = collect_parameters(eos)

    function (v::VolumeType, e0::EnergyType=0.0u"eV")
        x = (v / v0)^(1 / 3)
        xi = log(x)
        return e0 + 1 / 6 * b0 * v0 * xi^2 * ((bp0 + 2) * xi + 3)
    end
end

function eval_pressure(eos::PoirierTarantola3rd)::Function
    v0, b0, bp0 = collect_parameters(eos)

    function (v::VolumeType)
        x = (v / v0)^(1 / 3)
        xi = log(x)
        return -b0 * xi / (2 * x) * ((bp0 + 2) * xi + 2)
    end
end

function eval_energy(eos::PoirierTarantola4th)::Function
    v0, b0, bp0, bpp0 = collect_parameters(eos)

    function (v::VolumeType, e0::EnergyType=0.0u"eV")
        x = (v / v0)^(1 / 3)
        xi = log(x)
        h = b0 * bpp0 + bp0^2
        return e0 + 1 / 24 * b0 * v0 * xi^2 * ((h + 3 * bp0 + 3) * xi^2 + 4 * (bp0 + 2) * xi + 12)
    end
end

function eval_pressure(eos::PoirierTarantola4th)::Function
    v0, b0, bp0, bpp0 = collect_parameters(eos)

    function (v::VolumeType)
        x = (v / v0)^(1 / 3)
        xi = log(x)
        h = b0 * bpp0 + bp0^2
        return -b0 * xi / 6 / x * ((h + 3 * bp0 + 3) * xi^2 + 3 * (bp0 + 6) * xi + 6)
    end
end

function eval_energy(eos::Holzapfel)::Function
    v0, b0, bp0, z = collect_parameters(eos)

    function (v::VolumeType, e0::EnergyType=0.0u"eV")
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
    v0, b0, bp0, z = collect_parameters(eos)

    function (v::VolumeType)
        η = (v / v0)^(1 / 3)
        pfg0 = 3.8283120002509214 * (z / v0)^(5 / 3)
        c0 = -log(3 * b0 / pfg0)
        c2 = 3 / 2 * (bp0 - 3) - c0
        return p0 + 3 * b0 * (1 - η) / η^5 * exp(c0 * (1 - η)) * (1 + c2 * η * (1 - η))
    end
end

function eval_energy(eos::AntonSchmidt)::Function
    v0, β, n = collect_parameters(eos)

    function (v::VolumeType, e∞::EnergyType=0)
        x = v / v0
        η = n + 1
        return e∞ + β * v0 / η * x^η * (log(x) - 1 / η)
    end
end

function eval_pressure(eos::AntonSchmidt)::Function
    v0, β, n = collect_parameters(eos)

    function (v::VolumeType)
        x = v / v0
        return -β * x^n * log(x)
    end
end

function eval_pressure(eos::BreenanStacey)::Function
    v0, b0, γ0 = collect_parameters(eos)

    function (v::VolumeType)
        x = v0 / v
        return b0 / 2 / γ0 * x^(4 / 3) * (exp(2 * γ0 * (1 - x)) - 1)
    end
end

end