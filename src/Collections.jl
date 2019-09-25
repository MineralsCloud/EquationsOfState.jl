"""
# module Collections



# Examples

```jldoctest
julia>
```
"""
module Collections

using InteractiveUtils

using EquationsOfState

export apply,
       EquationOfState,
       FiniteStrainEquationOfState,
    #    PolynomialEquationOfState,
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

An abstraction of equations of state, where `T` specifies the elements' type.
"""
abstract type EquationOfState{T} end

"""
    FiniteStrainEquationOfState{T} <: EquationOfState{T}

An abstraction of finite strain equations of state.
"""
abstract type FiniteStrainEquationOfState{T} <: EquationOfState{T} end

# struct PolynomialEquationOfState{T<:Real,N} <: EquationOfState{T,N}
#     data::NTuple{N,T}
#     function PolynomialEquationOfState{T,N}(args::NTuple{N,T}) where {T,N}
#         @assert N ≤ 10
#         new(args)
#     end
# end
# function PolynomialEquationOfState(args...)
#     T = Base.promote_typeof(args...)
#     PolynomialEquationOfState{T,length(args)}(args)
# end

"""
    Murnaghan(v0, b0, bp0, e0=0)

Create a Murnaghan equation of state. The elements' type will be handled automatically.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `bp0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0=0`: the energy of solid at zero pressure. By default is `0`.
"""
struct Murnaghan{T<:Real} <: EquationOfState{T}
    v0::T
    b0::T
    bp0::T
    e0::T
end
function Murnaghan(v0::Real, b0::Real, bp0::Real, e0::Real)
    T = Base.promote_typeof(v0, b0, bp0, e0)
    Murnaghan{T}(v0, b0, bp0, e0)
end
Murnaghan(v0, b0, bp0) = Murnaghan(v0, b0, bp0, 0)

"""
    BirchMurnaghan2nd(v0, b0, e0=0)

Create a Birch–Murnaghan 2nd order equation of state. The elements' type will be handled automatically.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `e0=0`: the energy of solid at zero pressure. By default is `0`.
"""
struct BirchMurnaghan2nd{T<:Real} <: FiniteStrainEquationOfState{T}
    v0::T
    b0::T
    e0::T
end
function BirchMurnaghan2nd(v0::Real, b0::Real, e0::Real)
    T = Base.promote_typeof(v0, b0, e0)
    BirchMurnaghan2nd{T}(v0, b0, e0)
end
BirchMurnaghan2nd(v0, b0) = BirchMurnaghan2nd(v0, b0, 0)

"""
    BirchMurnaghan3rd(v0, b0, bp0, e0=0)

Create a Birch–Murnaghan 3rd order equation of state. The elements' type will be handled automatically.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `bp0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0=0`: the energy of solid at zero pressure. By default is `0`.
"""
struct BirchMurnaghan3rd{T<:Real} <: FiniteStrainEquationOfState{T}
    v0::T
    b0::T
    bp0::T
    e0::T
end
function BirchMurnaghan3rd(v0::Real, b0::Real, bp0::Real, e0::Real)
    T = Base.promote_typeof(v0, b0, bp0, e0)
    BirchMurnaghan3rd{T}(v0, b0, bp0, e0)
end
BirchMurnaghan3rd(v0, b0, bp0) = BirchMurnaghan3rd(v0, b0, bp0, 0)

"""
    BirchMurnaghan4th(v0, b0, bp0, bpp0, e0=0)

Create a Birch–Murnaghan 4th order equation of state. The elements' type will be handled automatically.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `bp0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `bpp0`: the second-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0=0`: the energy of solid at zero pressure. By default is `0`.
"""
struct BirchMurnaghan4th{T<:Real} <: FiniteStrainEquationOfState{T}
    v0::T
    b0::T
    bp0::T
    bpp0::T
    e0::T
end
function BirchMurnaghan4th(v0::Real, b0::Real, bp0::Real, bpp0::Real, e0::Real)
    T = Base.promote_typeof(v0, b0, bp0, bpp0, e0)
    BirchMurnaghan4th{T}(v0, b0, bp0, bpp0, e0)
end
BirchMurnaghan4th(v0, b0, bp0, bpp0) = BirchMurnaghan4th(v0, b0, bp0, bpp0, 0)

"""
    PoirierTarantola2nd(v0, b0, e0=0)

Create a Poirier–Tarantola order equation of state. The elements' type will be handled automatically.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `e0=0`: the energy of solid at zero pressure. By default is `0`.
"""
struct PoirierTarantola2nd{T<:Real} <: FiniteStrainEquationOfState{T}
    v0::T
    b0::T
    e0::T
end
function PoirierTarantola2nd(v0::Real, b0::Real, e0::Real)
    T = Base.promote_typeof(v0, b0, e0)
    PoirierTarantola2nd{T}(v0, b0, e0)
end
PoirierTarantola2nd(v0, b0) = PoirierTarantola2nd(v0, b0, 0)

"""
    PoirierTarantola3rd(v0, b0, bp0, e0=0)

Create a Poirier–Tarantola 3rd order equation of state. The elements' type will be handled automatically.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `bp0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0=0`: the energy of solid at zero pressure. By default is `0`.
"""
struct PoirierTarantola3rd{T<:Real} <: FiniteStrainEquationOfState{T}
    v0::T
    b0::T
    bp0::T
    e0::T
end
function PoirierTarantola3rd(v0::Real, b0::Real, bp0::Real, e0::Real)
    T = Base.promote_typeof(v0, b0, bp0, e0)
    PoirierTarantola3rd{T}(v0, b0, bp0, e0)
end
PoirierTarantola3rd(v0, b0, bp0) = PoirierTarantola3rd(v0, b0, bp0, 0)

"""
    PoirierTarantola4th(v0, b0, bp0, bpp0, e0=0)

Create a Birch–Murnaghan 4th order equation of state. The elements' type will be handled automatically.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `bp0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `bpp0`: the second-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0=0`: the energy of solid at zero pressure. By default is `0`.
"""
struct PoirierTarantola4th{T<:Real} <: FiniteStrainEquationOfState{T}
    v0::T
    b0::T
    bp0::T
    bpp0::T
    e0::T
end
function PoirierTarantola4th(v0::Real, b0::Real, bp0::Real, bpp0::Real, e0::Real)
    T = Base.promote_typeof(v0, b0, bp0, bpp0, e0)
    PoirierTarantola4th{T}(v0, b0, bp0, bpp0, e0)
end
PoirierTarantola4th(v0, b0, bp0, bpp0) = PoirierTarantola4th(v0, b0, bp0, bpp0, 0)

"""
    Vinet(v0, b0, bp0, e0=0)

Create a Vinet equation of state. The elements' type will be handled automatically.

# Arguments
- `v0`: the volume of solid at zero pressure.
- `b0`: the bulk modulus of solid at zero pressure.
- `bp0`: the first-order pressure-derivative bulk modulus of solid at zero pressure.
- `e0=0`: the energy of solid at zero pressure. By default is `0`.
"""
struct Vinet{T<:Real} <: EquationOfState{T}
    v0::T
    b0::T
    bp0::T
    e0::T
end
function Vinet(v0::Real, b0::Real, bp0::Real, e0::Real)
    T = Base.promote_typeof(v0, b0, bp0, e0)
    Vinet{T}(v0, b0, bp0, e0)
end
Vinet(v0, b0, bp0) = Vinet(v0, b0, bp0, 0)

struct AntonSchmidt{T<:Real} <: EquationOfState{T}
    v0::T
    β::T
    n::T
    e∞::T
end
function AntonSchmidt(v0::Real, β::Real, n::Real, e∞::Real)
    T = Base.promote_typeof(v0, β, n, e∞)
    AntonSchmidt{T}(v0, β, n, e∞)
end
AntonSchmidt(v0, β, n) = AntonSchmidt(v0, β, n, 0)

struct BreenanStacey{T<:Real} <: EquationOfState{T}
    v0::T
    b0::T
    γ0::T
    e0::T
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
"""
    apply(EnergyForm(), eos::EquationOfState)

Return a function that can take a volume as a parameter, suitable for batch-applying.

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
```
"""
apply(form::EnergyForm, eos::EquationOfState) = v -> apply(form, eos, v)
"""
    apply(EnergyForm(), eos::Murnaghan, v::Real)

Return the energy of a `Murnaghan` equation of state on volume `v`.
"""
function apply(::EnergyForm, eos::Murnaghan, v::Real)
    v0, b0, bp0, e0 = collect(eos)

    x = bp0 - 1
    y = (v0 / v)^bp0
    return e0 + b0 / bp0 * v * (y / x + 1) - v0 * b0 / x
end
"""
    apply(EnergyForm(), eos::BirchMurnaghan2nd, v::Real)

Return the energy of a `BirchMurnaghan2nd` equation of state on volume `v`.
"""
function apply(::EnergyForm, eos::BirchMurnaghan2nd, v::Real)
    v0, b0, e0 = collect(eos)

    f = (cbrt(v0 / v)^2 - 1) / 2
    return e0 + 9 / 2 * b0 * v0 * f^2
end
"""
    apply(EnergyForm(), eos::BirchMurnaghan3rd, v::Real)

Return the energy of a `BirchMurnaghan3rd` equation of state on volume `v`.
"""
function apply(::EnergyForm, eos::BirchMurnaghan3rd, v::Real)
    v0, b0, bp0, e0 = collect(eos)

    eta = cbrt(v0 / v)
    xi = eta^2 - 1
    return e0 + 9 / 16 * b0 * v0 * xi^2 * (6 + bp0 * xi - 4 * eta^2)
end
"""
    apply(EnergyForm(), eos::BirchMurnaghan4th, v::Real)

Return the energy of a `BirchMurnaghan4th` equation of state on volume `v`.
"""
function apply(::EnergyForm, eos::BirchMurnaghan4th, v::Real)
    v0, b0, bp0, bpp0, e0 = collect(eos)

    f = (cbrt(v0 / v)^2 - 1) / 2
    h = b0 * bpp0 + bp0^2
    return e0 + 3 / 8 * v0 * b0 * f^2 * ((9h - 63bp0 + 143) * f^2 + 12 * (bp0 - 4) * f + 12)
end
"""
    apply(EnergyForm(), eos::PoirierTarantola2nd, v::Real)

Return the energy of a `PoirierTarantola2nd` equation of state on volume `v`.
"""
function apply(::EnergyForm, eos::PoirierTarantola2nd, v::Real)
    v0, b0, e0 = collect(eos)

    return e0 + b0 / 2 * v0 * log(v / v0)^(2 / 3)
end
"""
    apply(EnergyForm(), eos::PoirierTarantola3rd, v::Real)

Return the energy of a `PoirierTarantola3rd` equation of state on volume `v`.
"""
function apply(::EnergyForm, eos::PoirierTarantola3rd, v::Real)
    v0, b0, bp0, e0 = collect(eos)

    x = cbrt(v / v0)
    xi = -3 * log(x)
    return e0 + b0 / 6 * v0 * xi^2 * ((bp0 - 2) * xi + 3)
end
"""
    apply(EnergyForm(), eos::PoirierTarantola4th, v::Real)

Return the energy of a `PoirierTarantola4th` equation of state on volume `v`.
"""
function apply(::EnergyForm, eos::PoirierTarantola4th, v::Real)
    v0, b0, bp0, bpp0, e0 = collect(eos)

    x = cbrt(v / v0)
    xi = log(x)
    h = b0 * bpp0 + bp0^2
    return e0 + b0 / 24v0 * xi^2 * ((h + 3bp0 + 3) * xi^2 + 4 * (bp0 + 2) * xi + 12)
end
"""
    apply(EnergyForm(), eos::Vinet, v::Real)

Return the energy of a `Vinet` equation of state on volume `v`.
"""
function apply(::EnergyForm, eos::Vinet, v::Real)
    v0, b0, bp0, e0 = collect(eos)

    x = cbrt(v / v0)
    xi = 3 / 2 * (bp0 - 1)
    return e0 + 9b0 * v0 / xi^2 * (1 + (xi * (1 - x) - 1) * exp(xi * (1 - x)))
end
"""
    apply(EnergyForm(), eos::AntonSchmidt, v::Real)

Return the energy of a `AntonSchmidt` equation of state on volume `v`.
"""
function apply(::EnergyForm, eos::AntonSchmidt, v::Real)
    v0, β, n, e∞ = collect(eos)

    x = v / v0
    η = n + 1
    return e∞ + β * v0 / η * x^η * (log(x) - 1 / η)
end
# ============================= Energy evaluation ============================ #


# ============================================================================ #
#                              Pressure evaluation                             #
# ============================================================================ #
"""
    apply(PressureForm(), eos::EquationOfState)

Return a function that can take a volume as a parameter, suitable for batch-applying.

# Examples
```jldoctest
julia> using EquationsOfState, EquationsOfState.Collections

julia> f = apply(PressureForm(), Vinet(1, 2, 3));

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
apply(::PressureForm, eos::EquationOfState) = v -> apply(PressureForm(), eos, v)
"""
    apply(PressureForm(), eos::Murnaghan, v::Real)

Return the pressure of a `Murnaghan` equation of state on volume `v`.
"""
function apply(::PressureForm, eos::Murnaghan, v::Real)
    v0, b0, bp0 = collect(eos)

    return b0 / bp0 * ((v0 / v)^bp0 - 1)
end
"""
    apply(PressureForm(), eos::BirchMurnaghan2nd, v::Real)

Return the pressure of a `BirchMurnaghan2nd` equation of state on volume `v`.
"""
function apply(::PressureForm, eos::BirchMurnaghan2nd, v::Real)
    v0, b0 = collect(eos)

    f = ((v0 / v)^(2 / 3) - 1) / 2
    return 3b0 * f * (1 + 2f)^(5 / 2)
end
"""
    apply(PressureForm(), eos::BirchMurnaghan3rd, v::Real)

Return the pressure of a `BirchMurnaghan3rd` equation of state on volume `v`.
"""
function apply(::PressureForm, eos::BirchMurnaghan3rd, v::Real)
    v0, b0, bp0 = collect(eos)

    eta = (v0 / v)^(1 / 3)
    return 3 / 2 * b0 * (eta^7 - eta^5) * (1 + 3 / 4 * (bp0 - 4) * (eta^2 - 1))
end
"""
    apply(PressureForm(), eos::BirchMurnaghan4th, v::Real)

Return the pressure of a `BirchMurnaghan4th` equation of state on volume `v`.
"""
function apply(::PressureForm, eos::BirchMurnaghan4th, v::Real)
    v0, b0, bp0, bpp0 = collect(eos)

    f = ((v0 / v)^(2 / 3) - 1) / 2
    h = b0 * bpp0 + bp0^2
    return b0 / 2 * (2f + 1)^(5 / 2) * ((9h - 63bp0 + 143) * f^2 + 9 * (bp0 - 4) * f + 6)
end
"""
    apply(PressureForm(), eos::PoirierTarantola2nd, v::Real)

Return the pressure of a `PoirierTarantola2nd` equation of state on volume `v`.
"""
function apply(::PressureForm, eos::PoirierTarantola2nd, v::Real)
    v0, b0 = collect(eos)

    x = (v / v0)^(1 / 3)
    return -b0 / x * log(x)
end
"""
    apply(PressureForm(), eos::PoirierTarantola3rd, v::Real)

Return the pressure of a `PoirierTarantola3rd` equation of state on volume `v`.
"""
function apply(::PressureForm, eos::PoirierTarantola3rd, v::Real)
    v0, b0, bp0 = collect(eos)

    x = v / v0
    xi = log(x)
    return -b0 * xi / 2x * ((bp0 - 2) * xi - 2)
end
"""
    apply(PressureForm(), eos::PoirierTarantola4th, v::Real)

Return the pressure of a `PoirierTarantola4th` equation of state on volume `v`.
"""
function apply(::PressureForm, eos::PoirierTarantola4th, v::Real)
    v0, b0, bp0, bpp0 = collect(eos)

    x = (v / v0)^(1 / 3)
    xi = log(x)
    h = b0 * bpp0 + bp0^2
    return -b0 * xi / 6 / x * ((h + 3bp0 + 3) * xi^2 + 3 * (bp0 + 6) * xi + 6)
end
"""
    apply(PressureForm(), eos::Vinet, v::Real)

Return the pressure of a `Vinet` equation of state on volume `v`.
"""
function apply(::PressureForm, eos::Vinet, v::Real)
    v0, b0, bp0 = collect(eos)

    x = (v / v0)^(1 / 3)
    xi = 3 / 2 * (bp0 - 1)
    return 3b0 / x^2 * (1 - x) * exp(xi * (1 - x))
end
"""
    apply(PressureForm(), eos::AntonSchmidt, v::Real)

Return the pressure of a `AntonSchmidt` equation of state on volume `v`.
"""
function apply(::PressureForm, eos::AntonSchmidt, v::Real)
    v0, β, n = collect(eos)

    x = v / v0
    return -β * x^n * log(x)
end
"""
    apply(PressureForm(), eos::BreenanStacey, v::Real)

Return the pressure of a `BreenanStacey` equation of state on volume `v`.
"""
function apply(::PressureForm, eos::BreenanStacey, v::Real)
    v0, b0, γ0 = collect(eos)

    x = v0 / v
    return b0 / 2 / γ0 * x^(4 / 3) * (exp(2γ0 * (1 - x)) - 1)
end
# ============================ Pressure evaluation =========================== #


# ============================================================================ #
#                            Bulk modulus evaluation                           #
# ============================================================================ #
"""
    apply(BulkModulusForm(), eos::EquationOfState)

Return a function that can take a volume as a parameter, suitable for batch-applying.

# Examples
```jldoctest
julia> using EquationsOfState, EquationsOfState.Collections

julia> f = apply(BulkModulusForm(), BirchMurnaghan3rd(1, 2, 3));

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
apply(::BulkModulusForm, eos::EquationOfState) = v -> apply(BulkModulusForm(), eos, v)
"""
    apply(BulkModulusForm(), eos::BirchMurnaghan2nd, v::Real)

Return the bulk modulus of a `BirchMurnaghan2nd` equation of state on volume `v`.
"""
function apply(::BulkModulusForm, eos::BirchMurnaghan2nd, v::Real)
    v0, b0 = collect(eos)

    f = ((v0 / v)^(2 / 3) - 1) / 2
    return b0 * (7f + 1) * (2f + 1)^(5 / 2)
end
"""
    apply(BulkModulusForm(), eos::BirchMurnaghan3rd, v::Real)

Return the bulk modulus of a `BirchMurnaghan3rd` equation of state on volume `v`.
"""
function apply(::BulkModulusForm, eos::BirchMurnaghan3rd, v::Real)
    v0, b0, bp0 = collect(eos)

    f = ((v0 / v)^(2 / 3) - 1) / 2
    return b0 / 2 * (2f + 1)^(5 / 2) * ((27 * f^2 + 6f) * (bp0 - 4) - 4f + 2)
end
"""
    apply(BulkModulusForm(), eos::BirchMurnaghan4th, v::Real)

Return the bulk modulus of a `BirchMurnaghan4th` equation of state on volume `v`.
"""
function apply(::BulkModulusForm, eos::BirchMurnaghan4th, v::Real)
    v0, b0, bp0, bpp0 = collect(eos)

    f = ((v0 / v)^(2 / 3) - 1) / 2
    h = b0 * bpp0 + bp0^2
    return b0 / 6 * (2f + 1)^(5 / 2) *
           ((99h - 693bp0 + 1573) * f^3 + (27h - 108bp0 + 105) * f^2 + 6f * (3bp0 - 5) + 6)
end
"""
    apply(BulkModulusForm(), eos::PoirierTarantola2nd, v::Real)

Return the bulk modulus of a `PoirierTarantola2nd` equation of state on volume `v`.
"""
function apply(::BulkModulusForm, eos::PoirierTarantola2nd, v::Real)
    v0, b0 = collect(eos)

    x = (v / v0)^(1 / 3)
    return b0 / x * (1 - log(x))
end
"""
    apply(BulkModulusForm(), eos::PoirierTarantola3rd, v::Real)

Return the bulk modulus of a `PoirierTarantola3rd` equation of state on volume `v`.
"""
function apply(::BulkModulusForm, eos::PoirierTarantola3rd, v::Real)
    v0, b0, bp0 = collect(eos)

    x = v / v0
    xi = log(x)
    return -b0 / 2x * (((bp0 - 2) * xi + 2 - 2bp0) * xi + 2)
end
"""
    apply(BulkModulusForm(), eos::PoirierTarantola4th, v::Real)

Return the bulk modulus of a `PoirierTarantola4th` equation of state on volume `v`.
"""
function apply(::BulkModulusForm, eos::PoirierTarantola4th, v::Real)
    v0, b0, bp0, bpp0 = collect(eos)

    x = (v / v0)^(1 / 3)
    xi = log(x)
    h = b0 * bpp0 + bp0^2
    return -b0 / (6x) *
           ((h + 3bp0 + 3) * xi^3 - 3 * xi^2 * (h + 2bp0 + 1) - 6xi * (bp0 + 1) - 6)
end
"""
    apply(BulkModulusForm(), eos::Vinet, v::Real)

Return the bulk modulus of a `Vinet` equation of state on volume `v`.
"""
function apply(::BulkModulusForm, eos::Vinet, v::Real)
    v0, b0, bp0 = collect(eos)

    x = (v / v0)^(1 / 3)
    xi = 3 / 2 * (bp0 - 1)
    return -b0 / (2 * x^2) * (3x * (x - 1) * (bp0 - 1) + 2 * (x - 2)) * exp(-xi * (x - 1))
end
"""
    apply(BulkModulusForm(), eos::AntonSchmidt, v::Real)

Return the bulk modulus of a `AntonSchmidt` equation of state on volume `v`.
"""
function apply(::BulkModulusForm, eos::AntonSchmidt, v::Real)
    v0, β, n = collect(eos)

    x = v / v0
    return β * x^n * (1 + n * log(x))
end
# ========================== Bulk modulus evaluation ========================= #


# ============================================================================ #
#                                 Miscellaneous                                #
# ============================================================================ #
function allsubtypes(T::Type, types = Type[])::Vector{Type}
    for S in subtypes(T)
        types = allsubtypes(S, push!(types, S))
    end
    types
end

nonabstract(T::Type)::Vector{Type} = filter(!isabstracttype, allsubtypes(T))

for E in nonabstract(EquationOfState)
    eval(quote
        similar_type(::Type{A}, ::Type{T}) where {A<:$E,T} = $E{T}
    end)
end

Base.iterate(eos::EquationOfState) = getfield(eos, 1), 2
function Base.iterate(eos::EquationOfState, state::Integer)
    state > length(eos) && return nothing
    return getfield(eos, state), state + 1
end # function Base.iterate
Base.length(eos::EquationOfState) = nfields(eos)
Base.size(eos::EquationOfState) = (length(eos),)
Base.eltype(::Type{<:EquationOfState{T}}) where {T} = T

Base.isapprox(a::T, b::T; kwargs...) where {T<:EquationOfState} = isapprox(collect(a), collect(b); kwargs...)
# Base.getindex(eos::PolynomialEquationOfState{T,N}, index::Int64) where {T,N} = getindex(eos.data, index)
# =============================== Miscellaneous ============================== #

end
