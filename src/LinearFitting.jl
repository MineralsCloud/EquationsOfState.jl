"""
This module provides some linear fitting methods.
"""
module LinearFitting

using Polynomials: Polynomial, fit, derivative, roots, coeffs

export FiniteStrain, Eulerian, Lagrangian, Natural, Infinitesimal, linearfit

abstract type FiniteStrain end
struct Eulerian <: FiniteStrain end
struct Lagrangian <: FiniteStrain end
struct Natural <: FiniteStrain end
struct Infinitesimal <: FiniteStrain end

struct StrainFromVolume{T<:FiniteStrain}
    v0::Any
end
(x::StrainFromVolume{Eulerian})(v) = (cbrt(x.v0 / v)^2 - 1) / 2
(x::StrainFromVolume{Lagrangian})(v) = (cbrt(v / x.v0)^2 - 1) / 2
(x::StrainFromVolume{Natural})(v) = log(v / x.v0) / 3
(x::StrainFromVolume{Infinitesimal})(v) = 1 - cbrt(x.v0 / v)

struct VolumeFromStrain{T<:FiniteStrain}
    v0::Any
end
(x::VolumeFromStrain{Eulerian})(f) = x.v0 / (2f + 1)^(3 / 2)
(x::VolumeFromStrain{Lagrangian})(f) = x.v0 * (2f + 1)^(3 / 2)
(x::VolumeFromStrain{Natural})(f) = x.v0 * exp(3f)
(x::VolumeFromStrain{Infinitesimal})(f) = x.v0 / (1 - f)^3

_islocalminimum(f, x, δx) = f(x) < f(x - δx) && f(x) < f(x + δx)

function _findglobalminimum(f, localminima, δx)
    if length(localminima) == 0
        error("no volume minimizes the energy!")
    else
        f0, i = findmin(f.(localminima))
        x0 = localminima[i]
        return x0, f0
    end
end # function _findglobalminimum

function linearfit(volumes, energies, deg = 3)
    poly = fit(volumes, energies, deg)
    poly1d = derivative(poly, 1)
    δx = minimum(diff(volumes)) / 10
    localminima = eltype(volumes)[]
    for x in real(filter(isreal, roots(poly1d)))  # Complex volume is meaningless
        if _islocalminimum(poly, x, δx)
            push!(localminima, x)
        end
    end
    v0, e0 = _findglobalminimum(poly, localminima, δx)
end # function linearfit

Base.inv(x::StrainFromVolume{T}) where {T} = VolumeFromStrain{T}(x.v0)
Base.inv(x::VolumeFromStrain{T}) where {T} = StrainFromVolume{T}(x.v0)

Base.:∘(a::StrainFromVolume{T}, b::VolumeFromStrain{T}) where {T} =
    a.v0 == b.v0 ? identity : error("operation undefined!")
Base.:∘(a::VolumeFromStrain{T}, b::StrainFromVolume{T}) where {T} =
    a.v0 == b.v0 ? identity : error("operation undefined!")

end
