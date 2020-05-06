"""
This module provides some linear fitting methods.
"""
module LinearFitting

using LinearAlgebra: dot
using Polynomials: degree, coeffs
using Polynomials.PolyCompat: polyfit, polyder, Poly

export FiniteStrain, EulerianStrain, LagrangianStrain, NaturalStrain, InfinitesimalStrain
export energy_strain_expansion,
    energy_strain_derivative,
    strain_volume_derivative,
    energy_volume_expansion,
    energy_volume_derivatives,
    energy_volume_derivative_at_order

abstract type FiniteStrain end

struct EulerianStrain <: FiniteStrain end
struct LagrangianStrain <: FiniteStrain end
struct NaturalStrain <: FiniteStrain end
struct InfinitesimalStrain <: FiniteStrain end

(::EulerianStrain)(v0, v) = (cbrt(v0 / v)^2 - 1) / 2
(::LagrangianStrain)(v0, v) = (cbrt(v / v0)^2 - 1) / 2
(::NaturalStrain)(v0, v) = log(v / v0) / 3
(::InfinitesimalStrain)(v0, v) = 1 - cbrt(v0 / v)

energy_strain_expansion(f::Vector{<:Real}, e::Vector{<:Real}, n::Int) = polyfit(f, e, n)

energy_strain_derivative(p::Poly, deg) = polyder(p, deg)

function strain_volume_derivative(s::EulerianStrain, v0, v, deg)
    deg == 1 && return -1 / (3v) * cbrt(v0 / v)^2
    -(3deg + 2) / (3v) * strain_volume_derivative(s, v0, v, deg - 1)
end
function strain_volume_derivative(s::LagrangianStrain, v0, v, deg)
    deg == 1 && return -1 / (3v) * cbrt(v / v0)^2
    -(3deg - 2) / (3v) * strain_volume_derivative(s, v0, v, deg - 1)
end
function strain_volume_derivative(s::NaturalStrain, v0, v, deg)
    deg == 1 && return 1 / (3v)
    -deg / v * strain_volume_derivative(s, v0, v, deg - 1)
end
function strain_volume_derivative(s::InfinitesimalStrain, v0, v, deg)
    deg == 1 && return (1 - s(v0, v))^4 / 3 / v0
    -(3deg + 1) / (3v) * strain_volume_derivative(s, v0, v, deg - 1)
end

function energy_volume_expansion(s::FiniteStrain, v0, v, p::Poly, highest_order = degree(p))
    # The zeroth order value plus values from the first to the ``highest_order`.
    p(v) + dot(
        energy_volume_derivatives(s, v0, v, p, highest_order),
        s(v0, v) .^ collect(1:highest_order),
    )
end

function energy_volume_derivatives(s::FiniteStrain, v0, v, p::Poly, highest_order)
    0 ≤ highest_order ≤ degree(p) ? (x = 1:highest_order) :
    throw(DomainError("The `highest_order` must be within 0 to $(degree(p))!"))
    strain_derivatives = map(m -> strain_volume_derivative(s, v0, v, m), x)
    energy_derivatives = map(f -> f(v), map(m -> energy_strain_derivative(p, m), x))
    map(
        m -> energy_volume_derivative_at_order(m)(strain_derivatives, energy_derivatives),
        x,
    )
end

function energy_volume_derivative_at_order(deg)
    function (f::Vector{<:Real}, e::Vector{<:Real})
        if deg == 1
            e[1] * f[1]
        elseif deg == 2
            e[2] * f[1]^2 + e[1] * f[1]
        elseif deg == 3
            e[3] * f[1]^3 + 3 * f[1] * f[2] * e[2] + e[1] * f[3]
        elseif deg == 4
            e[4] * f[1]^4 +
            6 * f[1]^2 * f[2] * e[3] +
            (4 * f[1] * f[3] + 3 * f[3]^2) * e[2] +
            e[1] * f[3]
        else
            error("Expansion is not defined at order = $(deg)!")
        end
    end
end

end
