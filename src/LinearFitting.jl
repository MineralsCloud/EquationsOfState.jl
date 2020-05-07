"""
This module provides some linear fitting methods.
"""
module LinearFitting

using LinearAlgebra: dot
using Polynomials: degree, coeffs
using Polynomials.PolyCompat: polyfit, polyder, Poly

export FiniteStrain, Eulerian, Lagrangian, Natural, Infinitesimal
export energy_strain_expansion,
    energy_strain_derivative,
    strain_volume_derivative,
    energy_volume_expansion,
    energy_volume_derivatives,
    energy_volume_nth_derivative

abstract type FiniteStrain end
struct Eulerian <: FiniteStrain end
struct Lagrangian <: FiniteStrain end
struct Natural <: FiniteStrain end
struct Infinitesimal <: FiniteStrain end

(::EulerianStrain)(v0, v) = (cbrt(v0 / v)^2 - 1) / 2
(::LagrangianStrain)(v0, v) = (cbrt(v / v0)^2 - 1) / 2
(::NaturalStrain)(v0, v) = log(v / v0) / 3
(::InfinitesimalStrain)(v0, v) = 1 - cbrt(v0 / v)

energy_strain_expansion(f::Vector{<:Real}, e::Vector{<:Real}, n::Int) = polyfit(f, e, n)

energy_strain_derivative(p::Poly, deg) = polyder(p, deg)

function strain_volume_derivative(s::Eulerian, v0, v, deg)
    if deg == 1
        return -1 / (3v) * cbrt(v0 / v)^2
    else  # Recursion
        return -(3deg + 2) / (3v) * strain_volume_derivative(s, v0, v, deg - 1)
    end
end
function strain_volume_derivative(s::Lagrangian, v0, v, deg)
    if deg == 1
        return -1 / (3v) * cbrt(v / v0)^2
    else  # Recursion
        return -(3deg - 2) / (3v) * strain_volume_derivative(s, v0, v, deg - 1)
    end
end
function strain_volume_derivative(s::Natural, v0, v, deg)
    if deg == 1
        return 1 / (3v)
    else  # Recursion
        return -deg / v * strain_volume_derivative(s, v0, v, deg - 1)
    end
end
function strain_volume_derivative(s::Infinitesimal, v0, v, deg)
    if deg == 1
        return (1 - s(v0, v))^4 / 3 / v0
    else  # Recursion
        return -(3deg + 1) / (3v) * strain_volume_derivative(s, v0, v, deg - 1)
    end
end

function energy_volume_expansion(s::FiniteStrain, v0, v, p::Poly, highest_order = degree(p))
    # The zeroth order value plus values from the first to the ``highest_order`.
    return p(v) + dot(
        energy_volume_derivatives(s, v0, v, p, highest_order),
        s(v0, v) .^ collect(1:highest_order),
    )
end

function energy_volume_derivatives(s::FiniteStrain, v0, v, p::Poly, highest_order)
    if 0 ≤ highest_order ≤ degree(p)
        x = 1:highest_order
        strain_derivatives = map(m -> strain_volume_derivative(s, v0, v, m), x)
        energy_derivatives = map(f -> f(v), map(m -> energy_strain_derivative(p, m), x))
        return map(x) do m
            energy_volume_nth_derivative(m)(strain_derivatives, energy_derivatives)
        end
    else
        throw(DomainError("The `highest_order` must be within 0 to $(degree(p))!"))
    end
end

function energy_volume_nth_derivative(deg)
    function (f::Vector{<:Real}, e::Vector{<:Real})
        if deg == 1
            return e[1] * f[1]
        elseif deg == 2
            return e[2] * f[1]^2 + e[1] * f[1]
        elseif deg == 3
            return e[3] * f[1]^3 + 3 * f[1] * f[2] * e[2] + e[1] * f[3]
        elseif deg == 4
            return e[4] * f[1]^4 +
                   6 * f[1]^2 * f[2] * e[3] +
                   (4 * f[1] * f[3] + 3 * f[3]^2) * e[2] +
                   e[1] * f[3]
        else
            error("Expansion is not defined at order = $(deg)!")
        end
    end
end

end
