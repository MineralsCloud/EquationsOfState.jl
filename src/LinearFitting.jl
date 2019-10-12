"""
This module provides some linear fitting methods.
"""
module LinearFitting

using LinearAlgebra: dot
using Polynomials: polyder, polyfit, degree, coeffs, Poly
using MLStyle: @match

using ..FiniteStrains

export energy_strain_expansion,
       energy_strain_derivative,
       strain_volume_derivative,
       energy_volume_expansion,
       energy_volume_derivatives,
       energy_volume_derivative_at_order

energy_strain_expansion(f::Vector{<:Real}, e::Vector{<:Real}, n::Int)::Poly =
    polyfit(f, e, n)

energy_strain_derivative(p::Poly, m::Int)::Poly = polyder(p, m)

function strain_volume_derivative(s::EulerianStrain, v0::Real, v::Real, m::Int)
    m == 1 && return -1 / (3v) * cbrt(v0 / v)^2
    -(3m + 2) / (3v) * strain_volume_derivative(s, v0, v, m - 1)
end
function strain_volume_derivative(s::LagrangianStrain, v0::Real, v::Real, m::Int)
    m == 1 && return -1 / (3v) * cbrt(v / v0)^2
    -(3m - 2) / (3v) * strain_volume_derivative(s, v0, v, m - 1)
end
function strain_volume_derivative(s::NaturalStrain, v0::Real, v::Real, m::Int)
    m == 1 && return 1 / (3v)
    -m / v * strain_volume_derivative(s, v0, v, m - 1)
end
function strain_volume_derivative(s::InfinitesimalStrain, v0::Real, v::Real, m::Int)
    m == 1 && return (1 - getstrain(s, v0, v))^4 / 3 / v0
    -(3m + 1) / (3v) * strain_volume_derivative(s, v0, v, m - 1)
end

function energy_volume_expansion(
    s::FiniteStrain,
    v0::Real,
    v::Real,
    p::Poly,
    highest_order::Int = degree(p),
)
    # The zeroth order value plus values from the first to the ``highest_order`.
    p(v) + dot(
        energy_volume_derivatives(s, v0, v, p, highest_order),
        getstrain(s, v0, v) .^ collect(1:highest_order),
    )
end

function energy_volume_derivatives(
    s::FiniteStrain,
    v0::Real,
    v::Real,
    p::Poly,
    highest_order::Int,
)
    0 ≤ highest_order ≤ degree(p) ? (x = 1:highest_order) :
    throw(DomainError("The `highest_order` must be within 0 to $(degree(p))!"))
    strain_derivatives = map(m -> strain_volume_derivative(s, v0, v, m), x)
    energy_derivatives = map(f -> f(v), map(m -> energy_strain_derivative(p, m), x))
    map(
        m -> energy_volume_derivative_at_order(m)(strain_derivatives, energy_derivatives),
        x,
    )
end

function energy_volume_derivative_at_order(m::Int)::Function
    function (f::Vector{<:Real}, e::Vector{<:Real})
        @match m begin
            1 => e[1] * f[1]
            2 => e[2] * f[1]^2 + e[1] * f[1]
            3 => e[3] * f[1]^3 + 3 * f[1] * f[2] * e[2] + e[1] * f[3]
            4 => e[4] * f[1]^4 + 6 * f[1]^2 * f[2] * e[3] +
                 (4 * f[1] * f[3] + 3 * f[3]^2) * e[2] + e[1] * f[3]
            _ => error("Expansion is not defined at order = $(m)!")
        end
    end
end

end
