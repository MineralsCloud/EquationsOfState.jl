"""
This module provides some linear fitting methods.
"""
module LinearFitting

using Polynomials: Polynomial, fit, derivative, coeffs
using PolynomialRoots: roots

using ..Collections: PolynomialEOS

export linfit

_islocalminimum(y, x, δx) = y(x) < y(x - δx) && y(x) < y(x + δx)

function _findlocalminima(y, xs)
    y′ = derivative(y, 1)
    δx = minimum(diff(xs)) / 10
    localminima = eltype(xs)[]
    for x in real(filter(isreal, roots(coeffs(y′))))  # Complex volumes are meaningless
        if _islocalminimum(y, x, δx)
            push!(localminima, x)
        end
    end
    return localminima
end # function _findlocalminima

function _findglobalminimum(y, localminima)
    # https://stackoverflow.com/a/21367608/3260253
    if isempty(localminima)
        error("no local minima found!")
    else
        y0, i = findmin(y.(localminima))
        x0 = localminima[i]
        return x0, y0
    end
end # function _findglobalminimum

function linfit(volumes, energies, deg = 3)
    poly = fit(volumes, energies, deg)
    localminima = _findlocalminima(poly, volumes)
    v0, e0 = _findglobalminimum(poly, localminima)
    bs = Tuple(derivative(poly, n)(v0) / factorial(n) for n in 1:deg)
    return PolynomialEOS(v0, bs, e0)
end # function linfit

end
