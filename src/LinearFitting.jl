"""
This module provides some linear fitting methods.
"""
module LinearFitting

using Polynomials: Polynomial, fit, derivative, coeffs
using PolynomialRoots: roots

using ..Collections: PolynomialEOS

export linfit

_islocalminimum(f, x, δx) = f(x) < f(x - δx) && f(x) < f(x + δx)

function _findlocalminima(f, xs)
    f′ = derivative(f, 1)
    δx = minimum(diff(xs)) / 10
    localminima = eltype(xs)[]
    for x in real(filter(isreal, roots(coeffs(f′))))  # Complex volumes are meaningless
        if _islocalminimum(f, x, δx)
            push!(localminima, x)
        end
    end
    return localminima
end # function _findlocalminima

function _findglobalminimum(f, localminima)
    # https://stackoverflow.com/a/21367608/3260253
    if isempty(localminima)
        error("no local minima found!")
    else
        f0, i = findmin(f.(localminima))
        x0 = localminima[i]
        return x0, f0
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
