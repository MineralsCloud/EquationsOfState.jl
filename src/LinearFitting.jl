"""
This module provides some linear fitting methods.
"""
module LinearFitting

using Polynomials: Polynomial, fit, derivative, roots, coeffs

using ..Collections: PolynomialEOS

export linearfit

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
    bs = Tuple(derivative(poly, n)(v0) / factorial(n) for n in 1:deg)
    return PolynomialEOS(v0, bs, e0)
end # function linearfit

end
