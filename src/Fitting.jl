module Fitting

using ConstructionBase: constructorof
using IterTools: fieldvalues
using LsqFit: curve_fit
using Polynomials: Polynomial, fit, derivative, coeffs
using PolynomialRoots: roots
using Unitful: AbstractQuantity, NoDims, upreferred, ustrip, unit, dimension, @u_str

using ..Collections: PhysicalProperty, EquationOfState, PolynomialEOS

export linfit, nonlinfit

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
        error("no real local minima found!")  # For some polynomials, could be all complex
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
    p0 = Tuple(derivative(poly, n)(v0) / factorial(n) for n in 1:deg)
    return PolynomialEOS(v0, p0, e0)
end # function linfit

"""
    nonlinfit(eos(prop), volumes, ydata; kwargs...)

Fit an equation of state using least-squares fitting method (with the Levenberg-Marquardt algorithm).

If `eos`, `volumes` and `ydata` are all unitless, `volumes` must have the same unit as
`eos.v0`, and `ydata`, if are energies, must have the same unit as `eos.e0`. The
`eos` parameter `b0` must have the unit of `e0 / v0`, and `b′′0` must have the unit of
`v0 / e0`, etc. Use with care. Better to use the `Unitful` version.

# Arguments
- `eos::EquationOfState`: a trial equation of state. If it has units, `volumes` and `ydata` must also have.
- `prop::PhysicalProperty`: a `PhysicalProperty` instance. If `Energy`, fit ``E(V)``; if `Pressure`, fit ``P(V)``; if `BulkModulus`, fit ``B(V)``.
- `volumes`: an array of volumes (``V``), with(out) units.
- `ydata`: an array of energies (``E``), pressures (``P``), or bulk moduli (``B``), with(out) units. It must be consistent with `prop`.
- `kwargs`: the rest keyword arguments are the same as that of `LsqFit.curve_fit`. See its [documentation](https://github.com/JuliaNLSolvers/LsqFit.jl/blob/master/README.md) and [tutorial](https://julianlsolvers.github.io/LsqFit.jl/latest/tutorial/).
"""
function nonlinfit(f, volumes, ydata; kwargs...)
    eos, property = fieldvalues(f)
    T = constructorof(typeof(eos))  # Get the `UnionAll` type
    params, volumes, ydata = _preprocess(float(eos), float(volumes), float(ydata))
    model = (x, p) -> map(T(p...)(property), x)
    fit = curve_fit(model, volumes, ydata, params; kwargs...)
    return _postprocess(T(fit.param...), eos)
end # function nonlinfit

_preprocess(eos, xdata, ydata) = (collect(fieldvalues(eos)), xdata, ydata)
function _preprocess(
    eos::EquationOfState{<:AbstractQuantity},
    xdata::AbstractArray{<:AbstractQuantity},
    ydata::AbstractArray{<:AbstractQuantity},
)
    values = fieldvalues(eos)
    original_units = unit.(values)  # Keep a record of `eos`'s units
    return _ustrip.(values), _ustrip.(xdata), _ustrip.(ydata)  # Convert to preferred then and strip the unit
end # function _preprocess

_postprocess(eos, args...) = eos
function _postprocess(eos, trial_eos::EquationOfState{<:AbstractQuantity})
    T = constructorof(typeof(trial_eos))  # Get the `UnionAll` type
    original_units = unit.(fieldvalues(trial_eos))  # Keep a record of `eos`'s units
    return T((
        x * _upreferred(dimension(u)) |> u for
        (x, u) in zip(fieldvalues(eos), original_units)
    )...)  # Convert back to original `eos`'s units
end # function _postprocess

_ustrip(quantity) = _ustrip(_upreferred(dimension(quantity)), quantity)
_ustrip(unit, quantity) = ustrip(unit, quantity)
_ustrip(::Int, quantity) = ustrip(quantity)

_upreferred(::typeof(dimension(u"1"))) = 1
_upreferred(::typeof(dimension(u"J"))) = u"eV"
_upreferred(::typeof(dimension(u"m^3"))) = u"angstrom^3"
_upreferred(::typeof(dimension(u"Pa"))) = u"eV/angstrom^3"
_upreferred(::typeof(dimension(u"1/Pa"))) = u"angstrom^3/eV"
_upreferred(::typeof(dimension(u"1/Pa^2"))) = u"angstrom^6/eV^2"

Base.float(eos::EquationOfState) =
    constructorof(typeof(eos))(map(float, fieldvalues(eos))...)

end # module Fitting
