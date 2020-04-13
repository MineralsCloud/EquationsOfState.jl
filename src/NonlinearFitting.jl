"""
This module provides `lsqfit` methods for fitting an equation of state
with(out) units.
"""
module NonlinearFitting

using ConstructionBase: constructorof
using LsqFit: curve_fit
using Unitful: AbstractQuantity, upreferred, ustrip, unit

using ..Collections: PhysicalProperty, EquationOfState, EquationOnVolume, fieldvalues

export lsqfit

"""
    lsqfit(eos(prop), xdata, ydata; debug = false, kwargs...)

Fit an equation of state using least-squares fitting method (with the Levenberg-Marquardt algorithm).

# Arguments
- `eos::EquationOfState`: a trial equation of state. If it has units, `xdata` and `ydata` must also have.
- `prop::PhysicalProperty`: a `PhysicalProperty` instance. If `Energy`, fit ``E(V)``; if `Pressure`, fit ``P(V)``; if `BulkModulus`, fit ``B(V)``.
- `xdata::AbstractVector`: a vector of volumes (``V``), with(out) units.
- `ydata::AbstractVector`: a vector of energies (``E``), pressures (``P``), or bulk moduli (``B``), with(out) units. It must be consistent with `prop`.
- `debug::Bool=false`: if `true`, then an `LsqFit.LsqFitResult` is returned, containing estimated Jacobian, residuals, etc.; if `false`, a fitted `EquationOfState` is returned. The default value is `false`.
- `kwargs`: the rest keyword arguments are the same as that of `LsqFit.curve_fit`. See its [documentation](https://github.com/JuliaNLSolvers/LsqFit.jl/blob/master/README.md)
    and [tutorial](https://julianlsolvers.github.io/LsqFit.jl/latest/tutorial/).
"""
function lsqfit(
    (eos, prop)::EquationOnVolume{<:Real},
    xdata::AbstractVector{<:Real},
    ydata::AbstractVector{<:Real};
    debug = false,
    kwargs...,
)
    E = constructorof(typeof(eos))  # Get the `UnionAll` type
    model = (x, p) -> map(E(p...)(prop), x)
    fitted = curve_fit(
        model,
        float(xdata),  # Convert `xdata` elements to floats
        float(ydata),  # Convert `ydata` elements to floats
        float(fieldvalues(eos));  # TODO: What if these floats are different types?
        kwargs...,
    )
    return debug ? fitted : E(fitted.param...)
end  # function lsqfit
function lsqfit(
    (eos, prop)::EquationOnVolume{<:AbstractQuantity},
    xdata::AbstractVector{<:AbstractQuantity},
    ydata::AbstractVector{<:AbstractQuantity};
    kwargs...,
)
    E = constructorof(typeof(eos))  # Get the `UnionAll` type
    values = fieldvalues(eos)
    original_units = unit.(values)  # Keep a record of `eos`'s units
    f = x -> map(ustrip ∘ upreferred, x)  # Convert to preferred units and strip the unit
    trial_params = f.(values)
    result = lsqfit(E(trial_params...)(prop), f.(xdata), f.(ydata); kwargs...)
    if result isa EquationOfState  # i.e., if `debug = false` and no error is thrown
        data = fieldvalues(result)
        # Convert back to original `eos`'s units
        return E((data[i] * upreferred(u) |> u for (i, u) in enumerate(original_units))...)
    end
    return result
end  # function lsqfit

end
