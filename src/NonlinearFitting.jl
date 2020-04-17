"""
This module provides `lsqfit` methods for fitting an equation of state
with(out) units.
"""
module NonlinearFitting

using ConstructionBase: constructorof
using IterTools: fieldvalues
using LsqFit: curve_fit
using Unitful: AbstractQuantity, upreferred, ustrip, unit

using ..Collections: PhysicalProperty, EquationOfState

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
    f::Function,
    xdata::AbstractVector{<:Real},
    ydata::AbstractVector{<:Real};
    debug = false,
    kwargs...,
)
    T = constructorof(typeof(f.eos))  # Get the `UnionAll` type
    model = (x, p) -> map(T(p...)(f.prop), x)
    fitted = curve_fit(
        model,
        float(xdata),  # Convert `xdata` elements to floats
        float(ydata),  # Convert `ydata` elements to floats
        float.(fieldvalues(f.eos));  # TODO: What if these floats are different types?
        kwargs...,
    )
    return debug ? fitted : T(fitted.param...)
end  # function lsqfit
function lsqfit(
    f::Function,
    xdata::AbstractVector{<:AbstractQuantity},
    ydata::AbstractVector{<:AbstractQuantity};
    kwargs...,
)
    T = constructorof(typeof(f.eos))  # Get the `UnionAll` type
    values = fieldvalues(f.eos)
    original_units = unit.(values)  # Keep a record of `eos`'s units
    g = x -> map(ustrip ∘ upreferred, x)  # Convert to preferred units and strip the unit
    trial_params = g.(values)
    result = lsqfit(T(trial_params...)(f.prop), g.(xdata), g.(ydata); kwargs...)
    if result isa EquationOfState  # i.e., if `debug = false` and no error is thrown
        return T((
            x * upreferred(u) |> u for (x, u) in zip(fieldvalues(result), original_units)
        )...)  # Convert back to original `eos`'s units
    else
        return result
    end
end  # function lsqfit

end
