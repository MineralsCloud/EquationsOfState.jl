"""
This module provides `lsqfit` methods for fitting an equation of state
with(out) units.
"""
module NonlinearFitting

using ConstructionBase: constructorof
using LsqFit: curve_fit
using Unitful: AbstractQuantity, upreferred, ustrip, unit

import ..EquationForm
using ..Collections

export lsqfit

"""
    lsqfit(form, eos, xdata, ydata; debug = false, kwargs...)

Fit an equation of state using least-squares fitting method (with the Levenberg-Marquardt algorithm).

# Arguments
- `form::EquationForm`: an `EquationForm` instance. If `EnergyForm`, fit ``E(V)``; if `PressureForm`, fit ``P(V)``; if `BulkModulusForm`, fit ``B(V)``.
- `eos::EquationOfState`: a trial equation of state. If it has units, `xdata` and `ydata` must also have.
- `xdata::AbstractVector`: a vector of volumes (``V``), with(out) units.
- `ydata::AbstractVector`: a vector of energies (``E``), pressures (``P``), or bulk moduli (``B``), with(out) units. It must be consistent with `form`.
- `debug::Bool=false`: if `true`, then an `LsqFit.LsqFitResult` is returned, containing estimated Jacobian, residuals, etc.; if `false`, a fitted `EquationOfState` is returned. The default value is `false`.
- `kwargs`: the rest keyword arguments are the same as that of `LsqFit.curve_fit`. See its [documentation](https://github.com/JuliaNLSolvers/LsqFit.jl/blob/master/README.md)
    and [tutorial](https://julianlsolvers.github.io/LsqFit.jl/latest/tutorial/).
"""
function lsqfit(
    form::EquationForm,
    eos::EquationOfState{<:Real},
    xdata::AbstractVector{<:Real},
    ydata::AbstractVector{<:Real};
    debug = false,
    kwargs...,
)
    T = promote_type(eltype(xdata), eltype(ydata), Float64)
    E = constructorof(typeof(eos))
    model = (x, p) -> map(apply(form, E(p...)), x)
    fitted = curve_fit(
        model,
        T.(xdata),
        T.(ydata),
        T.(Collections.fieldvalues(eos)),
        kwargs...,
    )
    return debug ? fitted : E(fitted.param...)
end  # function lsqfit
function lsqfit(
    form::EquationForm,
    eos::EquationOfState{<:AbstractQuantity},
    xdata::AbstractVector{<:AbstractQuantity},
    ydata::AbstractVector{<:AbstractQuantity};
    kwargs...,
)
    E = constructorof(typeof(eos))
    values = Collections.fieldvalues(eos)
    original_units = map(unit, values)
    trial_params, xdata, ydata = [map(ustrip ∘ upreferred, x) for x in (
        values,
        xdata,
        ydata,
    )]
    result = lsqfit(form, E(trial_params...), xdata, ydata; kwargs...)
    if result isa EquationOfState
        data = Collections.fieldvalues(result)
        return E([data[i] * upreferred(u) |> u for (i, u) in enumerate(original_units)]...)
    end
    return result
end  # function lsqfit

end
