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
- `xdata::AbstractArray`: an array of volumes (``V``), with(out) units.
- `ydata::AbstractArray`: an array of energies (``E``), pressures (``P``), or bulk moduli (``B``), with(out) units. It must be consistent with `prop`.
- `debug::Bool=false`: if `true`, then an `LsqFit.LsqFitResult` is returned, containing estimated Jacobian, residuals, etc.; if `false`, a fitted `EquationOfState` is returned. The default value is `false`.
- `kwargs`: the rest keyword arguments are the same as that of `LsqFit.curve_fit`. See its [documentation](https://github.com/JuliaNLSolvers/LsqFit.jl/blob/master/README.md)
    and [tutorial](https://julianlsolvers.github.io/LsqFit.jl/latest/tutorial/).
"""
function lsqfit(f, xdata, ydata; kwargs...)
    eos, prop = fieldvalues(f)
    T = constructorof(typeof(eos))  # Get the `UnionAll` type
    params, xdata, ydata = _preprocess(eos, xdata, ydata)
    model = (x, p) -> map(T(p...)(prop), x)
    fit = curve_fit(model, xdata, ydata, params; kwargs...)
    return _postprocess(T(fit.param...), eos)
end # function lsqfit

struct _Data{S,T}
    data::T
end
_Data(data::T) where {T} = _Data{eltype(data),T}(data)

_preprocess(eos, xdata, ydata) = _preprocess(_Data(eos), _Data(xdata), _Data(ydata))  # Holy trait
_preprocess(eos::_Data{<:Real}, xdata::_Data{<:Real}, ydata::_Data{<:Real}) =
    float.(fieldvalues(eos.data)), float(xdata.data), float(ydata.data)
function _preprocess(
    eos::_Data{<:AbstractQuantity},
    xdata::_Data{<:AbstractQuantity},
    ydata::_Data{<:AbstractQuantity},
)
    values = fieldvalues(eos.data)
    original_units = unit.(values)  # Keep a record of `eos`'s units
    f = x -> map(float ∘ ustrip ∘ upreferred, x)  # Convert to preferred units and strip the unit
    return map(f, (values, xdata.data, ydata.data))
end # function _preprocess

_postprocess(eos, trial_eos) = _postprocess(eos, _Data(trial_eos))  # Holy trait
_postprocess(eos, trial_eos::_Data{<:Real}) = eos
function _postprocess(eos, trial_eos::_Data{<:AbstractQuantity})
    T = constructorof(typeof(trial_eos.data))  # Get the `UnionAll` type
    original_units = unit.(fieldvalues(trial_eos.data))  # Keep a record of `eos`'s units
    return T((
        x * upreferred(u) |> u for (x, u) in zip(fieldvalues(eos), original_units)
    )...)  # Convert back to original `eos`'s units
end # function _postprocess

end
