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
    lsqfit(eos(prop), volumes, ydata; kwargs...)

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
function lsqfit(f, volumes, ydata; kwargs...)
    eos, property = fieldvalues(f)
    T = constructorof(typeof(eos))  # Get the `UnionAll` type
    params, volumes, ydata = _preprocess(eos, volumes, ydata)
    model = (x, p) -> map(T(p...)(property), x)
    fit = curve_fit(model, volumes, ydata, params; kwargs...)
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
