"""
# module NonlinearFitting



# Examples

```jldoctest
julia>
```
"""
module NonlinearFitting

using LsqFit: curve_fit
using Unitful
import Unitful: AbstractQuantity, 𝐋, 𝐌, 𝐓

import ..EquationForm
using ..Collections

export lsqfit

# This idea is borrowed from [SimpleTraits.jl](https://github.com/mauro3/SimpleTraits.jl/blob/master/src/SimpleTraits.jl).
abstract type Trait end
abstract type Not{T<:Trait} <: Trait end
struct HasUnit <: Trait end

_unit_trait(T::Type{<:Real}) = Not{HasUnit}
_unit_trait(T::Type{<:AbstractQuantity}) = HasUnit

"""
    lsqfit(form, eos, xdata, ydata; debug = false, kwargs...)

Fit an equation of state using least-squares fitting method (with the Levenberg-Marquardt algorithm).

# Arguments
- `form::EquationForm`: an `EquationForm` instance. If `EnergyForm`, fit ``E(V)``; if `PressureForm`, fit ``P(V)``; if `BulkModulusForm`, fit ``B(V)``.
- `eos::EquationOfState`: a trial equation of state.
- `xdata::AbstractVector`: a vector of volumes.
- `ydata::AbstractVector`: a vector of energies, pressures, or bulk moduli.
- `debug::Bool=false`: if `true`, then an `LsqFit.LsqFitResult` is returned, containing estimated Jacobian, residuals, etc.; if `false`, a fitted `EquationOfState` is returned. The default value is `false`.
- `kwargs`: the rest keyword arguments that will be sent to `LsqFit.curve_fit`. See its [documentation](https://github.com/JuliaNLSolvers/LsqFit.jl/blob/master/README.md).
"""
function lsqfit(
    form::EquationForm,
    eos::EquationOfState,
    xdata::AbstractVector,
    ydata::AbstractVector;
    kwargs...,
)
    T = eltype(eos)
    return lsqfit(_unit_trait(T), form, eos, xdata, ydata, kwargs...)
end # function lsqfit
function lsqfit(
    ::Type{Not{HasUnit}},
    form::EquationForm,
    eos::EquationOfState,
    xdata::AbstractVector,
    ydata::AbstractVector;
    debug = false,
    kwargs...,
)
    T = promote_type(eltype(eos), eltype(xdata), eltype(ydata), Float64)
    E = typeof(eos).name.wrapper
    model(x, p) = map(apply(form, E(p...)), x)
    fitted = curve_fit(model, T.(xdata), T.(ydata), T.(Collections.fieldvalues(eos)), kwargs...)
    return debug ? fitted : E(fitted.param...)
end  # function lsqfit
function lsqfit(
    ::Type{HasUnit},
    form::EquationForm,
    eos::EquationOfState,
    xdata::AbstractVector,
    ydata::AbstractVector;
    kwargs...,
)
    E = typeof(eos).name.wrapper
    units = unit.(Collections.fieldvalues(eos))
    trial_params = map(ustrip, Collections.fieldvalues(upreferred(eos)))
    xdata, ydata = map(ustrip ∘ upreferred, xdata), map(ustrip ∘ upreferred, ydata)
    result = lsqfit(form, E(trial_params...), xdata, ydata, kwargs...)
    if result isa EquationOfState
        E(
            [uconvert(u, Collections.fieldvalues(result)[i] * upreferred(u)) for (i, u) in enumerate(units)]...
        )
    end
end  # function lsqfit

end
