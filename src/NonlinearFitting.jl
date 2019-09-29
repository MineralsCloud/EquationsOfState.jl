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
    eos::E,
    xdata::AbstractVector,
    ydata::AbstractVector;
    debug = false,
    kwargs...
) where {E<:EquationOfState}

    # T = promote_type(eltype(eos), eltype(xdata), eltype(ydata))
    # T <: 
    P = Collections.similar_type(E, T)
    model(x, p) = map(apply(form, P(p...)), x)
    return lsqfit(model, xdata, ydata, debug = debug, kwargs...)
end  # function lsqfit
function lsqfit(
    form::EnergyForm,
    eos::EquationOfState,
    xdata::AbstractVector{A},
    ydata::AbstractVector{B};
    debug = false,
    kwargs...
) where {T<:Real,A<:AbstractQuantity{T,𝐋^3},B<:AbstractQuantity{T,𝐋^2*𝐌*𝐓^-2}}
    @assert(eltype(eos) <: AbstractQuantity, "The equation of state must have units!")
    xdata, ydata = uconvert.(u"angstrom^3", xdata), uconvert.(u"eV", ydata)
    E = typeof(eos).name.wrapper
    model(x, p) = map(apply(form, E(p...)), x)
    return lsqfit(model, xdata, ydata, debug = debug, kwargs...)
end  # function lsqfit
function lsqfit(
    model::Function,
    xdata::AbstractVector{A},
    ydata::AbstractVector{B},
    trial_params::AbstractVector{B};
    debug = false,
    kwargs...
) where {A<:AbstractFloat,B<:AbstractFloat}
    fitted = curve_fit(model, xdata, ydata, trial_params; kwargs...)
    return debug ? fitted : fitted.param
end  # function lsqfit

end
