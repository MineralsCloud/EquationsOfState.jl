"""
# module NonlinearFitting



# Examples

```jldoctest
julia>
```
"""
module NonlinearFitting

using LsqFit: curve_fit
using StaticArrays: similar_type

import ..EquationOfStateForm
using ..Collections

export lsqfit

function lsqfit(
    form::EquationOfStateForm,
    eos::E,
    xdata::Vector{T},
    ydata::Vector{T};
    debug::Bool = false, kwargs...
) where {T<:AbstractFloat,E<:EquationOfState{T}}
    model(x, p) = map(apply(form, E(p)), x)
    fitted = curve_fit(model, xdata, ydata, collect(eos); kwargs...)
    debug ? fitted : E(fitted.param)
end  # function lsqfit
"""
    lsqfit(form, eos, xdata, ydata; debug = false, kwargs...)

Fit an equation of state using least-squares fitting method (with the Levenberg-Marquardt algorithm).

# Arguments
- `form::EquationOfStateForm`: an `EquationOfStateForm` instance. If `EnergyForm`, fit ``E(V)``; if `PressureForm`, fit ``P(V)``; if `BulkModulusForm`, fit ``B(V)``.
- `eos::EquationOfState`: a trial equation of state.
- `xdata::AbstractVector`: a vector of volumes.
- `ydata::AbstractVector`: a vector of energies, pressures, or bulk moduli.
- `debug::Bool=false`: if `true`, then an `LsqFit.LsqFitResult` is returned, containing estimated Jacobian, residuals, etc.; if `false`, a fitted `EquationOfState` is returned. The default value is `false`.
- `kwargs`: the rest keyword arguments that will be sent to `LsqFit.curve_fit`. See its [documentation](https://github.com/JuliaNLSolvers/LsqFit.jl/blob/master/README.md).
"""
function lsqfit(
    form::EquationOfStateForm,
    eos::E,
    xdata::X,
    ydata::Y;
    kwargs...
) where {E<:EquationOfState,X<:AbstractVector,Y<:AbstractVector}
    T = promote_type(eltype(eos), eltype(xdata), eltype(ydata), Float64)
    lsqfit(form, convert(similar_type(E, T), eos), convert(Vector{T}, xdata), convert(Vector{T}, ydata); kwargs...)
end  # function lsqfit

end
