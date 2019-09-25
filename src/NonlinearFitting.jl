"""
# module NonlinearFitting



# Examples

```jldoctest
julia>
```
"""
module NonlinearFitting

using LsqFit: curve_fit

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
    T = promote_type(eltype(eos), eltype(xdata), eltype(ydata), Float64)
    P = Collections.similar_type(E, T)
    model(x, p) = map(apply(form, P(p...)), x)
    fitted = curve_fit(model, T.(xdata), T.(ydata), T.(collect(eos)); kwargs...)
    return debug ? fitted : P(fitted.param...)
end  # function lsqfit

end
