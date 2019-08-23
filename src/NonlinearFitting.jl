"""
# module NonlinearFitting



# Examples

```jldoctest
julia>
```
"""
module NonlinearFitting

using LsqFit: curve_fit

using EquationsOfState
using EquationsOfState.Collections

export lsqfit

function lsqfit(
    A::Type{<:EquationOfStateRelation},
    eos::E,
    xdata::Vector{T},
    ydata::Vector{T};
    debug::Bool = false, kwargs...
) where {T<:AbstractFloat,E<:EquationOfState{T}}
    model(x, p) = map(calculate(A, E(p)), x)
    fitted = curve_fit(model, xdata, ydata, collect(eos); kwargs...)
    debug ? fitted : E(fitted.param)
end  # function lsqfit
"""
    lsqfit(T, eos, xdata, ydata; debug = false, kwargs...)

Fit an equation of state using least-squares fitting method (with the Levenberg-Marquardt algorithm).

# Arguments
- `T::Type{<:EquationOfStateRelation}`: an `EquationOfStateRelation`. If it is `EnergyRelation`, fit ``E(V)``; if `PressureRelation`, fit ``P(V)``; if `BulkModulusRelation`, fit ``B(V)``.
- `eos::EquationOfState`: a trial equation of state.
- `xdata::AbstractVector`: a vector of volumes.
- `ydata::AbstractVector`: a vector of energies, pressures, or bulk moduli.
- `debug::Bool=false`: if `true`, then an `LsqFit.LsqFitResult` is returned, containing estimated Jacobian, residuals, etc.; if `false`, a fitted `EquationOfState` is returned. The default value is `false`.
- `kwargs`: the rest keyword arguments that will be sent to `LsqFit.curve_fit`. See its [documentation](https://github.com/JuliaNLSolvers/LsqFit.jl/blob/master/README.md).
"""
function lsqfit(
    A::Type{<:EquationOfStateRelation},
    eos::E,
    xdata::X,
    ydata::Y;
    kwargs...
) where {E<:EquationOfState,X<:AbstractVector,Y<:AbstractVector}
    T = promote_type(eltype(eos), eltype(xdata), eltype(ydata), Float64)
    lsqfit(A, convert(similar_type(E, T), eos), convert(Vector{T}, xdata), convert(Vector{T}, ydata); kwargs...)
end  # function lsqfit

end
