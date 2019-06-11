"""
# module NonlinearFitting



# Examples

```jldoctest
julia>
```
"""
module NonlinearFitting

using LsqFit: curve_fit

using EquationsOfState.Collections

export fit_energy,
    fit_pressure,
    fit_bulk_modulus,
    FitTrait,
    FitEnergy,
    FitPressure,
    FitBulkModulus,
    eval_of

abstract type FitTrait end
struct FitEnergy <: FitTrait end
struct FitPressure <: FitTrait end
struct FitBulkModulus <: FitTrait end

eval_of(::Type{FitEnergy}) = eval_energy
eval_of(::Type{FitPressure}) = eval_pressure
eval_of(::Type{FitBulkModulus}) = eval_bulk_modulus

function lsqfit(F::Type{<: FitTrait}, eos::EquationOfState{T}, xdata::AbstractVector{T}, ydata::AbstractVector{T}; kwargs...) where {T <: AbstractFloat}
    X = typeof(eos)
    model(x, p) = eval_of(F)(X(p), x)
    fitted = curve_fit(model, xdata, ydata, collect(eos); kwargs...)
    X(fitted.param)
end  # function lsqfit
function lsqfit(F::Type{<: FitTrait}, eos::E, xdata::X, ydata::Y; kwargs...) where {E <: EquationOfState,X <: AbstractVector,Y <: AbstractVector}
    T = promote_type(eltype(eos), eltype(xdata), eltype(ydata), Float64)
    lsqfit(F, convert(E{T}, eos), convert(X{T}, xdata), convert(Y{T}, ydata); kwargs...)
end  # function lsqfit

fit_energy(eos::EquationOfState, xdata::AbstractVector, ydata::AbstractVector; kwargs...) = lsqfit(FitEnergy, eos, xdata, ydata; kwargs...)
fit_pressure(eos::EquationOfState, xdata::AbstractVector, ydata::AbstractVector; kwargs...) = lsqfit(FitPressure, eos, xdata, ydata; kwargs...)
fit_bulk_modulus(eos::EquationOfState, xdata::AbstractVector, ydata::AbstractVector; kwargs...) = lsqfit(FitBulkModulus, eos, xdata, ydata; kwargs...)

end
