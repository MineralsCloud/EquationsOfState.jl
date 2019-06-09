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
    fit_bulk_modulus

convert_eltype(T::Type, a) = map(x -> convert(T, x), a)

abstract type FitTrait end
struct FitEnergy <: FitTrait end
struct FitPressure <: FitTrait end
struct FitBulkModulus <: FitTrait end

eval_of(::FitTrait) = error("Unsupported!")
eval_of(::FitEnergy) = eval_energy
eval_of(::FitPressure) = eval_pressure
eval_of(::FitBulkModulus) = eval_bulk_modulus

function lsqfit(F::FitTrait, eos::EquationOfState{T}, xdata::AbstractVector{T}, ydata::AbstractVector{T}; kwargs...) where {T <: AbstractFloat}
    @. model(x, p) = eval_of(F)(typeof(eos)(p), x)
    curve_fit(model, xdata, ydata, collect(eos); kwargs...)
end  # function lsqfit
function lsqfit(F::FitTrait, eos::EquationOfState, xdata::AbstractVector, ydata::AbstractVector; kwargs...)
    T = promote_type(eltype(eos), eltype(xdata), eltype(ydata), Float64)
    lsqfit(F, convert_eltype(T, eos), convert_eltype(T, xdata), convert_eltype(T, ydata); kwargs...)
end  # function lsqfit

fit_energy(eos::EquationOfState, xdata::AbstractVector, ydata::AbstractVector; kwargs...) = lsqfit(FitEnergy, eos, xdata, ydata; kwargs...)
fit_pressure(eos::EquationOfState, xdata::AbstractVector, ydata::AbstractVector; kwargs...) = lsqfit(FitPressure, eos, xdata, ydata; kwargs...)
fit_bulk_modulus(eos::EquationOfState, xdata::AbstractVector, ydata::AbstractVector; kwargs...) = lsqfit(FitBulkModulus, eos, xdata, ydata; kwargs...)

end