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

function fit_energy(eos::EquationOfState{T}, xdata::AbstractVector{T}, ydata::AbstractVector{T}; kwargs...) where {T <: AbstractFloat}
    @. model(x, p) = eval_energy(typeof(eos)(p), x)
    curve_fit(model, xdata, ydata, collect(eos); kwargs...)
end
function fit_energy(eos::EquationOfState, xdata::AbstractVector, ydata::AbstractVector; kwargs...)
    T = promote_type(eltype(eos), eltype(xdata), eltype(ydata), Float64)
    fit_energy(convert_eltype(T, eos), convert_eltype(T, xdata), convert_eltype(T, ydata); kwargs...)
end

function fit_pressure(eos::EquationOfState{T}, xdata::AbstractVector{T}, ydata::AbstractVector{T}; kwargs...) where {T <: AbstractFloat}
    @. model(x, p) = eval_pressure(typeof(eos)(p), x)
    curve_fit(model, xdata, ydata, collect(eos); kwargs...)
end
function fit_pressure(eos::EquationOfState, xdata::AbstractVector, ydata::AbstractVector; kwargs...)
    T = promote_type(eltype(eos), eltype(xdata), eltype(ydata), Float64)
    fit_pressure(convert_eltype(T, eos), convert_eltype(T, xdata), convert_eltype(T, ydata); kwargs...)
end

function fit_bulk_modulus(eos::EquationOfState{T}, xdata::AbstractVector{T}, ydata::AbstractVector{T}; kwargs...) where {T <: AbstractFloat}
    @. model(x, p) = eval_bulk_modulus(typeof(eos)(p), x)
    curve_fit(model, xdata, ydata, collect(eos); kwargs...)
end
function fit_bulk_modulus(eos::EquationOfState, xdata::AbstractVector, ydata::AbstractVector; kwargs...)
    T = promote_type(eltype(eos), eltype(xdata), eltype(ydata), Float64)
    fit_bulk_modulus(convert_eltype(T, eos), convert_eltype(T, xdata), convert_eltype(T, ydata); kwargs...)
end

end