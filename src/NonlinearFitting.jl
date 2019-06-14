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

function lsqfit(f::Function, eos::E, xdata::Vector{T}, ydata::Vector{T}; debug::Bool = false, kwargs...) where {T <: AbstractFloat,E <: EquationOfState{T}}
    model(x, p) = f(E(p), x)
    fitted = curve_fit(model, xdata, ydata, collect(eos); kwargs...)
    debug ? fitted : E(fitted.param)
end  # function lsqfit
function lsqfit(f::Function, eos::E, xdata::X, ydata::Y; kwargs...) where {E <: EquationOfState,X <: AbstractVector,Y <: AbstractVector}
    T = promote_type(eltype(eos), eltype(xdata), eltype(ydata), Float64)
    lsqfit(f, convert(similar_type(E, T), eos), convert(Vector{T}, xdata), convert(Vector{T}, ydata); kwargs...)
end  # function lsqfit

fit_energy(eos::EquationOfState, xdata::AbstractVector, ydata::AbstractVector; kwargs...) = lsqfit(eval_energy, eos, xdata, ydata; kwargs...)

function fit_pressure(eos::EquationOfState, xdata::AbstractVector, ydata::AbstractVector; silent::Bool = true, kwargs...)
    silent || @info "Fitting pressure... The parameter `e0 = $(eos.e0)` is not used and will be kept as input."
    lsqfit(eval_pressure, eos, xdata, ydata; kwargs...)
end

function fit_bulk_modulus(eos::EquationOfState, xdata::AbstractVector, ydata::AbstractVector; silent::Bool = true, kwargs...)
    silent || @info "Fitting bulk modulus... The parameter `e0 = $(eos.e0)` is not used and will be kept as input."
    lsqfit(eval_bulk_modulus, eos, xdata, ydata; kwargs...)
end

end
