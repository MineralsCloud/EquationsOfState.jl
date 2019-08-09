"""
# module NonlinearFitting



# Examples

```jldoctest
julia>
```
"""
module NonlinearFitting

using LsqFit: curve_fit

using EquationsOfState.Targets
using EquationsOfState.Collections

export lsqfit

function lsqfit(
    A::Type{<:EquationOfStateTarget},
    eos::E,
    xdata::Vector{T},
    ydata::Vector{T};
    debug::Bool = false, kwargs...
) where {T<:AbstractFloat,E<:EquationOfState{T}}
    model(x, p) = map(calculate(A, E(p)), x)
    fitted = curve_fit(model, xdata, ydata, collect(eos); kwargs...)
    debug ? fitted : E(fitted.param)
end  # function lsqfit
function lsqfit(
    A::Type{<:EquationOfStateTarget},
    eos::E,
    xdata::X,
    ydata::Y;
    kwargs...
) where {E<:EquationOfState,X<:AbstractVector,Y<:AbstractVector}
    T = promote_type(eltype(eos), eltype(xdata), eltype(ydata), Float64)
    lsqfit(A, convert(similar_type(E, T), eos), convert(Vector{T}, xdata), convert(Vector{T}, ydata); kwargs...)
end  # function lsqfit

end
