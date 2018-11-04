"""
# module Fitting



# Examples

```jldoctest
julia>
```
"""
module Fitting

using LsqFit

using EoS.CollectionsOfEoS

export fit_energy,
    fit_pressure

function fit_energy(eos::T, xdata::Vector{Float64}, ydata::Vector{Float64}; kwargs...) where {T <: EquationOfState}
    model(x, p) = eval_energy(T(p[1:(end - 1)])).(x, p[end])
    curve_fit(model, xdata, ydata, [eos.parameters; minimum(ydata)]; kwargs...)
end

function fit_pressure(eos::T, xdata::Vector{Float64}, ydata::Vector{Float64}; kwargs...) where {T <: EquationOfState}
    model(x, p) = eval_pressure(T(p)).(x)
    curve_fit(model, xdata, ydata, eos.parameters; kwargs...)
end

end