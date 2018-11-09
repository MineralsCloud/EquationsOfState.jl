"""
# module Fitting



# Examples

```jldoctest
julia>
```
"""
module Fitting

using LsqFit

using EOS.Collections

export fit_energy,
    fit_pressure,
    collect_fitting_parameters

function collect_fitting_parameters(eos::EquationOfState)
    filter(x -> !isa(x, NonFittingParameter), collect_parameters(eos))
end

function fit_energy(eos::T, xdata::Vector{Float64}, ydata::Vector{Float64}; kwargs...) where {T <: EquationOfState}
    model(x, p) = eval_energy(T(p[1:(end - 1)])).(x, p[end])
    curve_fit(model, xdata, ydata, [collect_fitting_parameters(eos); minimum(ydata)]; kwargs...)
end

function fit_pressure(eos::T, xdata::Vector{Float64}, ydata::Vector{Float64}; kwargs...) where {T <: EquationOfState}
    model(x, p) = eval_pressure(T(p)).(x)
    curve_fit(model, xdata, ydata, collect_fitting_parameters(eos); kwargs...)
end

function fit_pressure(eos::Holzapfel, xdata::Vector{Float64}, ydata::Vector{Float64}; kwargs...) where {T <: EquationOfState}
    model(x, p) = eval_pressure(T(p, eos.z)).(x)
    curve_fit(model, xdata, ydata, collect_fitting_parameters(eos); kwargs...)
end

end