"""
# module NonlinearFitting



# Examples

```jldoctest
julia>
```
"""
module NonlinearFitting

using LsqFit
using Transducers

using EquationsOfState.Collections

export fit_energy,
    fit_pressure,
    get_fitting_parameters

function get_fitting_parameters(eos::EquationOfState)
    collect(NotA(NonFittingParameter), get_parameters(eos))
end

function fit_energy(eos::EquationOfState, xdata::T, ydata::T; kwargs...) where {T <: AbstractVector}
    model(x, p) = map((eval_energy ∘ typeof(eos) ∘ collect)(DropLast(1), p), (x, last(p)))
    curve_fit(model, xdata, ydata, push!(get_fitting_parameters(eos), minimum(ydata)); kwargs...)
end

create_model(eos::EquationOfState) = (x::AbstractVector, p::AbstractVector) -> map((eval_pressure ∘ typeof(eos))(p), x)
create_model(eos::Holzapfel) = (x::AbstractVector, p::AbstractVector) -> map((eval_pressure ∘ typeof(eos))(push!(p, eos.z)), x)

function fit_pressure(eos::EquationOfState, xdata::T, ydata::T; kwargs...) where {T <: AbstractVector}
    model = create_model(eos)
    curve_fit(model, xdata, ydata, get_fitting_parameters(eos); kwargs...)
end

end