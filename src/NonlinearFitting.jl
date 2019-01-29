"""
# module NonlinearFitting



# Examples

```jldoctest
julia>
```
"""
module NonlinearFitting

using LsqFit

using EquationsOfState.Collections

export fit_energy,
    fit_pressure,
    get_fitting_parameters

function get_fitting_parameters(eos::EquationOfState)
    Iterators.filter(x -> !(x isa NonFittingParameter), get_parameters(eos)) |> Tuple
end

function fit_energy(eos::T, xdata::S, ydata::S; kwargs...) where {T <: EquationOfState, S <: AbstractVector}
    model(x, p) = eval_energy(T(p[1:(end - 1)])).(x, p[end])
    curve_fit(model, xdata, ydata, [get_fitting_parameters(eos); minimum(ydata)]; kwargs...)
end

create_model(eos::EquationOfState) = (x::AbstractVector, p::AbstractVector) -> eval_pressure(T(p)).(x)
create_model(eos::Holzapfel) = (x::AbstractVector, p::AbstractVector) -> eval_pressure(T([p; eos.z])).(x)

function fit_pressure(eos::EquationOfState, xdata::S, ydata::S; kwargs...) where {S <: AbstractVector}
    model = create_model(eos)
    curve_fit(model, xdata, ydata, get_fitting_parameters(eos); kwargs...)
end

end