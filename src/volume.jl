#=
volume:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-08-06
=#
using IntervalRootFinding

using EquationsOfState.Targets
using EquationsOfState.Collections

export find_volume

function find_volume(T::Type{<: EquationOfStateTarget}, eos, y, interval, method)
    f = v -> calculate(T, eos, v)
    solutions = roots(f, interval, method)
    length(solutions) != 1 ? error("Multiple roots find!") : return first(solutions)
end # function find_volume
