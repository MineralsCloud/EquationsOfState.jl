#=
volume:
- Julia version: 1.0
- Author: singularitti
- Date: 2019-08-06
=#
using IntervalRootFinding

using EquationsOfState.Collections

export find_volume

function find_volume(eos, energies, interval, method)
    f(v) = eval_energy(eos, v)
    solutions = roots(f, interval, method)
    length(solutions) != 1 ? error("Multiple roots find!") : return first(solutions)
end # function find_volume
