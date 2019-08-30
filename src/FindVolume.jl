"""
# module FindVolume



# Examples

```jldoctest
julia>
```
"""
module FindVolume

using IntervalRootFinding

using EquationsOfState
using EquationsOfState.Collections

export find_volume

function find_volume(T::Type{<:EquationOfStateForm}, eos::EquationOfState, y::Real, interval, method)
    f = v -> apply(T, eos, v) - y
    solutions = roots(f, interval, method)
    length(solutions) != 1 ? error("Multiple roots find!") : return first(solutions)
end # function find_volume

end
