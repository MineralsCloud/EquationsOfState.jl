"""
# module FindVolume



# Examples

```jldoctest
julia>
```
"""
module FindVolume

using Roots

using EquationsOfState
using EquationsOfState.Collections

export find_volume

function find_volume(form::EquationOfStateForm, eos::EquationOfState, y::Real, interval, method)
    f(v) = apply(form, eos, v) - y
    return find_zero(f, (minimum(interval), maximum(interval)), method)
end # function find_volume

end
