"""
# module FindVolume



# Examples

```jldoctest
julia>
```
"""
module FindVolume

using Roots: find_zero

using EquationsOfState
using EquationsOfState.Collections

export findvolume

function findvolume(form::EquationOfStateForm, eos::EquationOfState, y::Real, interval, method)
    f(v) = apply(form, eos, v) - y
    return find_zero(f, (minimum(interval), maximum(interval)), method)
end # function findvolume

end
