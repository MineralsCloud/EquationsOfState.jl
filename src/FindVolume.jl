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

function findvolume(form::EquationOfStateForm, eos::EquationOfState, y::Real, closepoint::Real, method)
    f(v) = apply(form, eos, v) - y
    return find_zero(f, closepoint, method)
end # function findvolume

end
