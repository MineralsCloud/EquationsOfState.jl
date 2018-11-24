"""
# module Linear



# Examples

```jldoctest
julia>
```
"""
module Linear

export finite_strain_at

function finite_strain_at(m::Int, v0::Float64, v::Float64)
    m == 1 && return -1 / 3 / v * (v0 / v)^(2 / 3)
    -(3 * m + 2) / (3 * v) * finite_strain_at(m - 1, v0, v)
end

end