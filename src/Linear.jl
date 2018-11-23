"""
# module Linear



# Examples

```jldoctest
julia>
```
"""
module Linear

function finite_strain_at_order(m::Int, reference_volume::Float64, v::Float64)
    m == 1 && return -(1 / 3 / v) * (reference_volume / v)
    -(3 * m + 2) / (3 * v) * finite_strain_at_order(m - 1, reference_volume, v)
end

end