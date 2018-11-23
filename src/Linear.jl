"""
# module Linear



# Examples

```jldoctest
julia>
```
"""
module Linear

function finite_strain_at_order(m::Val{N}, reference_volume::Float64, v::Float64) where {N}
    -(3 * m + 2) / (3 * v) * finite_strain_at_order(Val(m - 1), reference_volume, v)
end
finite_strain_at_order(::Val{1}, reference_volume::Float64, v::Float64) = -(1 / 3 / v) * (reference_volume / v)
finite_strain_at_order(m::Int, reference_volume::Float64, v::Float64) = finite_strain_at_order(Val(m), reference_volume, v)

end