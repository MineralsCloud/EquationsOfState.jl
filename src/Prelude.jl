"""
# module Prelude



# Examples

```jldoctest
julia>
```
"""
module Prelude

export Maybe,
    MaybeData

# Referenced from [@JeffreySarnoff's answer](https://discourse.julialang.org/t/aliases-for-union-t-nothing-and-union-t-missing/15402/4?u=singularitti).
const Maybe{T} = Union{T,Nothing}
const MaybeData{T} = Union{T,Missing}

end
