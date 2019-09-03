"""
# module FiniteStrains



# Examples

```jldoctest
julia>
```
"""
module FiniteStrains

export FiniteStrain, EulerianStrain, LagrangianStrain, NaturalStrain, InfinitesimalStrain, getstrain

abstract type FiniteStrain end

struct EulerianStrain <: FiniteStrain end
struct LagrangianStrain <: FiniteStrain end
struct NaturalStrain <: FiniteStrain end
struct InfinitesimalStrain <: FiniteStrain end

getstrain(::EulerianStrain, v0::Real, v::Real) = (cbrt(v0 / v)^2 - 1) / 2
getstrain(::LagrangianStrain, v0::Real, v::Real) = (cbrt(v / v0)^2 - 1) / 2
getstrain(::NaturalStrain, v0::Real, v::Real) = log(v / v0) / 3
getstrain(::InfinitesimalStrain, v0::Real, v::Real) = 1 - cbrt(v0 / v)

end
