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

getstrain(::EulerianStrain, v0::Real, v::Real) = ((v0 / v)^(2 / 3) - 1) / 2
getstrain(::LagrangianStrain, v0::Real, v::Real) = ((v / v0)^(2 / 3) - 1) / 2
getstrain(::NaturalStrain, v0::Real, v::Real) = log(v / v0) / 3
getstrain(::InfinitesimalStrain, v0::Real, v::Real) = 1 - (v0 / v)^(1 / 3)

end
