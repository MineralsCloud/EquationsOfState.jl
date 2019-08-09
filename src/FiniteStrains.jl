"""
# module FiniteStrains



# Examples

```jldoctest
julia>
```
"""
module FiniteStrains

export FiniteStrain, EulerianStrain, LagrangianStrain, NaturalStrain, InfinitesimalStrain, get_strain

abstract type FiniteStrain{T} end

const EulerianStrain = FiniteStrain{:Eulerian}
const LagrangianStrain = FiniteStrain{:Lagrangian}
const NaturalStrain = FiniteStrain{:Natural}
const InfinitesimalStrain = FiniteStrain{:Infinitesimal}

get_strain(::Type{EulerianStrain}, v0::Real, v::Real) = ((v0 / v)^(2 / 3) - 1) / 2
get_strain(::Type{LagrangianStrain}, v0::Real, v::Real) = ((v / v0)^(2 / 3) - 1) / 2
get_strain(::Type{NaturalStrain}, v0::Real, v::Real) = log(v / v0) / 3
get_strain(::Type{InfinitesimalStrain}, v0::Real, v::Real) = 1 - (v0 / v)^(1 / 3)

end
