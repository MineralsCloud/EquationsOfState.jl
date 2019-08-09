"""
# module Targets



# Examples

```jldoctest
julia>
```
"""
module Targets

export EquationOfStateTarget, EnergyTarget, PressureTarget, BulkModulusTarget

abstract type EquationOfStateTarget end

struct EnergyTarget <: EquationOfStateTarget end
struct PressureTarget <: EquationOfStateTarget end
struct BulkModulusTarget <: EquationOfStateTarget end

end
