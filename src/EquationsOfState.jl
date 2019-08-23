module EquationsOfState

export EquationOfStateRelation, EnergyRelation, PressureRelation, BulkModulusRelation

abstract type EquationOfStateRelation end

struct EnergyRelation <: EquationOfStateRelation end
struct PressureRelation <: EquationOfStateRelation end
struct BulkModulusRelation <: EquationOfStateRelation end

include("Collections.jl")
include("NonlinearFitting.jl")
include("FiniteStrains.jl")
include("LinearFitting.jl")
include("FindVolume.jl")

end # module
