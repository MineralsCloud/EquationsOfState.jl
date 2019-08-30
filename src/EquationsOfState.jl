module EquationsOfState

export EquationOfStateForm, EnergyRelation, PressureRelation, BulkModulusRelation

abstract type EquationOfStateForm end

struct EnergyRelation <: EquationOfStateForm end
struct PressureRelation <: EquationOfStateForm end
struct BulkModulusRelation <: EquationOfStateForm end

include("Collections.jl")
include("NonlinearFitting.jl")
include("FiniteStrains.jl")
include("LinearFitting.jl")
include("FindVolume.jl")

end # module
