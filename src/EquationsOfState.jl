module EquationsOfState

export EquationOfStateForm, EnergyForm, PressureRelation, BulkModulusRelation

abstract type EquationOfStateForm end

struct EnergyForm <: EquationOfStateForm end
struct PressureRelation <: EquationOfStateForm end
struct BulkModulusRelation <: EquationOfStateForm end

include("Collections.jl")
include("NonlinearFitting.jl")
include("FiniteStrains.jl")
include("LinearFitting.jl")
include("FindVolume.jl")

end # module
