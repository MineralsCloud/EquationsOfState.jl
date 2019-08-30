module EquationsOfState

export EquationOfStateForm, EnergyForm, PressureForm, BulkModulusRelation

abstract type EquationOfStateForm end

struct EnergyForm <: EquationOfStateForm end
struct PressureForm <: EquationOfStateForm end
struct BulkModulusRelation <: EquationOfStateForm end

include("Collections.jl")
include("NonlinearFitting.jl")
include("FiniteStrains.jl")
include("LinearFitting.jl")
include("FindVolume.jl")

end # module
