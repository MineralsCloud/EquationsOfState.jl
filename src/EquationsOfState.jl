module EquationsOfState

export EquationOfStateForm, EnergyForm, PressureForm, BulkModulusForm

abstract type EquationOfStateForm end

struct EnergyForm <: EquationOfStateForm end
struct PressureForm <: EquationOfStateForm end
struct BulkModulusForm <: EquationOfStateForm end

include("Collections.jl")
include("NonlinearFitting.jl")
include("FiniteStrains.jl")
include("LinearFitting.jl")
include("FindVolume.jl")

end # module
