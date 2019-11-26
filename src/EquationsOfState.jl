module EquationsOfState

export EquationForm, EnergyForm, PressureForm, BulkModulusForm

abstract type EquationForm end
struct EnergyForm <: EquationForm end
struct PressureForm <: EquationForm end
struct BulkModulusForm <: EquationForm end

include("Collections.jl")
include("NonlinearFitting.jl")
include("FiniteStrains.jl")
include("LinearFitting.jl")
include("Find.jl")

end # module
