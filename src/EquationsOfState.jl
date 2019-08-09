module EquationsOfState

using Reexport

include("Targets.jl")
include("Collections.jl")
include("NonlinearFitting.jl")
include("FiniteStrains.jl")
include("LinearFitting.jl")
include("NumericallyFindVolume.jl")

@reexport using .Targets

end # module
