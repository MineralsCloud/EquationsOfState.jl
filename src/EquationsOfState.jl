module EquationsOfState

using Reexport

include("Targets.jl")
include("Collections.jl")
include("NonlinearFitting.jl")
include("FiniteStrains.jl")
include("LinearFitting.jl")
include("volume.jl")

@reexport using .Collections

end # module
