module EquationsOfState

using Reexport

include("Collections.jl")
include("NonlinearFitting.jl")
include("FiniteStrains.jl")
include("LinearFitting.jl")
include("volume.jl")

@reexport using .Collections

end # module
