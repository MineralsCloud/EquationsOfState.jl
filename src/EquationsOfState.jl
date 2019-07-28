module EquationsOfState

using Reexport

include("Collections.jl")
include("NonlinearFitting.jl")
include("FiniteStrains.jl")
include("LinearFitting.jl")

@reexport using .Collections

end # module
