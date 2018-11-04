module EoS

using Reexport

include("CollectionsOfEoS.jl")
@reexport using .CollectionsOfEoS

include("Fitting.jl")
@reexport using .Fitting

end # module
