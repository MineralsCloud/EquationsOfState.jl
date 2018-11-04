module EOS

using Reexport

include("CollectionsOfEOS.jl")
@reexport using .CollectionsOfEOS

include("Fitting.jl")
@reexport using .Fitting

end # module
