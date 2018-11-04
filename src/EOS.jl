module EOS

using Reexport

include("Collections.jl")
@reexport using .Collections

include("Fitting.jl")
@reexport using .Fitting

end # module
