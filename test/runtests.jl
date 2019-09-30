using EquationsOfState
using Test

@testset "EquationsOfState.jl" begin
    include("FiniteStrains.jl")
    include("NonlinearFitting.jl")
    include("LinearFitting.jl")
    # include("Find.jl")
end
