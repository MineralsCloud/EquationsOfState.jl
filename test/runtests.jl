using EquationsOfState
using Test

@testset "EquationsOfState.jl" begin
    include("CollectionsTests.jl")
    include("NonlinearFittingTests.jl")
    include("LinearFittingTests.jl")
end
