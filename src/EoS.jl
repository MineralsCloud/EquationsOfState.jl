module EOS

using StaticArrays
using CurveFit

export EquationOfState,
    Birch

abstract type EquationOfState end

struct Birch <: EquationOfState
    parameters::SVector{3, Float64}
end

function free_energy(eos::Birch, v::Float64, f0::Float64=0)::Float64
    v0, b0, bp0 = eos.parameters
    x = (v0 / v)^(2 / 3) - 1
    xi = 9 / 16 * b0 * v0 * x^2
    return f0 + 2 * xi + (bp0 - 4) * xi * x
end

function pressure(eos::Birch, v::Float64)::Float64
    v0, b0, bp0 = eos.parameters
    x = v0 / v
    xi = x^(2 / 3) - 1
    return 3 / 8 * b0 * x^(5 / 3) * xi * (4 + 3 * (bp0 - 4) * xi)
end

end # module
