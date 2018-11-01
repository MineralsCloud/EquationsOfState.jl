module EOS

using StaticArrays
using CurveFit
using LsqFit: @., curve_fit

export fit_energy,
    fit_pressure,
    EquationOfState,
    Birch

abstract type EquationOfState end

struct Birch <: EquationOfState
    parameters::SVector{3, Float64}
end

function fit_energy(eos::T, xdata::Vector{Float64}, ydata::Vector{Float64}, initial_parameters::Vector{Float64}; kwargs...) where {T <: EquationOfState}
    @. model(x, p) = eval_energy(T(p))(x, p[end])
    curve_fit(model, xdata, ydata, initial_parameters; kwargs...)
end

function fit_pressure(eos::T, xdata::Vector{Float64}, ydata::Vector{Float64}, initial_parameters::Vector{Float64}; kwargs...) where {T <: EquationOfState}
    @. model(x, p) = eval_pressure(T(p))(x)
    curve_fit(model, xdata, ydata, initial_parameters; kwargs...)
end

function eval_energy(eos::Birch)::Function
    v0, b0, bp0 = eos.parameters

    function (v::Float64, f0::Float64=0)::Float64
        x = (v0 / v)^(2 / 3) - 1
        xi = 9 / 16 * b0 * v0 * x^2
        return f0 + 2 * xi + (bp0 - 4) * xi * x
    end
end

function eval_pressure(eos::Birch)::Function
    v0, b0, bp0 = eos.parameters

    function (v::Float64)::Float64
        x = v0 / v
        xi = x^(2 / 3) - 1
        return 3 / 8 * b0 * x^(5 / 3) * xi * (4 + 3 * (bp0 - 4) * xi)
    end
end

end # module
