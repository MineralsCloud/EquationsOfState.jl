module EOS

using StaticArrays
using CurveFit
using LsqFit: @., curve_fit

export fit_energy,
    fit_pressure,
    eval_energy,
    eval_pressure
    EquationOfState,
    Birch,
    Murnaghan,
    BirchMurnaghan2nd, BirchMurnaghan3rd, BirchMurnaghan4th,
    Vinet

abstract type EquationOfState end

struct Birch <: EquationOfState
    parameters::SVector{3, Float64}
end

struct Murnaghan <: EquationOfState
    parameters::SVector{3, Float64}
end

struct BirchMurnaghan2nd <: EquationOfState
    parameters::SVector{2, Float64}
end

struct BirchMurnaghan3rd <: EquationOfState
    parameters::SVector{3, Float64}
end

struct BirchMurnaghan4th <: EquationOfState
    parameters::SVector{4, Float64}
end

struct Vinet <: EquationOfState
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

function eval_energy(eos::Murnaghan)::Function
    v0, b0, bp0 = eos.parameters

    function (v::Float64, f0::Float64=0)::Float64
        x = bp0 - 1
        y = (v0 / v)^bp0
        return f0 + b0 / bp0 * v * (y / x + 1) - v0 * b0 / x
    end
end

function eval_pressure(eos::Murnaghan)::Function
    v0, b0, bp0 = eos.parameters

    function (v::Float64)::Float64
        return b0 / bp0 * ((v0 / v)^bp0 - 1)
    end
end

function eval_energy(eos::BirchMurnaghan2nd)::Function
    v0, b0 = eos.parameters

    function (v::Float64, f0::Float64=0)::Float64
        f = ((v0 / v)^(2 / 3) - 1) / 2
        return f0 + 9 / 2 * b0 * v0 * f^2
    end
end

function eval_pressure(eos::BirchMurnaghan2nd)::Function
    v0, b0 = eos.parameters

    function (v::Float64)::Float64
        f = ((v0 / v)^(2 / 3) - 1) / 2
        return 3 * b0 * f * (1 + 2 * f)^(5 / 2)
    end
end

function eval_energy(eos::BirchMurnaghan3rd)::Function
    v0, b0, bp0 = eos.parameters

    function (v::Float64, f0::Float64=0)::Float64
        eta = (v0 / v)^(1 / 3)
        xi = eta^2 - 1
        return f0 + 9 / 16 * b0 * v0 * xi^2 * (6 + bp0 * xi - 4 * eta^2)
    end
end

function eval_pressure(eos::BirchMurnaghan3rd)::Function
    v0, b0, bp0 = eos.parameters

    function (v::Float64)::Float64
        eta = (v0 / v)^(1 / 3)
        return 3 / 2 * b0 * (eta^7 - eta^5) * (1 + 3 / 4 * (bp0 - 4) * (eta^2 - 1))
    end
end

function eval_energy(eos::BirchMurnaghan4th)::Function
    v0, b0, bp0, bpp0 = eos.parameters

    function (v::Float64, f0::Float64=0)::Float64
        f = ((v0 / v)^(2 / 3) - 1) / 2
        h = b0 * bpp0 + bp0^2
        return f0 + 3 / 8 * v0 * b0 * f^2 * ((9 * h - 63 * bp0 + 143) * f^2 + 12 * (bp0 - 4) * f + 12)
    end
end

function eval_pressure(eos::BirchMurnaghan4th)::Function
    v0, b0, bp0, bpp0 = eos.parameters

    function (v::Float64)::Float64
        f = ((v0 / v)^(2 / 3) - 1) / 2
        h = b0 * bpp0 + bp0^2
        return 1 / 2 * b0 * (2 * f + 1)^(5 / 2) * ((9 * h - 63 * bp0 + 143) * f^2 + 9 * (bp0 - 4) * f + 6)
    end
end

function eval_energy(eos::Vinet)::Function
    v0, b0, bp0 = eos.parameters

    function (v::Float64, f0::Float64=0)::Float64
        x = (v / v0)^(1 / 3)
        xi = 3 / 2 * (bp0 - 1)
        return f0 + 9 * b0 * v0 / xi^2 * (1 + (xi * (1 - x) - 1) * np.exp(xi * (1 - x)))
    end
end

function eval_pressure(eos::Vinet)::Function
    v0, b0, bp0 = eos.parameters

    function (v::Float64)::Float64
        x = (v / v0)^(1 / 3)
        xi = 3 / 2 * (bp0 - 1)
        return 3 * b0 / x^2 * (1 - x) * np.exp(xi * (1 - x))
    end
end


end # module
