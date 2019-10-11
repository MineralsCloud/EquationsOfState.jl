"""
# module NonlinearFitting



# Examples

```jldoctest
julia>
```
"""
module NonlinearFitting

using ConstructionBase: constructorof
using LsqFit: curve_fit
using Unitful: AbstractQuantity, upreferred, ustrip, unit

import ..EquationForm
using ..Collections

export lsqfit

"""
    lsqfit(form, eos, xdata, ydata; debug = false, kwargs...)

Fit an equation of state using least-squares fitting method (with the Levenberg-Marquardt algorithm).

# Arguments
- `form::EquationForm`: an `EquationForm` instance. If `EnergyForm`, fit ``E(V)``; if `PressureForm`, fit ``P(V)``; if `BulkModulusForm`, fit ``B(V)``.
- `eos::EquationOfState`: a trial equation of state.
- `xdata::AbstractVector`: a vector of volumes.
- `ydata::AbstractVector`: a vector of energies, pressures, or bulk moduli.
- `debug::Bool=false`: if `true`, then an `LsqFit.LsqFitResult` is returned, containing estimated Jacobian, residuals, etc.; if `false`, a fitted `EquationOfState` is returned. The default value is `false`.
- `kwargs`: the rest keyword arguments that will be sent to `LsqFit.curve_fit`. See its [documentation](https://github.com/JuliaNLSolvers/LsqFit.jl/blob/master/README.md).

# Examples
```jldoctest

    # Data in the following tests are from
    # https://github.com/materialsproject/pymatgen/blob/1f0957b8525ddc7d12ea348a19caecebe6c7ff34/pymatgen/analysis/tests/test_eos.py

# no unit
volumes = [
    159.9086,
    162.5738,
    165.2389,
    167.9041,
    170.5692,
    173.2344,
    175.8995,
    178.5647,
    181.2298,
    183.8949,
    186.5601,
    189.2252,
    191.8904,
    194.5555,
    197.2207,
    199.8858,
    202.551,
    205.2161,
    207.8812,
    210.5464,
    213.2115,
    215.8767,
    218.5418,
    221.207,
    223.8721,
    226.5373,
    229.2024,
    231.8675,
    234.5327,
    237.1978,
    239.863,
    242.5281,
    245.1933,
    247.8584,
    250.5236,
    253.1887,
    255.8538,
    258.519,
    261.1841,
    263.8493,
    266.5144,
    269.1796,
    271.8447,
    274.5098,
    277.175,
    279.8401,
    282.5053,
    285.1704,
    287.8356,
    290.5007,
    293.1659,
    295.831,
    298.4961,
    301.1613,
    303.8264,
    306.4916,
    309.1567,
    311.8219,
    314.487,
    317.1522,
    319.8173
]

energies = [
    -323.4078898,
    -323.4089153,
    -323.4098546,
    -323.410722,
    -323.4115195,
    -323.4122481,
    -323.4129189,
    -323.413528,
    -323.4140871,
    -323.4145889,
    -323.4150471,
    -323.415459,
    -323.4158302,
    -323.4161579,
    -323.4164498,
    -323.4167071,
    -323.4169305,
    -323.4171194,
    -323.4172809,
    -323.4174144,
    -323.4175216
    -323.4176029,
    -323.417661,
    -323.4176975,
    -323.41771,
    -323.4177051,
    -323.417682,
    -323.4176375,
    -323.417579,
    -323.4175048,
    -323.4174142,
    -323.4173101,
    -323.4171922,
    -323.4170611,
    -323.4169184,
    -323.4167647,
    -323.4166002,
    -323.4164244,
    -323.4162386,
    -323.4160446,
    -323.4158421,
    -323.4156312,
    -323.4154125,
    -323.4151861,
    -323.4149528,
    -323.4147131,
    -323.414467,
    -323.414215,
    -323.4139583,
    -323.4136953,
    -323.4134285,
    -323.4131559,
    -323.4128797,
    -323.4125984,
    -323.4123147,
    -323.4120269,
    -323.411736,
    -323.4114399,
    -323.4111421,
    -323.4108418,
    -323.4105393
]

julia>lsqfit(EnergyForm(),BirchMurnaghan3rd(224, 9, 3.75, -161), volumes, energies)
BirchMurnaghan3rd{Float64}(224.4445651526474, 0.0006250620492341443, 3.740368611188446, -323.4177142148064)

#have units
volumes_u = [
    159.9086 a₀^3
    162.5738 a₀^3
    165.2389 a₀^3
    167.9041 a₀^3
    170.5692 a₀^3
    173.2344 a₀^3
    175.8995 a₀^3
    178.5647 a₀^3
    181.2298 a₀^3
    183.8949 a₀^3
    186.5601 a₀^3
    189.2252 a₀^3
    191.8904 a₀^3
    194.5555 a₀^3
    197.2207 a₀^3
    199.8858 a₀^3
    202.551 a₀^3
    205.2161 a₀^3
    207.8812 a₀^3
    210.5464 a₀^3
    213.2115 a₀^3
    215.8767 a₀^3
    218.5418 a₀^3
    221.207 a₀^3
    223.8721 a₀^3
    226.5373 a₀^3
    229.2024 a₀^3
    231.8675 a₀^3
    234.5327 a₀^3
    237.1978 a₀^3
    239.863 a₀^3
    242.5281 a₀^3
    245.1933 a₀^3
    247.8584 a₀^3
    250.5236 a₀^3
    253.1887 a₀^3
    255.8538 a₀^3
    258.519 a₀^3
    261.1841 a₀^3
    263.8493 a₀^3
    266.5144 a₀^3
    269.1796 a₀^3
    271.8447 a₀^3
    274.5098 a₀^3
    277.175 a₀^3
    279.8401 a₀^3
    282.5053 a₀^3
    285.1704 a₀^3
    287.8356 a₀^3
    290.5007 a₀^3
    293.1659 a₀^3
    295.831 a₀^3
    298.4961 a₀^3
    301.1613 a₀^3
    303.8264 a₀^3
    306.4916 a₀^3
    309.1567 a₀^3
    311.8219 a₀^3
    314.487 a₀^3
    317.1522 a₀^3
    319.8173 a₀^3
]

energies_u = [
    -323.4078898 Ry
    -323.4089153 Ry
    -323.4098546 Ry
    -323.410722 Ry
    -323.4115195 Ry
    -323.4122481 Ry
    -323.4129189 Ry
    -323.413528 Ry
    -323.4140871 Ry
    -323.4145889 Ry
    -323.4150471 Ry
    -323.415459 Ry
    -323.4158302 Ry
    -323.4161579 Ry
    -323.4164498 Ry
    -323.4167071 Ry
    -323.4169305 Ry
    -323.4171194 Ry
    -323.4172809 Ry
    -323.4174144 Ry
    -323.4175216 Ry
    -323.4176029 Ry
    -323.417661 Ry
    -323.4176975 Ry
    -323.41771 Ry
    -323.4177051 Ry
    -323.417682 Ry
    -323.4176375 Ry
    -323.417579 Ry
    -323.4175048 Ry
    -323.4174142 Ry
    -323.4173101 Ry
    -323.4171922 Ry
    -323.4170611 Ry
    -323.4169184 Ry
    -323.4167647 Ry
    -323.4166002 Ry
    -323.4164244 Ry
    -323.4162386 Ry
    -323.4160446 Ry
    -323.4158421 Ry
    -323.4156312 Ry
    -323.4154125 Ry
    -323.4151861 Ry
    -323.4149528 Ry
    -323.4147131 Ry
    -323.414467 Ry
    -323.414215 Ry
    -323.4139583 Ry
    -323.4136953 Ry
    -323.4134285 Ry
    -323.4131559 Ry
    -323.4128797 Ry
    -323.4125984 Ry
    -323.4123147 Ry
    -323.4120269 Ry
    -323.411736 Ry
    -323.4114399 Ry
    -323.4111421 Ry
    -323.4108418 Ry
    -323.4105393 Ry
]

julia>lsqfit(EnergyForm(),BirchMurnaghan3rd(224u"bohr^3", 9u"GPa", 3.75, -161u"Hy"), volumes_u, energies_u)
BirchMurnaghan3rd{Quantity{Float64,D,U} where U where D}(224.4445656763778 a₀^3, 9.194980249913018 GPa, 3.7403684211716297, -161.70885710742223 Eₕ)
```
"""
function lsqfit(
    form::EquationForm,
    eos::EquationOfState{<:Real},
    xdata::AbstractVector{<:Real},
    ydata::AbstractVector{<:Real};
    debug = false,
    kwargs...,
)
    T = promote_type(eltype(xdata), eltype(ydata), Float64)
    E = constructorof(typeof(eos))
    model = (x, p) -> map(apply(form, E(p...)), x)
    fitted = curve_fit(
        model,
        T.(xdata),
        T.(ydata),
        T.(Collections.fieldvalues(eos)),
        kwargs...,
    )
    return debug ? fitted : E(fitted.param...)
end  # function lsqfit
function lsqfit(
    form::EquationForm,
    eos::EquationOfState{<:AbstractQuantity},
    xdata::AbstractVector{<:AbstractQuantity},
    ydata::AbstractVector{<:AbstractQuantity};
    kwargs...,
)
    E = constructorof(typeof(eos))
    values = Collections.fieldvalues(eos)
    original_units = map(unit, values)
    trial_params, xdata, ydata = [map(ustrip ∘ upreferred, x) for x in (
        values,
        xdata,
        ydata,
    )]
    result = lsqfit(form, E(trial_params...), xdata, ydata; kwargs...)
    if result isa EquationOfState
        data = Collections.fieldvalues(result)
        return E([data[i] * upreferred(u) |> u for (i, u) in enumerate(original_units)]...)
    end
    return result
end  # function lsqfit

end
