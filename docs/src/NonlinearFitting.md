# Nonlinear fitting

```@meta
CurrentModule = EquationsOfState.NonlinearFitting
```

From Ref. 1,

> The equations of state depend nonlinearly of a collection of parameters, $E_0$, $V_0$, $B_0$, $B_0'$, ..., that represent physical properties of the solid at equilibrium and can, in principle, be obtained expermentally by independent methods. The use of a given analytical EOS may have significant influence on the results obtained, particularly because the parameters are far from being independent. The number of parameters has to be considered in comparing the goodness of fit of functional forms with different analytical flexibility. The possibility of using too many parameters, beyond what is physically justified by the information contained in the experimental data, is a serious aspect that deserves consideration.

In [`EquationsOfState`](https://github.com/MineralsCloud/EquationsOfState.jl), the nonlinear fitting is currently implemented by [`LsqFit`](https://github.com/JuliaNLSolvers/LsqFit.jl), a small library that provides basic least-squares fitting in pure Julia. It only utilizes the *Levenberg-Marquardt algorithm* for non-linear fitting. See its [documentation](https://github.com/JuliaNLSolvers/LsqFit.jl/blob/master/README.md) for more information.

## Usage

We provide API `lsqfit` currently.

```julia
using EquationsOfState
using EquationsOfState.Collections
using EquationsOfState.NonlinearFitting

volumes = [
    25.987454833,
    26.9045702104,
    27.8430241908,
    28.8029649591,
    29.7848370694,
    30.7887887064,
    31.814968055,
    32.8638196693,
    33.9353435494,
    35.0299842495,
    36.1477417695,
    37.2892088485,
    38.4543854865,
    39.6437162376,
    40.857201102,
    42.095136449,
    43.3579668329,
    44.6456922537,
    45.9587572656,
    47.2973100535,
    48.6614988019,
    50.0517680652,
    51.4682660281,
    52.9112890601,
    54.3808371612,
    55.8775030703,
    57.4014349722,
    58.9526328669
]
energies = [
    -7.63622156576,
    -8.16831294894,
    -8.63871612686,
    -9.05181213218,
    -9.41170988374,
    -9.72238224345,
    -9.98744832526,
    -10.210309552,
    -10.3943401353,
    -10.5427238068,
    -10.6584266073,
    -10.7442240979,
    -10.8027285713,
    -10.8363890521,
    -10.8474912964,
    -10.838157792,
    -10.8103477586,
    -10.7659387815,
    -10.7066179666,
    -10.6339907853,
    -10.5495538639,
    -10.4546677714,
    -10.3506386542,
    -10.2386366017,
    -10.1197772808,
    -9.99504030111,
    -9.86535084973,
    -9.73155247952
]

lsqfit(EnergyForm(), BirchMurnaghan3rd(40, 0.5, 4, 0), volumes, energies)
lsqfit(EnergyForm(), Murnaghan(41, 0.5, 4, 0), volumes, energies)
lsqfit(EnergyForm(), PoirierTarantola3rd(41, 0.5, 4, 0), volumes, energies)
lsqfit(EnergyForm(), Vinet(41, 0.5, 4, 0), volumes, energies)
```
Then 4 different equations of state will be fitted.

## Public interfaces

```@docs
lsqfit(::EquationForm, eos::E, xdata::X, ydata::Y; debug = false, kwargs...) where {E<:EquationOfState,X<:AbstractVector,Y<:AbstractVector}
```


## References

1. [A. Otero-De-La-Roza, V. Luaña, *Computer Physics Communications*. **182**, 1708–1720 (2011), doi:10.1016/j.cpc.2011.04.016.](https://www.sciencedirect.com/science/article/pii/S0010465511001470)
