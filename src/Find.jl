"""
    findvolume(form, eos, y, x, method)

Calculate volumes based on crude eos and input pressures

#Arguments
- `form::EquationForm`: an `EquationForm` instance. If `EnergyForm`, fit ``E(V)``; if `PressureForm`, fit ``P(V)``; if `BulkModulusForm`, fit ``B(V)``.
- `eos::EquationOfState`: a trial equation of state.
- `y::Real`: a real number of pressure value.
- `x0::Union{AbstractVector,Tuple}`: an AbstractVector or Tuple of volume range.
- `method <:AbstractUnivariateZeroMethod`: a method used to do calculation.

# Examples

```jldoctest

    # Data in the following tests are from
    # https://github.com/materialsproject/pymatgen/blob/1f0957b8525ddc7d12ea348a19caecebe6c7ff34/pymatgen/analysis/tests/test_eos.py

julia>pressures = collect(0:20:200) .* u"GPa"
julia>eos = BirchMurnaghan3rd(167u"angstrom^3", 2600u"kbar", 4.0u"1000mm/m")
julia>volumes = map(
        p -> findvolume(PressureForm(), eos, p, (eps() * u"bohr^3", eos.v0 * 1.3)),
        pressures,
    )

11-element Array{Quantity{Float64,�^3,Unitful.FreeUnits{(Å^3,),�^3,nothing}},1}:
167.0 Å^3
156.14036210727835 Å^3
147.99803635986558 Å^3
141.5109371379587 Å^3
136.13864615965332 Å^3
131.5678403193935 Å^3
127.60046278645824 Å^3
124.10332447387113 Å^3
120.9825768060646 Å^3
118.16962836248427 Å^3
115.61284838696815 Å^3
"""sf
module Find

using InteractiveUtils: subtypes
using Unitful: AbstractQuantity, ustrip

using Roots: find_zero,
             AbstractBracketing,
             AbstractNonBracketing,
             AbstractHalleyLikeMethod,
             AbstractNewtonLikeMethod,
             AbstractAlefeldPotraShi,
             AbstractBisection,
             AbstractSecant,
             Brent,
             Newton,
             ConvergenceFailed

import ..EquationForm
using ..Collections: EquationOfState, apply

export findvolume

function findvolume(form::EquationForm, eos::EquationOfState, y, x0, method)
    f = v -> apply(form, eos, v) - y
    return find_zero(f, x0, method)
end # function findvolume
function findvolume(
    form::EquationForm,
    eos::EquationOfState,
    y,
    x0::Union{AbstractVector,Tuple},
)
    for T in [subtypes(AbstractBisection); subtypes(AbstractAlefeldPotraShi)]
        @info("Using method \"$T\"...")
        try
            # `maximum` and `minimum` also works with `AbstractQuantity`s.
            return findvolume(form, eos, y, (minimum(x0), maximum(x0)), T())
        catch e
            @info("Method \"$T\" failed because of $e.")
            continue
        end
    end
    for T in [
        Brent
        subtypes(AbstractHalleyLikeMethod)
        Newton
        subtypes(AbstractSecant)
    ]
        @info("Using method \"$T\"...")
        try
            return findvolume(form, eos, y, (minimum(x0) + maximum(x0)) / 2, T())
        catch e
            @info("Method \"$T\" failed because of $e.")
            continue
        end
    end
end # function findvolume

end
