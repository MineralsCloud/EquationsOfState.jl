"""
This module provides `findvolume` methods to find the volume at a given
pressure, energy, or bulk modulus with(out) units.
"""
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

"""
    findvolume(form, eos, y, x0, method)
    findvolume(form, eos, y, x0::Union{AbstractVector,Tuple})

Find a volume which leads to the given pressure, energy, or bulk modulus based on an `eos`.

# Arguments
- `form::EquationForm`: an `EquationForm` instance.
- `eos::EquationOfState`: an equation of state. If it has units, `y` and `x0` must also have.
- `y`: a pressure, energy, or bulk modulus.
- `x0`: can be either a range of volumes (`Vector`, `Tuple`, etc.) or just a single volume.
    Units can be provided if necessary.
- `method::Roots.AbstractUnivariateZeroMethod`: a method used to find the root of an equation.
    If it is omitted, the algorithm will traverse all possible methods of 
    [Roots.jl](https://github.com/JuliaMath/Roots.jl). And the `x0` parameter must be
    an array or a tuple, of which only the maximum and minimum values will be used in the
    root-finding process.
"""
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
