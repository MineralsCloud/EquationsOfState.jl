"""
This module provides `findvolume` methods to find the volume at a given
pressure, energy, or bulk modulus with(out) units.
"""
module Find

using InteractiveUtils: subtypes
using Roots:
    find_zero,
    AbstractHalleyLikeMethod,
    AbstractAlefeldPotraShi,
    AbstractBisection,
    AbstractSecant,
    Brent,
    Newton
using Unitful: AbstractQuantity, ustrip

using ..Collections: EquationOnVolume, PhysicalProperty

export findvolume

"""
    findvolume(eos(prop), y, x0, method)
    findvolume(eos(prop), y, x0::Union{AbstractVector,Tuple})

Find a volume which leads to the given pressure, energy, or bulk modulus based on an `eos`.

# Arguments
- `eos::EquationOfState`: an equation of state. If it has units, `y` and `x0` must also have.
- `prop::PhysicalProperty`: a `PhysicalProperty` instance.
- `y`: a pressure, energy, or bulk modulus.
- `x0`: can be either a range of volumes (`Vector`, `Tuple`, etc.) or just a single volume.
    Units can be provided if necessary.
- `method::Roots.AbstractUnivariateZeroMethod`: a method used to find the root of an equation.
    If it is omitted, the algorithm will traverse all possible methods of
    [Roots.jl](https://github.com/JuliaMath/Roots.jl). And the `x0` parameter must be
    an array or a tuple, of which only the maximum and minimum values will be used in the
    root-finding process.
"""
findvolume(f::EquationOnVolume, y, x0, method) = find_zero(v -> f(v) - y, x0, method)
function findvolume(f::EquationOnVolume, y, x0::Union{AbstractVector,Tuple})
    for T in [subtypes(AbstractBisection); subtypes(AbstractAlefeldPotraShi)]
        @info("Using method \"$T\"...")
        try
            # `maximum` and `minimum` also works with `AbstractQuantity`s.
            return findvolume(f, y, (minimum(x0), maximum(x0)), T())
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
            return findvolume(f, y, (minimum(x0) + maximum(x0)) / 2, T())
        catch e
            @info("Method \"$T\" failed because of $e.")
            continue
        end
    end
end # function findvolume

end
