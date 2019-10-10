"""
# module Find



# Examples

```jldoctest
julia>
```
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

function findvolume(form::EquationForm, eos::EquationOfState, y, x0, method)
    f = v -> apply(form, eos, v) - y
    return find_zero(f, x0, method)
end # function findvolume
function findvolume(form::EquationForm, eos::EquationOfState, y, x0::Union{AbstractVector,Tuple})
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
