"""
# module Find



# Examples

```jldoctest
julia>
```
"""
module Find

using InteractiveUtils: subtypes
using Statistics: median

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

function findvolume(
    form::EquationForm,
    eos::EquationOfState,
    y::Real,
    domain::Union{AbstractVector,Tuple},
    method::AbstractBracketing,
)
    f(v) = apply(form, eos, v) - y
    return find_zero(f, (minimum(domain), maximum(domain)), method)
end # function findvolume
function findvolume(
    form::EquationForm,
    eos::EquationOfState,
    y::Real,
    domain::Union{AbstractVector,Tuple},
    method::Union{
        AbstractNonBracketing,
        AbstractHalleyLikeMethod,
        AbstractNewtonLikeMethod,
    },
)
    f(v) = apply(form, eos, v) - y
    return find_zero(f, median(domain), method)
end # function findvolume
function findvolume(
    form::EquationForm,
    eos::EquationOfState,
    y::Real,
    domain::Union{AbstractVector,Tuple},
)
    for T in [
        subtypes(AbstractAlefeldPotraShi)
        subtypes(AbstractBisection)
        Brent
        subtypes(AbstractHalleyLikeMethod)
        Newton
        subtypes(AbstractSecant)
    ]
        @info("Using method \"$T\"...")
        try
            return findvolume(form, eos, y, domain, T())
        catch e
            @info("Method \"$T\" failed because of $e.")
            continue
        end
    end
end # function findvolume

end
