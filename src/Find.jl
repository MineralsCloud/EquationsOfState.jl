"""
# module Find



# Examples

```jldoctest
julia>
```
"""
module Find

using InteractiveUtils: subtypes
using Unitful: AbstractQuantity, ustrip, upreferred

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

function _whose_zero(form::EquationForm, eos::EquationOfState, y::AbstractQuantity)
    @assert(eltype(eos) <: AbstractQuantity, "The elements type mismatched!")
    return v::AbstractQuantity -> ustrip(apply(form, eos, v) - y)
end # function _whose_zero
function _whose_zero(form::EquationForm, eos::EquationOfState, y::Real)
    @assert(eltype(eos) <: Real, "The elements type mismatched!")
    return v::Real -> apply(form, eos, v) - y
end # function _whose_zero

function findvolume(form::EquationForm, eos::EquationOfState, y, x0, method)
    f = _whose_zero(form, eos, y)
    return find_zero(f, x0, method)
end # function findvolume
function findvolume(form::EquationForm, eos::EquationOfState, y, x0)
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
            return findvolume(form, eos, y, x0, T())
        catch e
            @info("Method \"$T\" failed because of $e.")
            continue
        end
    end
end # function findvolume

end
