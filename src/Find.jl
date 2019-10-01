"""
# module Find



# Examples

```jldoctest
julia>
```
"""
module Find

using InteractiveUtils: subtypes
import Statistics
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

function _whose_zero(
    form::EquationForm,
    eos::EquationOfState,
    y::AbstractQuantity,
)
    @assert(eltype(eos) <: AbstractQuantity, "The elements type mismatched!")
    return v::AbstractQuantity -> ustrip(apply(form, eos, v) - y)
end # function _whose_zero
function _whose_zero(
    form::EquationForm,
    eos::EquationOfState,
    y::Real,
)
    @assert(eltype(eos) <: Real, "The elements type mismatched!")
    return v::Real -> apply(form, eos, v) - y
end # function _whose_zero

function _adapt_domain(domain::Union{AbstractVector,Tuple}, method::AbstractBracketing)
    return minimum(domain), maximum(domain)
end # function _adapt_domain
function _adapt_domain(
    domain::Union{AbstractVector,Tuple},
    method::Union{
        AbstractNonBracketing,
        AbstractHalleyLikeMethod,
        AbstractNewtonLikeMethod,
    },
)
    return Statistics.median(domain)
end # function _adapt_domain

function findvolume(
    form::EquationForm,
    eos::EquationOfState,
    y,
    domain::Union{AbstractVector,Tuple},
    method,
)
    f = _whose_zero(form, eos, y)
    return find_zero(f, _adapt_domain(domain), method)
end # function findvolume
function findvolume(
    form::EquationForm,
    eos::EquationOfState,
    y,
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

Statistics.middle(x::AbstractQuantity) = (x + zero(x)) / 1
Statistics.middle(a::T, b::T) where {T<:AbstractQuantity} = middle(ustrip(a), ustrip(b)) * unit(a)
function Statistics.middle(a::AbstractQuantity, b::AbstractQuantity)
    @assert(dimension(a) == dimension(b))
    a0, b0 = promote(map(ustrip, (a, b))...)
    a, b = a0 * unit(a), b0 * unit(b)
    a, b = map(upreferred, a, b)
    return middle(a, b)
end # Statistics.middle

end
