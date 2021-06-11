module EquationsOfState

abstract type EquationOfStateOfSolidsParameters{T} end
const Parameters = EquationOfStateOfSolidsParameters

abstract type EquationOfState{T<:Parameters} end

end
