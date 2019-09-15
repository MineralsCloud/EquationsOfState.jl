export EquationOfStateForm, EnergyForm, PressureForm, BulkModulusForm

abstract type EquationOfStateForm end
struct EnergyForm <: EquationOfStateForm end
struct PressureForm <: EquationOfStateForm end
struct BulkModulusForm <: EquationOfStateForm end
