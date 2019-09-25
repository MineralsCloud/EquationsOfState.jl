export EquationForm, EnergyForm, PressureForm, BulkModulusForm

abstract type EquationForm end
struct EnergyForm <: EquationForm end
struct PressureForm <: EquationForm end
struct BulkModulusForm <: EquationForm end
