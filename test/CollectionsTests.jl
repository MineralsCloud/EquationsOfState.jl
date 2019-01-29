#=
CollectionsTests:
- Julia version: 1.0
- Author: qz
- Date: Jan 29, 2019
=#
module CollectionsTests

using Test

using EquationsOfState

@testset "Test EOS promotion" begin
    typeof(Birch(1, 2, 3.0)) == Birch{Float64}
    typeof(Murnaghan(1, 2, 3.0)) == Murnaghan{Float64}
    typeof(BirchMurnaghan2nd(1, 2.0)) == BirchMurnaghan2nd{Float64}
    typeof(BirchMurnaghan3rd(1, 2, 3.0)) == BirchMurnaghan3rd{Float64}
    typeof(BirchMurnaghan4th(1, 2.0, 3, 4)) == BirchMurnaghan4th{Float64}
    typeof(Vinet(1, 2, 3.0)) == Vinet{Float64}
    typeof(PoirierTarantola2nd(1, 2.0)) == PoirierTarantola2nd{Float64}
    typeof(PoirierTarantola3rd(1, 2, 3.0)) == PoirierTarantola3rd{Float64}
    typeof(PoirierTarantola4th(1, 2, 3, 4)) == PoirierTarantola4th{Int}
    typeof(Holzapfel(1, 2, 3.0, 4.0)) == Holzapfel{Float64}
    typeof(AntonSchmidt(1, 2, 3.0)) == AntonSchmidt{Float64}
    typeof(BreenanStacey(1, 2, 3.0)) == BreenanStacey{Float64}
end

end