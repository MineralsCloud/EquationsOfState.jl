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
    @test typeof(Birch(1, 2, 3.0, 0)) == Birch{Float64}
    @test typeof(Murnaghan(1, 2, 3.0, 0)) == Murnaghan{Float64}
    @test typeof(BirchMurnaghan2nd(1, 2.0, 0)) == BirchMurnaghan2nd{Float64}
    @test typeof(BirchMurnaghan3rd(1, 2, 3.0, 0)) == BirchMurnaghan3rd{Float64}
    @test typeof(BirchMurnaghan4th(1, 2.0, 3, 4, 0)) == BirchMurnaghan4th{Float64}
    @test typeof(Vinet(1, 2, 3.0, 0)) == Vinet{Float64}
    @test typeof(PoirierTarantola2nd(1, 2.0, 0)) == PoirierTarantola2nd{Float64}
    @test typeof(PoirierTarantola3rd(1, 2, 3.0, 0)) == PoirierTarantola3rd{Float64}
    @test typeof(PoirierTarantola4th(1, 2, 3, 4, 0)) == PoirierTarantola4th{Int}
    @test typeof(AntonSchmidt(1, 2, 3.0, 0)) == AntonSchmidt{Float64}
    @test typeof(BreenanStacey(1, 2, 3.0, 0)) == BreenanStacey{Float64}
end

@testset "Test default EOS parameter `e0` and promotion" begin
    @test Birch(1, 2, 3.0) == Birch(1.0, 2.0, 3.0, 0.0)
    @test Murnaghan(1, 2, 3.0) == Murnaghan(1.0, 2.0, 3.0, 0.0)
    @test BirchMurnaghan2nd(1, 2.0) == BirchMurnaghan2nd(1.0, 2.0, 0.0)
    @test BirchMurnaghan3rd(1, 2, 3.0) == BirchMurnaghan3rd(1.0, 2.0, 3.0, 0.0)
    @test BirchMurnaghan4th(1, 2.0, 3, 4) == BirchMurnaghan4th(1.0, 2.0, 3.0, 4.0, 0.0)
    @test Vinet(1, 2, 3.0) == Vinet(1.0, 2.0, 3.0, 0.0)
    @test PoirierTarantola2nd(1, 2.0) == PoirierTarantola2nd(1.0, 2.0, 0.0)
    @test PoirierTarantola3rd(1, 2, 3.0) == PoirierTarantola3rd(1.0, 2.0, 3.0, 0.0)
    @test PoirierTarantola4th(1, 2, 3, 4) == PoirierTarantola4th(1, 2, 3, 4, 0)
    @test AntonSchmidt(1, 2, 3.0) == AntonSchmidt(1.0, 2.0, 3.0, 0.0)
    @test BreenanStacey(1, 2, 3.0) == BreenanStacey(1.0, 2.0, 3.0, 0.0)
end

@testset "Test collecting parameters" begin
    @test collect(Birch(1, 2, 3.0, 0)) == [1.0, 2.0, 3.0, 0.0]
    @test collect(Murnaghan(1, 2, 3.0, 0)) == [1.0, 2.0, 3.0, 0.0]
    @test collect(BirchMurnaghan2nd(1, 2.0, 0)) == [1.0, 2.0, 0.0]
    @test collect(BirchMurnaghan3rd(1, 2, 3.0, 0)) == [1.0, 2.0, 3.0, 0.0]
    @test collect(BirchMurnaghan4th(1, 2.0, 3, 4, 0)) == [1.0, 2.0, 3.0, 4.0, 0.0]
    @test collect(Vinet(1, 2, 3.0, 0)) == [1.0, 2.0, 3.0, 0.0]
    @test collect(PoirierTarantola2nd(1, 2.0, 0)) == [1.0, 2.0, 0.0]
    @test collect(PoirierTarantola3rd(1, 2, 3.0, 0)) == [1.0, 2.0, 3.0, 0.0]
    @test collect(PoirierTarantola4th(1, 2, 3, 4, 0)) == [1, 2, 3, 4, 0]
    @test collect(AntonSchmidt(1, 2, 3.0, 0)) == [1.0, 2.0, 3.0, 0.0]
    @test collect(BreenanStacey(1, 2, 3.0, 0)) == [1.0, 2.0, 3.0, 0.0]
end

@testset "Test converting elements' type" begin
    @test convert(Birch{Float64}, Birch(1, 2, 3, 0)) == Birch(1.0, 2.0, 3.0, 0.0)
    @test convert(Murnaghan{Float64}, Murnaghan(1, 2, 3, 0)) == Murnaghan(1.0, 2.0, 3.0, 0.0)
    @test convert(BirchMurnaghan2nd{Float64}, BirchMurnaghan2nd(1, 2, 0)) == BirchMurnaghan2nd(1.0, 2.0, 0.0)
    @test convert(BirchMurnaghan3rd{Float64}, BirchMurnaghan3rd(1, 2, 3, 0)) == BirchMurnaghan3rd(1.0, 2.0, 3.0, 0.0)
    @test convert(BirchMurnaghan4th{Float64}, BirchMurnaghan4th(1, 2, 3, 4, 0)) == BirchMurnaghan4th(1.0, 2.0, 3.0, 4.0, 0.0)
    @test convert(Vinet{Float64}, Vinet(1, 2, 3, 0)) == Vinet(1.0, 2.0, 3.0, 0.0)
    @test convert(PoirierTarantola2nd{Float64}, PoirierTarantola2nd(1, 2, 0)) == PoirierTarantola2nd(1.0, 2.0, 0.0)
    @test convert(PoirierTarantola3rd{Float64}, PoirierTarantola3rd(1, 2, 3, 0)) == PoirierTarantola3rd(1.0, 2.0, 3.0, 0.0)
    @test convert(PoirierTarantola4th{Float64}, PoirierTarantola4th(1, 2, 3, 4, 0)) == PoirierTarantola4th(1, 2, 3, 4, 0)
    @test convert(AntonSchmidt{Float64}, AntonSchmidt(1, 2, 3, 0)) == AntonSchmidt(1.0, 2.0, 3.0, 0.0)
    @test convert(BreenanStacey{Float64}, BreenanStacey(1, 2, 3, 0)) == BreenanStacey(1.0, 2.0, 3.0, 0.0)
end

end
