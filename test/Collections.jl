using Test

using Unitful

using EquationsOfState.Collections

@testset "Test EOS promotion" begin
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
    @test Murnaghan(1, Int32(2), Int8(3), 0) == Murnaghan{Int}(1, 2, 3, 0)
    @test Murnaghan(1, 2//1, Int8(3), 0) == Murnaghan{Rational{Int}}(1//1, 2//1, 3//1, 0//1)
    @test typeof(Murnaghan(1u"angstrom^3", 2u"eV/angstrom^3", 3.0, 4u"eV")) == Murnaghan{Quantity{Float64}}
    @test typeof(Murnaghan(1u"angstrom^3", 2u"eV/angstrom^3", 3//2, 4u"eV")) == Murnaghan{Quantity{Rational{Int}}}
end

@testset "Test default EOS parameter `e0` and promotion" begin
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
    @test typeof(Murnaghan(1u"angstrom^3", 2u"eV/angstrom^3", 3)) == Murnaghan{Quantity{Int}}
    @test typeof(Murnaghan(1u"angstrom^3", 2u"eV/angstrom^3", 3.0)) == Murnaghan{Quantity{Float64}}
    @test typeof(Murnaghan(1.0u"angstrom^3", 2u"eV/angstrom^3", 3.0)) == Murnaghan{Quantity{Float64}}
end
