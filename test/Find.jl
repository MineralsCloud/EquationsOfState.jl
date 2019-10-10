using Test

using Roots
using Unitful
using UnitfulAtomic

using EquationsOfState
using EquationsOfState.Collections
using EquationsOfState.Find

# Data in the following tests are from
# https://github.com/materialsproject/pymatgen/blob/1f0957b8525ddc7d12ea348a19caecebe6c7ff34/pymatgen/analysis/tests/test_eos.py
@testset "Test data from Pymatgen" begin
    results = [
        25.987454833,
        26.9045702104,
        27.8430241908,
        28.8029649591,
        29.7848370694,
        30.7887887064,
        31.814968055,
        32.8638196693,
        33.9353435494,
        35.0299842495,
        36.1477417695,
        37.2892088485,
        38.4543854865,
        39.6437162376,
        40.857201102,
        42.095136449,
        43.3579668329,
        44.6456922537,
        45.9587572656,
        47.2973100535,
        48.6614988019,
        50.0517680652,
        51.4682660281,
        52.9112890601,
        54.3808371612,
        55.8775030703,
        57.4014349722,
        58.9526328669,
    ] .* u"angstrom^3"
    energies = [
        -7.63622156576,
        -8.16831294894,
        -8.63871612686,
        -9.05181213218,
        -9.41170988374,
        -9.72238224345,
        -9.98744832526,
        -10.210309552,
        -10.3943401353,
        -10.5427238068,
        -10.6584266073,
        -10.7442240979,
        -10.8027285713,
        -10.8363890521,
        -10.8474912964,
        -10.838157792,
        -10.8103477586,
        -10.7659387815,
        -10.7066179666,
        -10.6339907853,
        -10.5495538639,
        -10.4546677714,
        -10.3506386542,
        -10.2386366017,
        -10.1197772808,
        -9.99504030111,
        -9.86535084973,
        -9.73155247952,
    ] .* u"eV"
    isapprox(
        map(
            e -> findvolume(
                EnergyForm(),
                BirchMurnaghan3rd(
                    40.98926572528106u"angstrom^3",
                    0.5369258245417454u"eV/angstrom^3",
                    4.178644235500821u"1000mm/m",
                    -10.842803908240892u"eV",
                ),
                e,
                (eps(), 100) .* u"angstrom^3",
            ),
            energies,
        ),
        results,
    )
end

@testset "Test `findvolume` with random unit" begin
    pressures = collect(0:20:200) .* u"GPa"
    eos = BirchMurnaghan3rd(167u"angstrom^3", 2600u"kbar", 4.0u"1000mm/m")
    volumes = map(
        p -> findvolume(PressureForm(), eos, p, (eps() * u"bohr^3", eos.v0 * 1.3)),
        pressures,
    )
    @test isapprox(
        ustrip.(map(apply(PressureForm(), eos), volumes) - pressures),
        zeros(11),
        atol = 1e-5,
    )
end # testset
