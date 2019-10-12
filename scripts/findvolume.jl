# This file plots the `findvolume.png` in docs
using Unitful, UnitfulAtomic, Plots
using EquationsOfState, EquationsOfState.Collections, EquationsOfState.Find

plotlyjs()

pressures = collect(0:20:200) .* u"GPa"
eos = BirchMurnaghan3rd(167 * u"angstrom^3", 2600 * u"kbar", 4.0)
volumes = map(
    p -> findvolume(PressureForm(), eos, p, (eps() * u"bohr^3", eos.v0 * 1.3)),
    pressures,
)
plot(ustrip.(volumes), ustrip.(pressures), label = "pressures")
scatter!(
    ustrip.(volumes),
    ustrip.(u"GPa", apply(PressureForm(), eos).(volumes)),
    label = "P(V)",
    dpi = 400,
)
xlabel!("volume (angstrom³)")
ylabel!("pressure (GPa)")
savefig("docs/src/assets/findvolume.png")
