# EquationsOfState.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://MineralsCloud.github.io/EquationsOfState.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://MineralsCloud.github.io/EquationsOfState.jl/dev)
[![Build Status](https://travis-ci.com/MineralsCloud/EquationsOfState.jl.svg?branch=master)](https://travis-ci.com/MineralsCloud/EquationsOfState.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/singularitti/EquationsOfState.jl?svg=true)](https://ci.appveyor.com/project/singularitti/EquationsOfState-jl)
[![Codecov](https://codecov.io/gh/MineralsCloud/EquationsOfState.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/MineralsCloud/EquationsOfState.jl)
[![Coveralls](https://coveralls.io/repos/github/MineralsCloud/EquationsOfState.jl/badge.svg?branch=master)](https://coveralls.io/github/MineralsCloud/EquationsOfState.jl?branch=master)
[![Build Status](https://api.cirrus-ci.com/github/MineralsCloud/EquationsOfState.jl.svg)](https://cirrus-ci.com/github/MineralsCloud/EquationsOfState.jl)

This package implements some _equations of state_ (EOS) of solids which are useful in research. It currently includes:

1. `Birch` EOS
2. `Murnaghan` EOS
3. Birch–Murnaghan EOS family:
    1. `BirchMurnaghan2nd`
    2. `BirchMurnaghan3rd`
    3. `BirchMurnaghan4th`
4. `Vinet` EOS
5. Poirier–Tarantola EOS family:
    1. `PoirierTarantola2nd`
    2. `PoirierTarantola3rd`
    3. `PoirierTarantola4th`
6. `AntonSchmidt` EOS (experimental)
7. `BreenanStacey` EOS (experimental)

The formula are referenced from Ref. 1.

This package also includes linear and nonlinear fitting methods, also referenced from Ref. 1.

It is built for Julia v1.0+.

## TODOs

- [ ] Implement nonlinear fitting using [CMPFit.jl](https://github.com/gcalderone/CMPFit.jl).

## Related packages

1. [CommandLineEquationsOfState.jl](https://github.com/MineralsCloud/CommandLineEquationsOfState.jl)
2. [ExtendedEquationsOfState.jl](https://github.com/MineralsCloud/ExtendedEquationsOfState.jl)

## References

1. A. Otero-De-La-Roza, V. Luaña, *Comput. Phys. Commun.* **182**, 1708–1720 (2011).
