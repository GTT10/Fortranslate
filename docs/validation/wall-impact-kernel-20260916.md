# Ahamed wall-impact kernel — 2026-09-16

Source: parent `8477871` plus this change; GNU Fortran 13.3 Debug,
`build/port-spray-debug`, MPI disabled.

## Scope

The complete six-regime classification and the reference splash-velocity and
azimuthal mass-fraction correlations are ported to a standalone Fortran module.
The application is not connected yet: wall-film storage, impact geometry,
child-parcel creation and conservative wall exchange/restart remain required.
No wall-impact or spray physical validation is claimed.

PelePhysics commit `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`,
`Source/Spray/BreakupSplash/AhamedSplash.H`, SHA-256
`37ddd95a06bc96d6c2b1e740cf583a872335efc787b56a4cce80211f4381cb92`.
The existing pinned spray license applies; provenance is in
`references/spray_wall_impact.json`.

## Focused evidence

The pinned C++ functions `splash_criteria`, `get_splash_vels`, `get_ms_theta`
were extracted unchanged into a C++17 wrapper with scalar AMReX type aliases
and macro/enum stubs. The compiled output supplies the fixed Fortran comparison
values for dry and wet cases. For Kv=150, Tstar=1, alpha=pi/6, Kcrit=80,
film ratios 0/.03/.07/.2 produced regimes 2/5/5/4 and splashed fractions
0.14206181552235969/0.26793114657822259/0.66775589161253757/0.
The native test compares dry/thin/intermediate-film velocity components and
dry/thin-film mass distributions, checks all six regimes and rejects invalid
inputs. One CTest entry passed in 0.02 s, zero failures/skips.

The reference thick-film distribution at alpha=pi/6, splashed fraction=.5,
film ratio=.2 returns a negative fourth fraction (-0.13034999999999991).
The Fortran routine explicitly rejects that candidate and zeros mass outputs.
It also rejects a negative classified splashed fraction rather than silently
renormalizing it. This is a declared admissibility boundary, not reference
parity for unphysical outputs or completion of wall-film dynamics.

Build, fixed-mechanism provenance, project contract and whitespace checks
passed. No Python source changed. No broad matrix was run.

Artifacts: `build/wall-impact-reference/{AhamedSplash.H,oracle.cpp,oracle,oracle.txt}`;
`build/port-spray-debug/wall-impact-build.log` and
`build/port-spray-debug/Testing/Temporary/LastTest.log`.

```sh
c++ -std=c++17 build/wall-impact-reference/oracle.cpp -o build/wall-impact-reference/oracle
build/wall-impact-reference/oracle
cmake --build build/port-spray-debug --target test_spray_core --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_wall_impact$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```
