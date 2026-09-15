# Temperature-dependent liquid viscosity — 2026-09-15

Source: parent `0669b8b` plus this change; GNU Fortran 13.3 Debug,
`build/port-spray-debug`.

## Implemented scope

The optional inverse-temperature polynomial follows `MPLiqProps::mu_i` in
`Source/Spray/SprayPropertiesPeleMP.H` at pinned PelePhysics commit
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`, inspected with `git show`.
The file SHA256 is
`c75d5f7acc53c58ce2bf8ea56aadd603bef965bfa54b89ca0b817e7cb2a979dd`, also recorded
in the [liquid-property source manifest](../../references/spray_antoine.json).

ETAB damping and KH-RT's Ohnesorge number now use mu=a+b/T+c/T^2+d/T^3 at the
current droplet temperature when selected. Temperature remains fixed within
an individual breakup call. Both models already use temperature-dependent
density. Setup validates finite coefficients and an ordered interval; each
actual viscosity evaluation must be within 1e-10--1e3 Pa s. Out-of-range
states reject the breakup/coupling transaction. Positivity throughout the
entire fit interval is not inferred from endpoint or setup checks.

The fit requires breakup to be enabled, and its coefficients/range participate
in restart identity. Constant-viscosity behavior and checkpoint context are
retained by default. This does not add liquid thermal conduction, a surface
tension fit, boiling, wall impact or multicomponent droplets.

## Focused evidence

- Four core cases passed: viscosity, density, ETAB and KH-RT. The new case
  checks eight temperatures against the expanded formula, both breakup models
  against their constant-viscosity evaluation at the same temperature, range
  rollback, negative viscosity and NaN rejection.
- The CLI suite passed eleven tests in 8.61 s, including combined evaporation,
  KH-RT, conservation, exact gas/parcel restart, changed viscosity coefficient
  rejection and invalid initial fit-range rejection. Five CTest entries passed
  in 8.65 s.
- The supplied combined density/viscosity case completed 20 steps to 2e-6 s
  on a 4^3 grid with four water parcels, 900 K nitrogen at 4 MPa and a 200 m/s
  injection. Sauter mean diameter decreased from 80 to 70.4277398 micrometers.
  Liquid mass changed from 1e-10 to 9.9999110724570894e-11 kg; most of the
  diameter change is representative breakup, whose multiplicity conserves
  liquid mass. Maximum relative system mass and energy errors were
  6.6279291715e-16 and 1.0620497517e-15.
- Provenance preflight, project contract, Python syntax and whitespace passed.
  Coefficients are illustrative; no property fit, physical breakup experiment
  or MPI/AMR-particle qualification was performed.

Artifacts: `build/liquid-viscosity-water/` contains the example output,
checkpoint, log and `result.json`; build log is
`build/port-spray-debug/viscosity-build.log`. CLI fixtures use temporary directories.

```sh
cmake --build build/port-spray-debug --target test_spray_core pelef_spray_3d --parallel 8
ctest --test-dir build/port-spray-debug --output-on-failure -R '^spray_(viscosity|density|etab|khrt|application_cli)$'
mkdir -p build/liquid-viscosity-water
cd build/liquid-viscosity-water
../port-spray-debug/pelef_spray_3d ../../cases/spray_3d/water_viscosity.nml
```

Use a fresh output directory. Additional checks from the repository root:

```sh
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```
