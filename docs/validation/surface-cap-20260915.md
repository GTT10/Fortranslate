# Near-boiling vapor-state temperature cap — 2026-09-15

Source: parent `8f576ce` plus this change; GNU Fortran 13.3 Debug,
`build/port-spray-debug`. Pinned source files and SHA256 identities are in
[`references/spray_surface_cap.json`](../../references/spray_surface_cap.json).

## Implemented scope

The pinned PelePhysics `Drag.H` limits the liquid temperature used to compute
vapor pressure to 0.99 of the pressure-dependent boiling temperature.
`SprayFuelData.H::calcBoilT` computes that temperature using Clausius--Clapeyron
and caps it at the critical temperature. Both were read at commit
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b` using `git show`, not the external
checkout's current HEAD. That pinned update does not delete droplets merely
because the updated temperature exceeds the boiling point.

The new optional `surface-cap` path uses the same 0.99 convention with the
existing constant reference latent heat. It changes only the temperature used
for saturation pressure, keeping the conservative parcel-to-gas exchange.
The actual droplet temperature, liquid energy, film transport temperature and
property-fit limits are not overwritten. Setup requires the Clausius law and
Tc>Tref. A nonpositive inverse boiling temperature fails explicitly.

Default `reject` behavior and checkpoint context remain; the optional mode and
critical temperature participate in restart identity. This is not a full
physical boiling model. The upstream Watson latent-property correlation,
real-fluid behavior and multicomponent vapor-state mapping remain outside this
increment.

## Focused evidence

- Three core cases passed: surface cap, Antoine and evaporation. The new case
  checks default rejection at 400 K/101325 Pa, equality with evaluation at
  0.99*Tref, pressure dependence over 0.1/1/10 reference atmospheres, and
  invalid critical temperature rejection.
- Twelve CLI cases passed in 9.08 s, including evaporation from 400 K droplets,
  coupled conservation, exact gas/parcel restart, changed critical-temperature
  rejection and failure of the same initial temperature under the default law.
  Four CTest entries passed in 9.11 s. After placing new type components last
  and generalizing an input error message, the build and core cap test passed
  again (0.02 s); no broad matrix was rerun.
- The supplied 4^3-grid/four-parcel example completed 20 steps to 2e-6 s with
  400 K water droplets in 900 K nitrogen at 101325 Pa. Liquid mass decreased
  from 1e-10 to 9.9534331885794867e-11 kg. Maximum relative system mass and
  energy errors were 1.2215819090e-15 and 1.4520968824e-15.
- Mechanism provenance, project contract, Python syntax and whitespace passed.
  No physical boiling experiment, MPI/AMR particle or supercritical claim is
  supported by this result.

Artifacts: `build/water-surface-cap/` contains the example log, outputs,
checkpoint and `result.json`; build logs are
`build/port-spray-debug/surface-cap-build.log` and `surface-cap-final-build.log`.
CLI fixtures use temporary directories.

```sh
cmake --build build/port-spray-debug --target test_spray_core pelef_spray_3d --parallel 8
ctest --test-dir build/port-spray-debug --output-on-failure -R '^spray_(surface_cap|antoine|evaporation|application_cli)$'
mkdir -p build/water-surface-cap
cd build/water-surface-cap
../port-spray-debug/pelef_spray_3d ../../cases/spray_3d/water_surface_cap.nml
```

Use a fresh output directory. Additional checks from the repository root:

```sh
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```
