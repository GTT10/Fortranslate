# Temperature-dependent liquid density — 2026-09-15

Source: parent `006a86f` plus this change; GNU Fortran 13.3 Debug,
`build/port-spray-debug`.

## Implemented scope

The optional cubic density law follows `MPLiqProps::rho_i` in
`Source/Spray/SprayPropertiesPeleMP.H` at pinned PelePhysics commit
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`. The source was read with `git show`
at that commit; the file SHA256 is
`c75d5f7acc53c58ce2bf8ea56aadd603bef965bfa54b89ca0b817e7cb2a979dd`, also recorded
in the [liquid-property source manifest](../../references/spray_antoine.json).

The runtime evaluates rho=a+T*(b+T*(c+T*d)) with explicit temperature bounds.
Injection, diameter reconstruction, drag relaxation, frozen-film mass loss,
and ETAB/KH-RT share this density. Evaporation holds it fixed over its substep;
the next diameter reflects the updated temperature. Heating without evaporation
changes radius while preserving mass. The liquid energy/heat-capacity model is
unchanged and still neglects p/rho and surface-energy corrections.

The setup checks finite bounded coefficients and an ordered interval; each
actual evaluation requires density between 1e-3 and 1e5 kg/m3. A heating update
that leaves the valid range is rejected before transfer is committed. There is
no fit positivity proof over the whole interval or automatic extrapolation.
The default constant model retains its original context. Polynomial coefficients
and bounds are included in checkpoint identity.

## Focused evidence

- Four core cases passed: new density, existing evaporation, ETAB and KH-RT.
  The density case checks injection mass, mass-preserving thermal expansion,
  exact transfer accounting, both breakup models against constant density at
  the same temperature, nonpositive density rejection and range rollback.
- The CLI suite passed ten tests in 8.13 s. Its new case confirms evaporation,
  coupled conservation, exact gas/parcel restart, changed-coefficient rejection
  and invalid initial density rejection. Five CTest entries passed in 8.17 s.
- `water_density.nml`, an illustrative linear density fixture (1200-0.5*T),
  completed 20 steps to 2e-6 s on a 4^3 grid with four water parcels in nitrogen.
  Liquid mass decreased from 1e-10 to 9.9955843611206153e-11 kg; maximum relative
  system mass error 1.0858505858e-15 and energy error 7.8485860355e-16.
- Provenance preflight, project contract, Python syntax and whitespace passed.
  This is implementation/conservation evidence, not a calibrated liquid fit,
  experiment, boiling model or MPI/AMR particle qualification.

Artifacts: `build/liquid-density-water/` holds the example's outputs, checkpoint,
log and `result.json`; `build/port-spray-debug/density-build.log` records the build.
CLI fixtures use temporary directories.

```sh
cmake --build build/port-spray-debug --target test_spray_core pelef_spray_3d --parallel 8
ctest --test-dir build/port-spray-debug --output-on-failure -R '^spray_(density|evaporation|etab|khrt|application_cli)$'
mkdir -p build/liquid-density-water
cd build/liquid-density-water
../port-spray-debug/pelef_spray_3d ../../cases/spray_3d/water_density.nml
```

Use a fresh output directory. Additional checks from the repository root:

```sh
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```
