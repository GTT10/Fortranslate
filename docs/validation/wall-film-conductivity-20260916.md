# Wall-film conductivity — 2026-09-16

## Change and source

Parent: `ec7116a27d308c3a7a7c0494d05921ef68a4d1b6`, plus the conductivity
implementation and tests in this commit. Exact tested source SHA-256 values,
compiler identity and logs are retained in `build/wall-film-conductivity-checks/`.

The cubic `lambda_i` expression comes from PelePhysics
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`,
`Source/Spray/SprayPropertiesPeleMP.H` (SHA-256
`c75d5f7acc53c58ce2bf8ea56aadd603bef965bfa54b89ca0b817e7cb2a979dd`).
The reference identity and license are recorded in `references/spray_wall_impact.json`.

Film conductivity now supports a global constant or cubic and per-port model
and scalar overrides. Overrides select the deposited film's material ID.
The complete fit interval is checked at endpoints and stationary points;
substeps evaluate conductivity at current film temperature and reject an
out-of-range candidate before committing. Restart identity includes the new
settings, and loaded film temperatures must lie in the selected fit domain.
Legacy constant settings retain their existing checkpoint context.

## Focused verification

GNU Fortran 13.3.0, existing `build/port-spray-debug`, Debug, spray/tests enabled,
MPI disabled. Commands completed:

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d test_spray_core test_spray_integration --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(film_evolution|application_cli|integration_wall_coupling)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Three CTest entries passed in 27.06 s, including 22 CLI unittests; zero failures
and skips. Build completed with real-comparison and existing bounds-analysis
warnings. Provenance, Python syntax and project contract passed.

The added checks cover the cubic value, per-material constant override, invalid
material/temperature, and a polynomial with positive endpoints but a negative
interior. Application cases check global and per-port fits, distinct film
heating for two otherwise identical liquids, conservation, byte-identical final
gas/parcel/film CSVs after restart, final history equality, changed-coefficient
restart rejection, out-of-domain temperature and unused-port rejection.
The existing film evolution and wall coupling tests also passed.

## Remaining boundary

This is a stationary pure-liquid film model with fixed footprint. Conductivity
uses SI coefficients and film temperature; reference mixture averaging and
boiling-temperature clamping are not implemented. Illustrative coefficients
are not calibrated liquid properties. Lateral film flow, moving walls,
entrainment and physical spray qualification remain open. This focused check
does not establish full PeleC parity or full-suite qualification.
