# Mixed-liquid planar wall crossing — 2026-09-16

Tested parent: `f3f7df0` plus this change. Source/test hashes and build/test/preflight
logs are in `build/mixture-crossing-checks/`.

## Implementation and source

The mixed dispatcher uses the earliest crossed box face, current mixture density,
prescribed viscosity/surface tension, existing film geometry and local pressure to
compute impact parameters. It preserves component masses in retained/splash drops
and deposition, using per-drop and whole-film kg respectively. A failed event leaves
all supplied state unchanged.

Pinned PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`,
`SprayPropertiesPeleMP.H` lines 124–128, defines `boilT_mix` as the mass-weighted
component boiling temperatures. `AhamedSplash.H` lines 231–236 uses that value for
the wall-temperature ratio. These source locations were read from the pinned local
Git object. The new helper follows this weighting, with the existing PeleF Clausius
or bounded Antoine inversion at local pressure. It does not claim physical mixture
bubble-point or boiling-rate parity. Grid/CLI/mixed-film restart remains open.

## Focused verification

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks enabled:

```sh
cmake --build build/port-spray-debug --target test_spray_core --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(mixture_wall_crossing|mixture_wall_event|wall_event)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Three entries passed in 0.03 s total, no failures/skips. New checks cover reference
pressure weighted boiling temperature, wall crossing and deposition composition,
wall location/normal, no-crossing preservation and invalid-pressure rollback.
Existing mixed/pure event cases retain six-regime conservation and failure coverage.
Provenance, contract and whitespace checks passed; no Python source changed. No CLI,
full configuration matrix or physical validation was run.
