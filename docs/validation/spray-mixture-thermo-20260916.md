# Mixed-liquid thermodynamics — 2026-09-16

## Implementation and provenance

Tested parent `06d13bd` plus this change; exact source/test hashes and logs are
retained in `build/mixture-thermo-checks/`. `references/spray_mixture.json` pins
PelePhysics SprayPropertiesPeleMP.H at
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`; its SHA256 was checked from the local
pinned Git blob.

The new module evaluates additive-volume density and mass-weighted heat capacity
and stored energy for up to 16 pure-component records. It recovers temperature
from energy inside the active fit intersection and applies fixed-composition
parcel heating transactionally, returning the gas energy debit. Component energy
references and polynomial cp integrals reuse the existing liquid implementation.

The preceding wet-impact limitation was also checked against the pinned
AhamedSplash.H: its thick-film forward mass polynomial has the same positive
1.9e-4 quadratic coefficient as the Fortran port. The invalid nearly-normal
mass partition therefore remains a reference-correlation limitation; no coefficient
or rejection guard was changed in this increment.

## Focused verification

GNU Fortran 13.3.0, existing Debug spray build, MPI disabled, bounds/FPE enabled:

```sh
cmake --build build/port-spray-debug --target test_spray_core --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_mixture$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

The new test passed (1/1, 0.03 s total), without failures or skips. Provenance,
project contract and whitespace checks passed. No Python source changed.
It checks independently calculated two-component density/cp/energy, inversion
from a distant guess, nonlinear polynomial-cp heating, liquid/gas energy closure,
fixed mass, unreachable-state rollback, composition rejection and a pure-component
limit whose inactive second component has a narrower fit range.
The H2O/H2O2 values are synthetic implementation fixtures, not a property dataset
or experimental mixture validation. CMake automatically regenerated the existing
build to include the new module and test registration.

## Remaining work

This is a thermodynamic API, not a complete mixed-droplet solver. Component mass
storage, vapor partitioning, composition-changing evaporation, CFD/CLI integration
and mixture restart remain open. Boiling-temperature property caps, excess mixing
energy and physical validation are not implemented here. No full configuration
matrix or external solver comparison was run.
