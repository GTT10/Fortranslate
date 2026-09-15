# Mixed-liquid KH children — 2026-09-16

Tested parent: `128a533` plus this change. Source/test SHA-256 hashes and build,
test and preflight logs are retained in `build/mixture-kh-checks/`.

## Implementation

The regular-grid mixture route now uses the existing KH shed-volume accumulator
and four-child constructor. Each child inherits the parent component fractions;
parent component masses scale with its remaining per-drop mass. Transverse kinetic
energy is debited from the gas through the existing interpolation weights. Children
advance through the remaining interval, and boundary losses retain component budgets.
Capacity, invalid boundary state or inadmissible gas energy rejects the whole update.
Computational splitting partitions KH inventories and copies per-drop composition.

The CLI accepts mixed-liquid KH children and retains both KH and component state
through the existing combined checkpoint fields. No checkpoint layout or evolution
tag changes: this option combination was previously rejected. Mixture wall impact
remains unsupported; viscosity and surface tension are prescribed mixture properties.
The four symmetric children are a PeleF choice, not upstream stochastic parity.

## Focused verification

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks enabled:

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d test_spray_integration --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_integration_mixture_grid$' --output-on-failure
python3 tests/spray/test_spray_application.py --exe build/port-spray-debug/pelef_spray_3d SprayApplication.test_mixture_kh_children_restart
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```

One integration entry passed (0.04 s); one CLI method with release and accumulation
subcases passed (2.060 s), without failures or skips. Integration checks cover
component/total inventory conservation and capacity-failure rollback. CLI checks
cover inherited fractions, transverse motion, bounded computational multiplicity,
retained accumulation, and exact restarted checkpoint, gas/parcel CSV and final
history agreement. Provenance, syntax, contract and whitespace checks passed.
The full CLI suite, configuration matrix and physical validation were not rerun.
