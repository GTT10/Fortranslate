# Combined EB/exterior particle drift — 2026-09-18

Implementation `e87b6ab`, base `c72f581`.

The optional exterior-boundary group extends specular EB drift to root-domain
outflow/inflow removal, specular wall reflection and periodic wrapping. Each
segment is clipped at its first exterior crossing before the original level-set
intersection is queried; an earlier EB hit takes precedence. Remaining flight
continues after each event. The existing interface retains its root-contained
behavior and ledger calculation when exterior options are absent.

The caller supplies per-parcel specific liquid internal energy. EB and exterior
ledgers are separate; escaping inventory includes internal plus kinetic energy.
Escaped parcels retain records with zero mass. Component inventory retirement
and species budgets belong to the future coupled-source connection. A common
event limit and private trial population preserve complete rollback. Periodic
faces must be paired and their nodal geometry must match exactly. Noncrossing
directions skip crossing-time division, avoiding extreme small-step ratios.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Executed test was serial.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_drift --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_drift$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **1/1 passed**, total **0.04 s**. Existing oblique reflection,
thin-cavity repeated hits, zero interval and rollback checks remain. Added cases
verify an EB hit followed by exterior escape, an EB hit followed by exterior
reflection, three periodic wraps, common hit-budget failure, mismatched periodic
geometry rejection and a later invalid start rolling back an earlier escape.

The escape case has represented mass 0.03 kg, unit speed and specific internal
energy 10 J/kg: EB x-impulse is +0.06 kg m/s; outflow removes 0.03 kg,
0.03 kg m/s and 0.315 J. Exterior reflection instead cancels the EB impulse.
Numerical assertions use absolute 1e-14 where applicable. Periodic drift has
zero inventory change. Failed calls retain parcels, both ledgers and hit count.

Fixed/selected library builds, provenance, Python syntax, contract and whitespace
checks passed. No failures or skips. Selected runtime, source/application
coupling, non-specular EB impact/film flight, MPI execution, broad matrix,
installed audit, external parity and physical validation were not run.

Artifacts: `build/spray-eb-outer-drift-checks/`, including `provenance.log`,
`build-final.log`, `ctest-final.log`, `build-selected-final.log`.
