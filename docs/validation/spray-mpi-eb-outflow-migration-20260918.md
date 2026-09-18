# MPI EB outflow with final ownership migration — 2026-09-18

Implementation `4b9c07f`, base `6035e3d`.

`advance_mpi_eb_tree_impacts` now accepts optional `outflow_change(5+nc)`.
It collectively validates output size, forwards private storage to the MPI
impact interval, then publishes the globally reduced exterior inventory only
after final particle and film migration succeed. The existing interval enforces
consistent optional presence across ranks. Failed migration preserves caller
particles, component masses, films, ID watermark, hit count and both ledgers.

## Completed checks

GNU Fortran 13.3.0/Open MPI 4.1.6, Debug/Ninja, reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Selected library compiled only.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_impact --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_impact_(2|4)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **2/2 passed**, 0 failures/skips, 0.98 s total with configured
`/usr/bin/mpiexec`. Builds, provenance, syntax, contract and diff checks passed.
Each rank starts with one parcel; its four splash children escape. Final
population is empty, watermark and global hit count both equal five times the
rank count. Films concentrate on the finest fluid-side owner (last rank).
Global film plus wall receipt minus initial inventory equals exterior change
within existing relative 3e-12 tolerance and 1e-20 floor. Component loss equals
mass loss for this single-water fixture.

Setting per-rank film capacity to one rejects final concentration after the
private splash/outflow work, preserving original local population, positions,
component masses, empty films, watermark and budget sentinels. Capacity 100
succeeds. Existing impact, no-outflow migration, repeated-contact and mismatched
outflow-mode checks remain in the same executable.

Artifacts: `build/spray-mpi-eb-outflow-migration-checks/{provenance.log,build.log,ctest.log,build-selected.log}`.
This covers the interval and migration wrapper. MPI source/reactive and
application wiring remain pending; no physical or external-solver parity claim.
