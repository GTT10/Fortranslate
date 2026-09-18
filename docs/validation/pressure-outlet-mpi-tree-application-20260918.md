# MPI patch-tree pressure outlet application — 2026-09-18

Implementation `4e5a195`, base `108c251`.

Fixed/selected MPI tree frontends enable pressure-outlet input and the shared
driver forwards targets to coupled common or recursive stepping. Existing
canonical configuration serialization includes the targets in checkpoint
identity; no checkpoint format change is required. MPI operator entry points
check target agreement collectively.

## Focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6, `/usr/bin/mpiexec`; reused Debug/Ninja fixed
tests-enabled and selected H2O2 tests-disabled configurations.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_tree_reactive_3d test_selected_elementary_mpi_amr_tree_3d pelef_amr_tree_reactive_3d test_selected_elementary_amr_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_tree_reactive_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^regression_mpi_amr_tree_pressure_restart_3d$' --output-on-failure
python3 -m py_compile tools/check_mpi_amr_tree_restart_3d.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The scheduled hotspot fixture adds x-high 102325 Pa, moving no-slip y walls,
active chemistry and molecular transport. Common steps traverse three topology
changes. Fixed and generated selected elementary mechanisms each compare state
and budget CSV bytes for uninterrupted MPI/serial execution, 2-to-3-rank restart,
MPI-to-serial and serial-to-1-rank continuation, and root-only checkpoint reload
with empty ranks. Serial/MPI checkpoint bytes are also compared exactly.
Changed/removed pressure, changed physics/schedule, truncation, inconsistent
ledger, malformed input and checkpoint path aliasing are rejection cases.

Recursive pressure splitting was checked separately at library level; this
application regression uses common timesteps. No physical-validation,
nonreflection, scaling or installed-binary claim is made. Tolerances are unchanged.

Build/check logs: `build/pressure-outlet-mpi-tree-app-checks/`.
Runtime fixtures/logs:
`build/spray-amr-reactive-debug/tests/mpi_amr_tree_pressure_restart_3d/{fixed,selected}/`.

Final CTest **1/1 passed**, zero skips, 55.07 s total (test 55.03 s).
Fixed/selected builds and their MPI ABI-40 dependency audits passed, as did
provenance, Python syntax, project contract and diff checks. No failed or
interrupted runs.
