# Dynamic LES adaptation and selected continuation — 2026-09-25

Test registration `ebc42d02613fee771ddae86e4bd98d2478e8beee`, based on `b90f1e1`.
The implementation is unchanged; this closes two specific unverified combinations.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing fixed tests-enabled and
selected H2O2 tests-disabled configurations. Both application targets built.
MPI launches used `HWLOC_COMPONENTS=-x86` and
`OMPI_MCA_hwloc_base_binding_policy=none`, following the recorded launcher recovery.

- CTest `spray_mpi_amr_dynamic_adaptation`: **1/1 passed**, 26.74 seconds (total
  26.78). Two-port injection, temperature-triggered refinement and measured-work
  rebalancing run with Dynamic Gaussian LES. Continuous and 2-to-4-rank restart
  gas/particle results agree. Output contains both levels and more than two
  patches; logs verify adaptation resumes at the correct step and rebalancing
  executes after continuation. Changed filter context and malformed/immutable
  checkpoint cases reject as expected.
- Direct selected H2O2 checker: **passed, exit 0**. Dynamic LES, scheduled
  regridding and 2-to-4-rank continuation use the selected-mechanism executable
  and named composition. Context/truncation/overwrite rejection checks pass.
  This selected run does not exercise adaptive regridding/rebalancing.

Existing comparison tolerance is relative `2e-12`, absolute `1e-25`; no tolerance
or solver implementation changed. No unexpected failures/skips occurred.
Provenance preflight, Python syntax, project contract and whitespace pass.
These are bounded implementation combinations, not detailed-fuel Dynamic LES,
physical validation, scaling or a full configuration matrix.

## Commands and artifacts

Logs: `build/les-dynamic-adaptive-checks/` (`provenance.log`, `build.log`,
`build-selected.log`, `ctest.log`, `selected.log`, `contract.log`).
Fixed artifacts:
`build/spray-amr-reactive-debug/tests/spray_mpi_amr_dynamic_adaptation/restart_jswyjxlp/`.
Selected artifacts:
`build/les-dynamic-adaptive-checks/selected/restart_qwfqh2ly/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_dynamic_adaptation$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tests/application/check_mpi_amr_spray_restart.py --exe build/spray-amr-reactive-selected-debug/pelef_mpi_amr_spray_3d_selected --mpiexec /usr/bin/mpiexec --source . --work build/les-dynamic-adaptive-checks/selected --selected --les --dynamic
python3 -m py_compile tools/check_project_contract.py tests/application/check_mpi_amr_spray_restart.py
python3 tools/check_project_contract.py --source .
git diff --check
```
