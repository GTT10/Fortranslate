# Selected-level SDC redistribution output — 2026-09-27

Implementation `06f8c664`, based on `9ad29599` on local `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

Serial and MPI selected-level SDC now accept `redistribution='flux_redist'`
when the caller supplies `redistribution_transfer` along with the existing
required integrated-flux output. The selected level advances using its local
part of each redistributed rate. Other gas levels retain their state while
the corresponding time-integrated changes are returned separately.

The accepted transfer is the final hydro transfer plus half the sum of the
old and final-new transport transfers. Intermediate SDC iterations replace
provisional values rather than accumulate them. Disabled transport contributes
zero. All outputs remain private until the stage succeeds. The MPI interface
collectively agrees transfer-output presence along with the existing controls.
Requesting transfer output outside selected-level flux redistribution is rejected.

Common-step routing remains unchanged. Recursive drivers do not yet accumulate
or apply the returned transfers; recursive application redistribution remains
disabled. This is a stage interface, not completed recursive integration.

## Completed focused checks

GNU Fortran 13.3.0, `/usr/bin/mpifort`, `/usr/bin/mpiexec`
(OpenMPI/OpenRTE 4.1.6), Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance and the touched MPI full-step test target built successfully.
- `unit_mpi_eb_tree_sdc_redist_level_3d_2`: 1/1 passed,
  7.65 s test / 7.70 s total.
- Three levels, five patches, selected depth 1, all transport operators,
  reactions and external force, two SDC iterations, dt 1e-10 s.
- The serial state change plus deferred transfers balances the divergence of
  exported face fluxes and the wall, reaction and forcing inventories. The
  component-wise discrepancy divided by max(1, abs(flux inventory)) is
  3.0799545883542763e-13 (test bound 2e-9). Coarse-level transfer is nonzero;
  same-level transfer is exactly zero and other gas levels remain unchanged.
- For round-robin owners and all patches on rank zero, MPI state, temperature,
  reaction estimates and transfer fields match serial exactly. Transfers allocate
  only owned fields; ledgers agree within 2e-12 relative and limiter is exact.
- A hydro-stage rejection after initial transport evaluation preserves gas,
  temperature, reaction estimates, transfer fields, ledgers and limiter.
- Python syntax, project contract and `git diff --check` passed. No failed,
  skipped or interrupted checks were counted.

Artifacts: `build/sdc-eb-redist-level-loop-checks/`, with `provenance.log`,
`build.log`, `ctest.log`, `detail.log` and `contract.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug \
  --target test_mpi_eb_tree_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 \
  OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug \
  -R '^unit_mpi_eb_tree_sdc_redist_level_3d_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

No application trajectory, selected-mechanism build, installed audit, full matrix,
recursive transfer schedule, temporal convergence or physical validation was run.
