# PCM tree boundary checkpoint repair — 2026-09-26

Implementation `ff469d3`, parent `411b04b`, local main, version 0.246.0.
The prior explicit-wall KH fixture advanced with PCM but failed to save because
tree patch records reused a single-level configuration check that requires CTU
for physical boundaries. Serial/MPI tree readers and writers now pass an
explicit `tree_patch` context to their embedded records. That context permits
PCM boundary configuration without changing the stored identity or file schema.
It requires a borrowed checkpoint unit. Default single-level calls remain strict.
Finite inlet/thermal data, paired/known faces, enabled thermal transport, solver
identity, geometry, state, composition and diagnostics still undergo validation.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3, Ninja Debug, fixed H2/O2,
existing tests-enabled `build/spray-amr-reactive-debug`.

The previously failing explicit no-slip-wall PCM/HLLC KH fixture now passes
**3 executions and 1 expected context rejection**: MPI4 full, MPI2 checkpoint,
MPI4 continuation. It includes three regrid stages, measured-rebalance control,
particle reflection and KH children. Gas, particle and budget continuation agree
within the unchanged 3e-12 bound, with lineage and pending KH state retained.
The PCM checker mode and CTest registration preserve this reproducer.

Existing single-level configuration and checkpoint unit tests pass **2/2**,
0.52 s and 0.15 s (total 0.72 s). These include default rejection of PCM physical
boundaries and the existing malformed-record/rollback checks. No unexpected
failures or skips occurred in this repair turn. The original failure is retained
under `build/mpi-eb-wall-khrt-checks/application-initial.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_eb_3d_config test_reactive_eb_3d_checkpoint --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tests/application/check_mpi_eb_impact.py --exe build/spray-amr-reactive-debug/pelef_mpi_eb_spray_3d --mpiexec /usr/bin/mpiexec --source . --work build/eb-tree-pcm-checkpoint-checks/application --khrt --domain-wall --pcm
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_eb_3d_(config|checkpoint)$' --output-on-failure
python3 -m py_compile tests/application/check_mpi_eb_impact.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, project contract and whitespace passed. Build warnings remain.
No numerical update, tolerance, golden data or checkpoint schema changed.
Artifacts under `build/eb-tree-pcm-checkpoint-checks`: `build.log`,
`build-unit.log`, `application.log`, `ctest.log`, `LastTest.log`,
`application/impact_x_o5i247`. The application checker ran directly; its new
CTest registration was not counted as another execution.
Mechanism SHA: `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
PeleC: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
PelePhysics: `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.

No selected-mechanism, full matrix, installed, external-solver or physical
qualification ran. Serial tree forwarding was built; a separate PCM serial-tree
continuation was not rerun. Complete PeleC parity remains unproven.
