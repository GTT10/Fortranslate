# Dynamic AMR LES inflow and outflow — 2026-09-25

Implementation `4cdd1cb`, based on `733c332`. Serial/MPI AMR timestep selection,
SSPRK2, reacting-spray coupling and the application now accept prescribed
reservoir inflow and extrapolated outflow through the existing wide donor halos.
Walls and nonzero pressure targets remain rejected for Dynamic LES.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug. Fixed tests-enabled targets and
the selected H2O2 tests-disabled application built successfully. The selected
application was not executed in this check.

Focused CTest: **2/2 passed**, zero failures, no skips, total 96.52 seconds:

- `spray_mpi_amr_dynamic_open`: 14.98 seconds. Gaussian-filter Dynamic LES with
  reservoir inflow, extrapolated outflow, scheduled regridding and 2-to-4-rank
  spray restart. Existing gas/particle continuation and invalid-checkpoint
  checks pass.
- `spray_les_mpi_amr_tree_2`: 81.50 seconds. Periodic, all-outflow and mixed
  inflow/outflow Dynamic cases on root-only and nested trees; serial/MPI state,
  timestep and diffusivity agreement, composite boundary-flux conservation,
  nonzero open-boundary inventory transfer and rollback. Both concentrated
  ownership (including an empty rank) and distributed ownership are exercised.
  Existing static LES boundary cases remain in the same test.

No numerical tolerances changed. Provenance preflight, Python syntax, project
contract and whitespace checks passed. The standalone serial test target was
built; serial operators were exercised by the MPI comparison test.
This is focused implementation evidence, not external PeleC parity, physical
validation, scaling or a full configuration matrix. Wall, pressure-target and
EB Dynamic LES remain outside this change.

## Commands and artifacts

Logs: `build/les-dynamic-open-checks/` (`provenance.log`, `build.log`,
`configure.log`, `build-selected.log`, `ctest.log`, `contract.log`).
Application artifacts:
`build/spray-amr-reactive-debug/tests/spray_mpi_amr_dynamic_open/restart_w2kg3g5_/`.
MPI launches use the host workaround recorded in
[the launcher recovery](mpi-dynamic-amr-runtime-20260925.md).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_amr_tree test_les_mpi_amr_tree pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_les_mpi_amr_tree_2|spray_mpi_amr_dynamic_open)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py tests/application/check_mpi_amr_spray_restart.py
python3 tools/check_project_contract.py --source .
git diff --check
```
