# SDC outflow application routing — 2026-09-26

Implementation `854b6a11`, parent `5a5d2202`, local main, version 0.246.0.

Fixed and selected regular gas applications now forward resolved boundary
faces to SDC and accumulate its boundary inventory. Input accepts periodic
and zero-gradient outflow faces with PPM/PeleC, including reactions and
external forcing. Physical-boundary transport and prescribed outlet pressure
are explicitly rejected. Inflow, walls and AMR/MPI/EB SDC remain open.
Example: `cases/reactive_boundaries_3d/sdc_outflow.nml`.

## Completed checks

- GNU Fortran 13.3.0 Debug fixed tests-enabled and selected tests-disabled
  builds passed in `build/spray-amr-reactive-debug` and
  `build/spray-amr-reactive-selected-debug`, respectively. Pinned H2/O2 bundle;
  existing compiler warnings remain.
- **2 application cases passed**, fixed/selected with matching compositions:
  a 4³ reacting hotspot, all-outflow boundaries, external force, no molecular
  transport, two SDC iterations and one step to 1e-9 s. Temperature spans
  **1108.555019–1271.822881 K**. Maximum normalized Euler inventory residual
  is **3.61349e-18** after boundary/source accounting.
- All fixed/selected CSV entries match exactly (maximum scaled difference 0).
  This checks shared-driver routing, not an independent solver comparison.
- **4 expected input rejections passed**: physical-boundary transport and
  prescribed outlet pressure for both executables. No unexpected failures.
- Provenance, Python syntax, contract and whitespace passed. No CTest rerun
  for this dispatch change; library evidence is in `sdc-outflow-20260926.md`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_3d_selected --parallel 2
python3 build/sdc-outflow-app-checks/check.py
python3 -m py_compile tools/check_project_contract.py build/sdc-outflow-app-checks/check.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/sdc-outflow-app-checks/` contains checker, inputs, CSVs,
build/run/rejection logs, provenance and report. One short nonuniform step
does not establish temporal convergence, ignition accuracy or physical
boundary-flow validation.
