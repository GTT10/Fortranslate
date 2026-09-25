# Caller-coordinated SDC source estimates — 2026-09-26

Implementation `1af831e2`, parent `38f4f2db`, local main, version 0.246.0.

The SDC library accepts optional `new_source_estimate` with shape
(nvar,nx,ny,nz,source_count). It initializes the new-time estimate from this
array and returns the final evaluated rates only after success. Shape and
finiteness are checked before work. Existing calls without the array retain
old-source initialization.

For one-iteration calls, the coordinator initially seeds this array with old
source rates, retains it and any reaction-source array, and resets state and
temperature to the same interval-start values before each call. Time is not
advanced between these calls. This is a prerequisite for hierarchy control;
coherent intermediate neighbor data still require stage coordination.

## Completed checks

- GNU Fortran 13.3.0 Debug fixed tests-enabled build
  `build/spray-amr-reactive-debug`, target `test_reactive_ppm_step_3d`, passed.
- `unit_reactive_ppm_step_3d`: **1/1 passed**, 21.82 s test / 21.86 s total.
  Added three externally coordinated one-iteration calls with state-dependent
  source families. They match the existing three-iteration call within 1e-14
  scaled state error. A failed new-source evaluation preserves both caller
  state and the source estimate. Existing SDC tests remain active.
- Provenance, Python syntax, project contract and whitespace passed. No failed
  or skipped tests; existing compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_step_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_step_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/sdc-iteration-state-checks/` contains build, provenance,
CTest and LastTest logs. No hierarchy driver, coarse/fine interpolation,
reflux or MPI exchange is implemented by this change; those remain open.
