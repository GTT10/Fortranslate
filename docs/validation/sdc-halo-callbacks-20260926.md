# Caller-supplied SDC halos — 2026-09-26

Implementation `45013b16`, parent `9df7f815`, local main, version 0.246.0.

SDC accepts paired optional state and source halo callbacks. The state callback
receives time, halo width and a viscous flag: old-time hydro uses four layers,
while old/new transport evaluations use two viscous layers. The callback owns
all ghost treatment, including pressure conditions when applicable. Halo CTU
retains the supplied face masks; ghosted physical transport uses the supplied
state halos. The source callback fills four layers for the predictor rates at
the midpoint time. Callers must populate the full arrays including interiors.
Absent callbacks retain existing local halo behavior.

The paired interface allows a hierarchy or MPI driver to supply neighbors
and coarse interpolation without silently wrapping a patch's own data.
It does not implement hierarchy stage synchronization itself; callers must
provide coherent neighboring/coarse data for each evaluation. Built-in
transport continues to export full-face fluxes and boundary inventories.

## Completed checks

- GNU Fortran 13.3.0 Debug fixed tests-enabled build
  `build/spray-amr-reactive-debug`, target `test_reactive_ppm_step_3d`, passed.
- `unit_reactive_ppm_step_3d`: **1/1 passed**, 21.70 s test / 21.74 s total.
  Added periodic callback test with nonuniform state and molecular transport.
  Two SDC iterations produce four state-halo calls (one hydro, three transport)
  and two source-halo calls. Result matches native periodic state within
  1e-12 scaled error and temperature within 1e-7 K. A failed new-time halo
  evaluation preserves the complete input state and temperature.
- Provenance, Python syntax, project contract and whitespace passed. No
  failed/skipped tests; existing compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_step_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_step_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/sdc-halo-callback-checks/` contains build, provenance, CTest
and LastTest logs. Actual coarse/fine interpolation, hierarchy iteration,
reflux integration and MPI exchange remain open. This focused callback
check is not AMR/MPI qualification or physical validation.
