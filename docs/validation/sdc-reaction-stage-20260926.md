# Independent SDC reaction stage — 2026-09-26

Implementation `a2600bd1`, parent `320d9c71`, local main, version 0.246.0.

`react_sdc_state_3d` exposes the existing cell-local SDC reaction update as a
separate stage. Inputs are interval-start state/temperature, nonreacting rates
(hydro plus averaged sources), timestep and chemistry tolerances. It performs
the existing forced full-interval integration, updates density and total
energy, and returns reaction-source estimates. Shape, finiteness, timestep,
old-state thermodynamics and final-state recovery are checked. All outputs
are published only on complete success. Ordinary SDC now calls this routine;
the final temperature is recovered once, retaining the existing ordering.

This lets a future hierarchy coordinator evaluate all hydro candidates and
new-time transport/source rates before invoking reaction updates across
patches. It does not itself provide a hierarchy stage barrier or interpolation.

## Completed checks

- GNU Fortran 13.3.0 Debug fixed tests-enabled build
  `build/spray-amr-reactive-debug`, target `test_reactive_ppm_step_3d`, passed.
- `unit_reactive_ppm_step_3d`: **1/1 passed**, 22.01 s test / 22.06 s total.
  Added direct one-cell stage comparison against forced full-interval
  integration (species within 1e-12 absolute, temperature within 1e-6 K).
  Invalid timestep preserves state, temperature and reaction-source outputs.
  Existing complete reacting SDC, callback and flux-export checks remain.
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

Artifacts: `build/sdc-reaction-stage-checks/` contains build, provenance,
CTest and LastTest logs. These focused checks do not establish complete
AMR/MPI SDC, external solver parity or physical validation.
