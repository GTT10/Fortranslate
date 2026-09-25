# SDC full-face flux export — 2026-09-26

Implementation `cd8ce348`, parent `08e786de`, local main, version 0.246.0.

The SDC library now optionally returns all x/y/z face fluxes, including lower
domain faces. Outputs contain the final iteration's hydro flux plus the
average of old and final-new built-in transport fluxes. Periodic lower faces
receive the corresponding wrapped transport flux. Output is flux density;
callers multiply by timestep and face area for a flux register. External and
reaction volume sources are not included; custom source callbacks do not
provide face fluxes. All three output arrays must be supplied with full-face
shapes. Outputs remain zero on failure and state publication is transactional.

## Completed checks

- GNU Fortran 13.3.0 Debug, fixed tests-enabled build
  `build/spray-amr-reactive-debug`, target `test_reactive_ppm_step_3d`, passed.
- `unit_reactive_ppm_step_3d`: **1/1 passed**, 21.52 s test / 21.56 s total.
  Added nonuniform periodic and mixed outflow/periodic cases with molecular
  transport and no volume sources. The returned full-face divergence
  reconstructs every conserved component in every cell within 1e-12 scaled
  error. A malformed face shape rejects with unchanged state/temperature
  and zero provided outputs. Existing reacting/source/boundary checks remain.
- Provenance, Python syntax, contract and whitespace passed. No failed or
  skipped tests; existing compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_step_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_step_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/sdc-flux-export-checks/` contains build, CTest, LastTest and
provenance logs. This provides fluxes needed by AMR conservation correction;
it does not yet connect AMR halo interpolation, reflux, MPI or EB SDC, nor
establish complete solver equivalence or physical validation.
