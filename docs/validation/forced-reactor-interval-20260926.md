# Forced reactor full interval — 2026-09-26

Implementation `5915aa99`, parent `6cbe939f`, local main, version 0.246.0.

`advance_forced_reactor_interval` repeats adaptive implicit substeps until the
requested duration is covered. It carries elapsed time against the original
internal-energy density, limits the final substep to the remaining interval,
and publishes species densities and temperature only on full success.
An optional initial step and maximum substep count control integration;
the default maximum is 100000. Failure returns the number of accepted local
substeps for diagnosis while preserving caller state. A valid zero interval
is a no-op. Stalled floating-point time advancement is rejected.

## Completed checks

- GNU Fortran 13.3.0 Debug, fixed tests-enabled build
  `build/spray-amr-reactive-debug`, target `test_implicit_h2o2_reactor` passed.
- `unit_implicit_h2o2_reactor`: **1/1 passed**, 0.09 s test / 0.14 s total.
  Added checks cover rollback after one accepted substep when the configured
  cap is reached, four-step inert integration against the analytic species
  and temperature result, and reacting full-interval integration with mass
  and internal-energy inventory checks. Existing substep/RHS checks remain.
- Mechanism provenance, Python syntax, project contract and whitespace passed.
  No failed or skipped tests; existing compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_implicit_h2o2_reactor --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_implicit_h2o2_reactor$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/forced-reactor-interval-checks/` contains build, provenance,
CTest and LastTest logs. SDC cell/application coupling and comparison with an
external reacting solver remain open. This focused library check does not
establish physical validation or complete PeleC equivalence.
