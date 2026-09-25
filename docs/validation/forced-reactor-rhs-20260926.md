# Forced reactor right-hand side — 2026-09-26

Implementation `ef157d79`, parent `3ac5b01a`, local main, version 0.246.0.

`forced_reactor_rhs` evaluates chemistry and fixed external species-density
rates using density = sum(rhoY) and internal-energy density
rhoe(t) = rhoe(0) + t * external_energy. It recovers temperature from the
resulting specific internal energy and composition. An empty reaction list
is explicitly inert. Invalid inputs leave outputs zero. This is a library
building block; adaptive forced integration and SDC reaction coupling remain open.

## Completed checks

- GNU Fortran 13.3.0 Debug build, existing fixed tests-enabled directory
  `build/spray-amr-reactive-debug`, target `test_implicit_h2o2_reactor` passed.
- `unit_implicit_h2o2_reactor`: **1/1 passed**, 0.01 s test / 0.05 s total.
  Added checks cover zero forcing versus the existing chemistry RHS,
  simultaneous density/energy addition in the inert limit, chemistry plus
  external species rates, and rejection of negative species density.
  The existing implicit constant-volume integration check remains active.
- Initial test failed because the underlying kinetics routine rejects an
  empty mechanism. The new forced RHS explicitly handles that inert limit;
  no tolerance changes were made. No final failures or skips.
- Mechanism provenance, Python syntax, project contract and whitespace passed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_implicit_h2o2_reactor --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_implicit_h2o2_reactor$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/forced-reactor-rhs-checks/` contains provenance/build logs,
`initial-failure.log` and the successful `LastTest.log`. The original
`ctest.log` also records the initial failure. No application or external solver
comparison was run for this change.
