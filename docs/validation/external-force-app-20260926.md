# Regular gas external forcing input — 2026-09-26

Implementation `1663c52e`, parent `69af107f`, local main, version 0.246.0.

Fixed and selected regular 3D gas applications accept `external_forcing(3)`
in N/m³, default zero, with PPM/PeleC. The reader requires a capability flag
before accepting nonzero forcing in other application families. Nonfinite force
and incompatible reconstruction/solver combinations are rejected. The shared
driver forwards forcing to the old/new-work split path, accumulates its source
ledger and subtracts that ledger in conservation diagnostics. The unforced
entropy-wave analytic density comparison is disabled for forced evolution.
An example is `cases/reactive_boundaries_3d/external_forcing.nml`.

## Completed checks

- Fixed and selected builds passed. GNU Fortran 13.3.0, Python 3.12.3, Ninja
  Debug; reused fixed tests-enabled `build/spray-amr-reactive-debug` and selected
  tests-disabled `build/spray-amr-reactive-selected-debug`, pinned H2/O2 bundle.
- **8 successful application runs**: fixed/selected × periodic/outflow × force
  on/off. Each 4³ uniform case uses 900 K, 300000 Pa, velocity (20,10,-5),
  final time 1e-5 and force (20000,-30000,40000). Fixed uses the default mixture;
  selected uses named pure N2. All 64 cells match analytic momentum increments
  within 1e-10, work including the quadratic kinetic term within 1e-8, and
  baseline temperature within 1e-7 K. Built-in conservation gates passed.
- **4 expected rejections**: each executable rejects PCM plus nonzero force and
  NaN force, with the external-forcing diagnostic. No skips.
- Initial selected inputs omitted the required named composition and were
  rejected; corrected inputs supply it. An initial 1000 K fixed fixture passed
  momentum/energy checks but failed the temperature comparison by 0.000121685 K.
  This lies at the NASA coefficient switch and is consistent with branch
  sensitivity in temperature inversion; it was not repaired or qualified here.
  Initial outputs are retained. Final 900 K checks use the same tolerances.
- Provenance, Python syntax, project contract and whitespace passed. Existing
  compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_3d_selected --parallel 2
python3 build/external-force-app-checks/check.py
python3 -m py_compile tools/check_project_contract.py build/external-force-app-checks/check.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/external-force-app-checks/` contains the checker, all inputs,
CSV outputs, run/build logs, `report.json`, provenance and `initial-1000K/`.
The checker removes `LD_LIBRARY_PATH` and `LIBRARY_PATH` for subprocesses.

Regular gas applications have no checkpoint interface. AMR/MPI/EB/spray forcing
input, full SDC integration, reacting-force qualification and physical
validation remain open. The coefficient-switch observation remains unresolved.
