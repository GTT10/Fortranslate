# Per-port breakup properties — 2026-09-16

Source: parent `f475ea0` plus this change; GNU Fortran 13.3 Debug,
`build/port-spray-debug`, fixed H2/O2 mechanism, MPI disabled.

## Implemented behavior

Ports can override surface tension, dynamic viscosity and the existing
inverse-temperature viscosity fit. The application resolves inheritance and
validates fit bounds before injection. Both coupling half steps select the
breakup options through parcel material identity. Global breakup model,
calibration constants, integration controls and splitting limits remain shared
in the application. Resolved properties enter checkpoint identity. Core ETAB
and KH-RT equations, conservation transfer and rollback rules are unchanged.

## Focused evidence

Six CTest entries passed, zero failures/skips, in 15.06 s: viscosity, KH-RT,
ETAB, integration rollback, integration materials and application CLI.
The CLI contains 18 tests (CTest wall time 14.77 s). Its added case exercises
ETAB and KH-RT separately with identical initial droplets but different port
surface tensions and constant/inverse-temperature viscosities. It checks that
actual diameter/deformation states differ, all-phase balances hold, and gas
and parcel CSVs plus final history exactly match after restart. Changed surface
tension rejects restart; zero surface tension, injection outside the viscosity
range and unused port viscosity inputs reject.

An initial six-entry run also passed; after explicitly retaining global minimum
diameter and breakup substep limits in the core API, the same six entries were
rerun with the final results above.

Build, mechanism provenance, Python syntax, project contract and whitespace
checks passed. No failures or broad configuration matrix were encountered.
These are implementation/conservation checks, not calibrated breakup validation.

Artifacts: `build/port-spray-debug/port-breakup-build.log` and
`build/port-spray-debug/Testing/Temporary/LastTest.log`. CLI fixtures use temporary
directories. Mechanism source SHA-256 remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d test_spray_integration --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(application_cli|integration_(materials|rollback)|etab|khrt|viscosity)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```
