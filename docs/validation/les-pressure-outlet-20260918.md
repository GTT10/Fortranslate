# Serial LES pressure outlet — 2026-09-18

Implementation `af5b6ddd4b209e5f64573c8ad35f6946e9631a83`, based on `57112a3`.

The serial LES RHS and SSPRK2 APIs accept optional trailing `outlet_pressure(6)`.
They validate finite values in `[0,1e15]` Pa and allow positive values only on
outflow faces. Omission/zero retains ordinary boundary behavior. After the
existing viscous ghost fill, `fill_pressure_outlet_ghosts_3d` supplies the two
characteristic pressure layers used by the shared face-gradient kernel. Both
SSPRK2 stages follow this path, with existing interval checks and transactional
EOS recovery. This connects an existing pressure boundary algorithm; it does
not introduce a new characteristic formulation. MPI and coupled/application
selection remain open. The ghosted RHS API still expects prepared halos.

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing configurations unchanged.
Fixed tests-enabled `test_spray_core` and selected H2O2 tests-disabled
`pelef_selected_spray_runtime` built. Selected runtime was not run.
Focused CTest **3/3 passed**, 0.43 seconds, zero failures/skips:
`spray_les_boundary`, `spray_les_pressure_outlet`, `spray_les_isotropic`.

The new case covers all three static closures with CI 0.1. It checks exact
zero-pressure RHS/diffusivity compatibility, nonzero flux response to a
100000 Pa outlet from a 101325 Pa fixture, RHS and integrated boundary ledgers,
and a successful nontrivial SSPRK2 update. Ledger tolerance is `1e-12` with unit
floor. Nonzero pressure on a periodic face, NaN pressure even at zero dt, and
negative pressure are rejected; rejected steps preserve state/temperature.
No tolerances changed. Provenance, Python syntax, project contract and whitespace
checks passed. No builds/tests failed. No MPI execution, application run, full
matrix, external solver parity or physical validation was performed.

## Commands and artifacts

Logs: `build/les-pressure-outlet-checks/{provenance,build,build-selected,ctest}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_(pressure_outlet|boundary|isotropic)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
