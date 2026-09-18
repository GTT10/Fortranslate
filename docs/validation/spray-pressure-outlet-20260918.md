# Serial spray pressure outlet — 2026-09-18

Implementation `552a54d538ade7fd42f7c1893b5d7be6a1ee0f9d`, based on `21089b8`.

Serial spray coupling and fixed/selected serial application namelists accept
`outlet_pressure(6)` in Pa, default zero. Finite values in `[0,1e15]` are allowed,
with positive values restricted to outflow faces. The coupling forwards pressure
to hydro/molecular transport and both static LES half steps. Particle outflow
semantics are unchanged. Nonzero settings enter checkpoint identity, while zero
preserves the previous context. MPI coupled/application controls remain pending.

## Completed verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing configurations unchanged.
Fixed tests-enabled `pelef_spray_3d` and selected H2O2 tests-disabled
`pelef_spray_3d_selected` built. Selected executable was not run.
Focused CTest **2/2 passed**, 19.34 seconds, zero test failures/skips:
`spray_isotropic_restart` (7.97 s), `spray_pressure_outlet_restart` (11.33 s).

The pressure case enables molecular transport, WALE with coefficient 0.17 and
CI 0.1, drag and parcel outflow, with a 100000 Pa z-high outlet. Chemistry,
heating and evaporation are disabled in this fixture. It checks exact gas/parcel
CSV continuation and final history after stop/restart, existing mass/energy/
momentum budgets, parcel exit, changed/removed pressure checkpoint rejection
without history publication, and negative pressure rejection. Existing CI/model
guards also execute. Tolerances were unchanged. The zero-pressure existing case
also passes. No physical spray or reacting-outlet validation is claimed.

The first fixed build failed on a newly extended line exceeding 132 columns;
the argument list was wrapped and the build rerun successfully. Initial log:
`build-line-length.log`. Provenance, Python syntax, project contract and whitespace
checks passed. No full matrix, MPI application run or external solver parity ran.

## Commands and artifacts

Logs: `build/spray-pressure-outlet-checks/`
(`provenance.log`, `build.log`, `build-line-length.log`, `build-selected.log`, `ctest.log`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_(pressure_outlet_restart|isotropic_restart)$' --output-on-failure
python3 -m py_compile tools/check_spray_boundaries.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
