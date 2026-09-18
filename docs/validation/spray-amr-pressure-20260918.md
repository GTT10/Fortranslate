# MPI AMR spray pressure outlet — 2026-09-18

Implementation `79236161b7aff18f7b3fa84e2e45c6088efe58e9`, based on `b48b8a1`.

Fixed and selected MPI AMR spray applications now accept the existing gas
`outlet_pressure` controls. The coupled API validates finite pressures in
`[0,1e15]`, allows positive values only on outflow faces, and requires exact
rank agreement before particle or gas candidate work. It forwards pressures
to the common/recursive gas split and both static LES halves. The application
also passes them to the ghost-aware LES interval estimate. Existing gas
checkpoint context already serializes these controls. No checkpoint schema or
particle outflow rule changes.

## Completed verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing configurations unchanged.
Fixed `pelef_mpi_amr_spray_3d` and selected H2O2
`pelef_mpi_amr_spray_3d_selected` built. Selected executable was not run.
CTest **1/1 passed**, zero failures/skips: `spray_mpi_amr_pressure_application`
(7.06 s; 7.10 s total).

The fixture combines a low-x slip wall, high-x pressure outlet at 100000 Pa
(initial gas pressure 101325 Pa), two injection ports, scheduled coarsening,
Vreman 0.17 and CI 0.1. A continuous 2-rank run agrees with 2-to-4 rank
checkpoint continuation for gas, parcels and cumulative boundary inventory at
the existing `2e-12` relative / `1e-25` absolute tolerances. Changed outlet
pressure rejects restart before CSV publication. Existing LES/context change,
payload truncation and checkpoint overwrite rejection checks also pass.
Chemistry and molecular transport are disabled in this short fixture; their
pressure paths were not separately rerun. Rank-disagreement rejection was
reviewed in code, not independently exercised by this application test.

Mechanism provenance, Python syntax, project contract and whitespace checks
passed. No failed builds/tests or changed tolerances. No full matrix, external
solver parity, scaling or physical validation was performed.

## Commands and artifacts

Logs: `build/spray-amr-pressure-checks/{provenance,build,build-selected,ctest}.log`.
Run artifacts: `build/spray-amr-reactive-debug/tests/spray_mpi_amr_pressure_application/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_pressure_application$' --output-on-failure
python3 -m py_compile tests/application/check_mpi_amr_spray_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
