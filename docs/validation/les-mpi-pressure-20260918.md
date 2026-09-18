# MPI LES pressure outlets — 2026-09-18

Implementation `90077a2b431985b493a93bdca4ce293f587c7988`, based on `6fa738d`.

`fill_mpi_spray_ghosts` accepts optional species metadata and six outlet pressures.
Positive pressure requires species and an outflow face. Pressure range, finiteness
and exact rank agreement are checked collectively. Pressure corrections act on
private ghost candidates after ordinary boundary filling and before exchanging
x halos, including their transverse layers. Interior x faces receive no pressure
correction; their halos come from the owning neighbor. Results publish only after
successful filling and exchange. `advance_mpi_spray_les` forwards pressure to both
stages and checks agreement even at zero dt/eddy coefficient. Supplying nonzero
pressure without physical faces is rejected. Coupled/application controls remain
pending. This reuses the existing characteristic boundary formula.

## Completed verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; configurations unchanged.
Fixed tests-enabled `test_les_mpi_models` and selected H2O2 tests-disabled
`pelef_selected_spray_runtime` built. Selected runtime was not run.
Final CTest **2/2 passed**, 6.79 seconds, zero final failures/skips:
`spray_les_mpi_models_2` (3.57 s), `spray_les_mpi_models_4` (3.18 s).

Each run covers three closures, CI zero/positive/CI-only variants and four
boundary configurations: periodic, outflow plus moving isothermal walls,
pressure x-outlet plus those walls, and distinct pressures on all six outflow
faces. Existing serial parity bounds remain unchanged: state relative `2e-12`
with unit floor, temperature absolute `1e-8 K`, boundary ledger relative `1e-10`
with `1e-12` scale floor. Rejections cover rank-dependent pressure at a disabled
step, one-rank NaN pressure and missing physical faces, preserving state and
temperature exactly. Existing CI/model disagreement and dt rollback checks pass.

The initial 2/2 run failed in the serial reference because the new all-outflow
fixture retained nonzero wall velocity. The fixture now sets wall velocity only
for cases containing walls; no numerical guard or tolerance was relaxed.
Initial failures are retained in `ctest-invalid-fixture.log`.
Provenance, Python syntax, project contract and whitespace checks passed.
No full matrix, external solver parity or physical validation was performed.

## Commands and artifacts

Logs: `build/les-mpi-pressure-checks/`
(`provenance.log`, `build.log`, `build-selected.log`, `ctest.log`,
`ctest-invalid-fixture.log`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_models --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_mpi_models_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
