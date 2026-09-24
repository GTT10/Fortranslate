# Regular MPI spray rate files — 2026-09-25

Implementation `24a45a139eb5f07dc1a2430d00b0d11f9e6b4bba`, based on `7c4555f`.

The fixed/selected regular MPI application accepts optional shared and per-port
rate files. Rank 0 loads and validates them before serializing the existing input
for broadcast. Other ranks receive the numerical tables, not file paths.
Existing event integration and content-bound checkpoint identity apply.
File/inline conflicts and inactive-port paths reject. No new physical closure,
live file reloading or automatic schedule generation is introduced.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing fixed tests-enabled
and selected H2O2 tests-disabled builds. Both application targets built; selected
execution was not run. CTest `spray_mpi_rate_files_application` **1/1 passed**
(6.86 seconds, total 6.90), zero failures/skips.

The checker first runs the existing default case with gas/species, mass/energy
and particle checks. It then compares inline and external shared/per-port rate
tables and a 2-to-4-rank restart. Gas and parcel fields match at existing relative
`2e-11`, absolute `1e-25` tolerance. Integrated injection is `3e-10 kg` at relative
`1e-14`. Renaming/comment changes preserve restart; changed numerical rates reject
restart. Negative rate input rejects before history output. No tolerances changed.

Provenance, Python syntax, project contract and whitespace passed. No failed
builds/tests occurred. No full matrix, physical validation, external solver
parity, installation audit or scaling check was run.

## Commands and artifacts

Logs: `build/spray-regular-rate-files-checks/` (`provenance.log`, `build.log`,
`build-selected.log`, `ctest.log`, `contract.log`). Inputs/checkpoints/CSVs:
`build/spray-amr-reactive-debug/tests/spray_mpi_rate_files_application/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_rate_files_application$' --output-on-failure
python3 -m py_compile tests/spray/test_mpi_spray_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
