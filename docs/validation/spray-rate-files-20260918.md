# External AMR spray rate files — 2026-09-18

Implementation `ad57ce97fa42efaf739b0aefb6d81606b81a9459`, based on `83e58d4`.

The shared profile module now reads bounded two-column tables transactionally.
It accepts spaces/tabs, blank lines, `#`/`!` comments and numeric exponent forms;
rejects extra fields, nonnumeric/nonfinite values, long records and excess rows;
then applies the existing monotone-time/rate validator. Arrays and count publish
only after complete successful validation, with unused entries zero-padded.

The fixed/selected MPI AMR application loads optional per-port files before
integrating events. File/inline conflicts and unused port paths reject. Loaded
tables serialize through the existing rate namelist into checkpoint identity;
paths and comments do not. Shared consistent input files remain the caller's
MPI contract. This does not introduce live file reloading or automatic event
subdivision.

## Completed verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing configurations unchanged.
Fixed `pelef_mpi_amr_spray_3d` and `test_spray_rate_file`, and selected H2O2
`pelef_mpi_amr_spray_3d_selected` built. Selected executable was not run.
CTest **3/3 passed**, 9.02 seconds, zero failures/skips:
`spray_mpi_amr_rate_files_application` (4.91 s), `spray_mpi_amr_restart` (4.06 s),
`spray_rate_file` (0.01 s).

The application reads two distinct profiles with comments/tabs and preserves
the analytic injected mass `8.125e-11 kg`. Continuous 2-rank gas/parcel results
match 2-to-4 rank restart from another run directory at existing `2e-12` relative
and `1e-25` absolute tolerances. Changed curve values reject restart; negative
rates reject input. Existing context/truncation/overwrite checks pass, as does
the default-input case. The reader test checks valid zero-padding and exact
output preservation for an extra field, duplicate time, NaN, overlong record,
65 rows and a missing file. No tolerances changed.

Provenance, Python syntax, project contract and whitespace checks passed.
No failed builds/tests, full matrix, external solver parity, scaling or physical
validation occurred. Initial/final builds preceded the final test-target build;
reported tests used the completed reader implementation.

## Commands and artifacts

Logs: `build/spray-rate-files-checks/` (`provenance.log`, `build.log`,
`build-final.log`, `build-tests.log`, `build-selected.log`, `ctest.log`).
Application artifacts:
`build/spray-amr-reactive-debug/tests/spray_mpi_amr_rate_files_application/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d test_spray_rate_file --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_rate_file|spray_mpi_amr_(rate_files_application|restart))$' --output-on-failure
python3 -m py_compile tests/application/check_mpi_amr_spray_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
