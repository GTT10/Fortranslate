# Dynamic MPI AMR spray coupling and input — 2026-09-25

Implementation `d41ff58261cd2e3d22c0694c5e6b0e019e4ac78f`, based on `50b81e0`.

The coupled MPI AMR spray driver now enables dynamic LES by model selection,
including when static coefficients are zero. Filter settings and periodic-only
boundaries validate before particle work; effective model/filter options agree
collectively. Existing private gas/particle/film state and LES half-step ledger
retain the full coupled transaction. The application now forwards dynamic options
to timestep selection and coupled stepping. Separate `amr_spray_dynamic` input
and context serialization preserve existing static checkpoint formatting while
binding dynamic defaults and nondefaults on restart.

## Completed checks and pending verification

GNU Fortran 13.3.0 / Open MPI 4.1.6, Ninja Debug; existing fixed tests-enabled and
selected H2O2 tests-disabled builds. Fixed `pelef_mpi_amr_spray_3d` and
`test_spray_mpi_amr_reactive`, and selected `pelef_mpi_amr_spray_3d_selected` built.
The initial fixed build failed on an overlong argument-list line; wrapping it
fixed the build. No numerical tolerance was changed.

Provenance preflight, Python syntax, project contract and whitespace passed.
**No application or MPI numerical test was executed.** The earlier CTest and
launcher PIDs 2018710/2018713 remained present on initial observation; no duplicate
MPI launch was started. The registered `spray_mpi_amr_dynamic_application` case
is prepared, not passed. It selects dynamic LES with zero static coefficients,
uses a Gaussian filter, compares continuous and changed-rank restart, and changes
the grid filter ratio to exercise context rejection. Existing malformed/context
checks remain in that checker.

This change supplies implementation/build evidence only. Do not claim validated
dynamic MPI AMR spray behavior from it. Nonperiodic dynamic boundaries, LES
subcycling, external parity and physical validation remain outside this result.

## Commands and artifacts

Logs: `build/les-dynamic-amr-application-checks/` (`provenance.log`, `build.log`
initial failure, `build-final.log`, `build-selected.log`, `contract.log`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d test_spray_mpi_amr_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
python3 -m py_compile tools/check_project_contract.py tests/application/check_mpi_amr_spray_restart.py
python3 tools/check_project_contract.py --source .
git diff --check
```
