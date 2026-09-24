# External MPI EB spray rate files — 2026-09-25

Implementation `3f8df2a8c3613861e008d6762a6d1f008a01c388`, based on `e12d48a`.

The fixed/selected MPI EB spray application now loads optional per-port files
through the shared transactional rate reader. A separate `eb_spray_rate_files`
namelist leaves existing input/context serialization unchanged when omitted.
Loaded tables enter the existing `eb_spray` context; paths and comments do not.
File/inline conflicts and paths beyond the configured port count are rejected.
Files are read at startup and must be consistent across MPI ranks. This change
does not connect file loading to the separate serial single-level EB application.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug;
existing fixed tests-enabled and selected H2O2 tests-disabled configurations.
Both application targets below built successfully; selected execution was not run.
CTest **2/2 passed**, zero failures/skips, 77.65 seconds total:

- `spray_mpi_eb_rate_files_application`: 45.52 seconds.
- `spray_mpi_eb_rate_velocity_application`: 32.09 seconds.

The new checker verifies two distinct rate files, inline/file gas and particle
agreement, 2/4-rank agreement and 2-to-4-rank checkpoint continuation. Renaming a
file and changing its comment preserves continuation; changing a numerical rate
rejects the checkpoint. A negative rate rejects startup before CSV output.
Existing analytic final-event mass and discharge-speed checks remain active,
including the zero-scale event. CSV tolerance remains relative `3e-12` with
scale floor `1e-30`; no tolerances changed.

Provenance, Python syntax, project contract and whitespace checks passed.
No failed builds or tests occurred. No full matrix, external solver parity,
physical validation, installed-binary audit or scaling check was run.

## Commands and artifacts

Logs: `build/spray-eb-rate-files-checks/` (`provenance.log`, `build.log`,
`build-selected.log`, `ctest.log`, `contract.log`). New application artifacts:
`build/spray-amr-reactive-debug/tests/spray_mpi_eb_rate_files_application/nozzle_8oz_pom7/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_(rate_files|rate_velocity)_application$' --output-on-failure
python3 -m py_compile tests/application/check_mpi_eb_nozzle.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
