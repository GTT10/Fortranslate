# Dynamic LES harmonic application controls — 2026-09-18

Implementation `18386ec5d21ed5c23d3e15ffd7600dfcbe48be8a`, based on
`83099e032c428827420389df3aeb9aae96fc8b1c`.

Serial/MPI spray coupling forwards optional `les_harmonic_mean` to both dynamic
SGS half steps. Fixed/selected application namelists expose the same logical,
defaulting to false. True requires the dynamic model. MPI coupling checks
agreement collectively, including when SGS is disabled. True appends a versioned
checkpoint identity suffix; false preserves the previous context. Restart
rejects changing true to false or removing the setting.

## Completed verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing configurations unchanged.
Fixed tests-enabled serial/MPI applications and selected H2O2 tests-disabled
serial/MPI applications built successfully. Selected binaries were not run.

Focused CTest **4/4 passed**, 22.75 seconds, zero failures/skips:

- `spray_dynamic_boundary_restart` (1.98 s)
- `spray_mpi_dynamic_restart` (8.88 s)
- `spray_dynamic_harmonic_restart` (2.16 s)
- `spray_mpi_dynamic_harmonic_restart` (9.70 s)

These check arithmetic/harmonic continuation, same-rank exact restart,
MPI 2-to-4-rank continuation within the existing tolerance, conservation budgets,
model/filter guards, and harmonic setting change/removal rejection without
history publication. Serial fixture chemistry is disabled. Provenance, Python
syntax, project contract and whitespace checks passed. No failed build/test ran.
No full matrix, external solver parity or physical validation was performed.
Nonperiodic and AMR dynamic LES remain outside this increment.

## Commands and artifacts

Logs: `build/les-harmonic-application-checks/`
(`provenance.log`, `build.log`, `build-selected.log`, `ctest.log`, `contract.log`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_spray_3d pelef_mpi_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_3d_selected pelef_mpi_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_(dynamic_boundary_restart|mpi_dynamic_restart|dynamic_harmonic_restart|mpi_dynamic_harmonic_restart)$' --output-on-failure
python3 -m py_compile tools/check_spray_boundaries.py tests/spray/test_mpi_spray_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
