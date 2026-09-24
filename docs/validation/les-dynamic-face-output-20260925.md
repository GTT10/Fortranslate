# Dynamic LES periodic face outputs — 2026-09-25

Implementation `bd1a7ec5bddbbb037147c3708b177089ad658d21`, based on `52d87f1`.

`compute_dynamic_les_rhs_3d` can return all three full periodic face-flux arrays.
It retains the existing numerical fluxes and RHS evaluation. Array shapes are
checked before work; partial requests and cropped `active_x` workspaces reject.
Private face storage publishes with RHS only after successful evaluation.
Mass/species fluxes remain zero and periodic low/high end faces agree.
This is a prerequisite interface for conservative coupling, not an AMR Dynamic
LES implementation: coarse/fine halos, filtering and reflux are still absent.

## Completed checks and pending MPI attempt

GNU Fortran 13.3.0, Ninja Debug, existing configurations. Fixed
`test_spray_core`/`test_les_mpi_dynamic` and selected H2O2
`pelef_selected_spray_runtime` built successfully. Selected execution was not run.

Two direct serial invocations completed with exit 0: `les_dynamic_flux` and
`les_dynamic_step`. They verify unchanged RHS, exact face-divergence reconstruction,
periodic endpoints, zero species fluxes, rejection of partial requests and
preservation of RHS/faces on invalid gas state, together with existing conservation,
translation, dissipation and SSPRK2 checks. No tolerances changed.

The combined CTest attempt had not completed at this record: it was waiting in
`/usr/bin/mpiexec -n 2 .../test_les_mpi_dynamic`, with no child ranks observed.
CTest PID was 2018710 and launcher PID 2018713. TERM and then KILL were sent only
to that launcher; its process still appeared running with a pending KILL signal
at the last observation. Do not count this as MPI success, completed termination,
or a completed CTest run. Recheck these handles before any repeat launch.

Provenance, Python syntax, project contract and whitespace checks passed.
No full matrix, external solver parity, AMR execution or physical validation was
performed. The MPI runtime issue is not established as a numerical-code failure.

## Commands and artifacts

Logs: `build/les-dynamic-faces-checks/` (`provenance.log`, `build.log`,
`build-selected.log`, `ctest.log`, `serial-flux.log`, `serial-step.log`, `contract.log`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core test_les_mpi_dynamic --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
# Attempted; pending at this record:
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_(dynamic_(flux|step)|mpi_dynamic_2)$' --output-on-failure
# Completed independently:
env -u LD_LIBRARY_PATH -u LIBRARY_PATH build/spray-amr-reactive-debug/tests/test_spray_core les_dynamic_flux
env -u LD_LIBRARY_PATH -u LIBRARY_PATH build/spray-amr-reactive-debug/tests/test_spray_core les_dynamic_step
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
