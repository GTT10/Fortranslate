# MPI AMR film checkpoint and changed-rank restore — 2026-09-18

Implementation: `fedbc786b1758f1f7bdfe582d0c82e038d5fe2ee`, based on `078cd90`.

Optional film arrays/component masses/receive capacity form a collective mode
in `mpi_amr_spray_checkpoint`. Film mode binds the gas generation context and
uses particle payload V2: parcel records, film count, explicit film records,
then the existing footer. The manifest/file set is unchanged; no-film calls
retain particle V1. Film material, coordinates, normal, thermal/geometry state
and components round-trip through es27.18 text. Private collection is bounded
to one million films globally. Restored films use finest-cell ownership on the
reading communicator, preserving exact physical wall coordinates.

All particle/film records, metadata, footer/EOF and migration capacities validate
before publication. A film mode mismatch rejects instead of discarding films.
Existing gas storage carries cumulative boundary inventory. Immutable prefix
publication remains the existing POSIX hard-link protocol, not an fsync claim.

## Completed focused verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja, `/usr/bin/gfortran`,
`/usr/bin/mpifort`, `/usr/bin/mpiexec`; unchanged fixed tests-enabled and selected
H2/O2 tests-disabled configurations reused.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_film_checkpoint pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_(film_checkpoint_(2|4)|restart)$' --output-on-failure --output-log build/spray-amr-film-checkpoint-checks/ctest.log
# After fixing the fixture's temperature preparation:
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_film_checkpoint --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_film_checkpoint_(2|4)$' --output-on-failure --output-log build/spray-amr-film-checkpoint-checks/ctest-final.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final film checkpoint tests **2/2 passed**: two ranks **0.50 s**, four ranks
**0.53 s**, total **1.06 s**. Film-free application/restart regression passed
**4.10 s** in the first combined invocation. No skips. Fixed/selected builds,
provenance, syntax, contract and whitespace checks passed. Selected evidence is
compilation only; existing compiler warnings remain.

The three-level fixture writes two films at opposite walls, then each rank
independently reads that checkpoint with MPI_COMM_SELF. Every gas/temperature
value and film field/component matches exactly; saved time/steps/event/ID and
nonzero boundary energy restore. A subsequent world-communicator read checks
redistribution again. Immutable overwrite, zero film receive capacity,
rank-disagreed/omitted film mode and trailing archive corruption reject without
publishing state or sentinel metadata. Test-generated archives are retained;
the final corruption check intentionally leaves trailing data in `.particles`.

The first film tests failed the gas writer's strict temperature consistency
check because they supplied a literal temperature guess. Adding the same EOS
recovery used by the application during initialization resolved both tests;
writer tolerances and production code were not changed for this failure.
Failed prefix reservations remain as evidence.

Artifacts: `build/spray-amr-film-checkpoint-checks/{build,build-test,build-final,build-selected,ctest,ctest-final}.log`.
Test archives: `build/spray-amr-reactive-debug/tests/film_checkpoint_{2,4}_<tick>*`.
Executable: `build/spray-amr-reactive-debug/tests/test_spray_mpi_amr_film_checkpoint`.

This is library restart evidence. The new fixture has no airborne parcels;
existing particle application restart is checked separately above. Film-aware
application controls/continuation, EB impact, physical validation, convergence
and production scaling remain open.
