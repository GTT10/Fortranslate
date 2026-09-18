# MPI KH and splash child ID ordering — 2026-09-18

Implementation `026dac631e79fe1231adafec398d77e4e040c93a`, based on `86685ee`.

The shared MPI KH-release and supplied-wall-event splash routines now use the
same parent/sibling ordering as computational splitting. The private metadata
helper is renamed `assign_parent_child_ids`. Initial reservation still checks
watermark agreement and overflow. Children construct in private candidates,
then receive owner-independent IDs before publication. Construction was inspected:
child velocities, mass fractions and energy/work are independent of assigned
IDs, so reassignment does not alter the physical quadrature. Allocation or
assignment rejection preserves caller state. No physical-model formula changed.

This guarantees ID ordering for the same set of parent release events. It does
not make different event/substep groupings identical, change injection ordering,
or retroactively renumber saved IDs. Future child IDs can differ from older
builds. Other independent child constructors are outside this increment.

## Completed verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing configurations unchanged.
Fixed `test_spray_mpi` and selected H2O2 `pelef_selected_spray_tree_runtime`
built. Selected runtime was not run.
CTest **2/2 passed**, 17.20 seconds, zero failures/skips:
`spray_mpi_migration_2` (8.28 s), `spray_mpi_migration_4` (8.88 s).

KH and supplied-wall-impact fixtures place parent IDs in reverse rank order
with an empty final owner. They check exact child/parent IDs, four-child release,
global watermark, component inheritance and existing mass/momentum/energy
budgets. Existing capacity, late-event and whole-interval rollback assertions
also pass. At four ranks the new formula differs from the old rank-prefix
formula. Existing computational parent-order checks still run in the same
executable. No tolerance changed. No separate AMR/EB application restart,
full matrix, external solver parity, scaling or physical validation was run.

Provenance, Python syntax, project contract and whitespace checks passed.
No failed builds/tests occurred. An initial KH-only build preceded the final
KH-plus-splash build; reported tests used the final executable.

## Commands and artifacts

Logs: `build/spray-kh-id-checks/{provenance,build,build-final,build-selected,ctest}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_migration_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
