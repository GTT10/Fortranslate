# Rank-neutral EB particle payload — 2026-09-18

Implementation `0b960097cac5f3f9002bd41c2447325fd3e84642`, base `6901977`.

`mpi_eb_particle_checkpoint` is a particle-only persistence API. It preserves all
parcel fields (including parent/material IDs, deformation, RT clock and KHRT
state), component masses, time, accepted step count, next event, ID watermark and
injected inventory. Context/path, operation and dimensions agree collectively;
write metadata must agree exactly across ranks. Private particles are validated
against the current nodal geometry and gathered in ID order on rank zero.
Formatted records use 19 significant decimal digits. A closed `.pending` file
is linked to an immutable final path. Failed publication leaves its staging
file and never replaces an existing final archive.

Reads check magic/context, bounds, finite metadata, particle/component validity,
strictly increasing IDs, watermark, footer and EOF. Private restored records
migrate to the supplied current EB tree owners; state and metadata publish only
after successful collective migration. The caller owns context binding to the
gas generation, geometry, liquid tables, schedule and model. No gas state is
saved by this API, and application restart remains disabled pending integration.

## Completed focused checks

GNU Fortran 13.3.0 `/usr/bin/gfortran`, Open MPI 4.1.6 `/usr/bin/mpifort` and
`/usr/bin/mpiexec`; Debug/Ninja. Existing unchanged fixed tests-enabled and
selected H2/O2 tests-disabled configurations were reused.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_checkpoint --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_particle_checkpoint$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tests/application/check_mpi_eb_particle_checkpoint.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **1/1 passed**, test **1.50 s**, total **1.54 s**. The script saves on two
ranks, reads on four ranks with a changed fine owner, and rejects a truncated
copy on four ranks. It verifies every parcel field exactly, both component
masses, time/step/event, watermark and injected ledger. Nonzero parent/material,
deformation/RT/KHRT state is included. It also rejects overwrite and wrong
context; failed reads preserve empty incoming particles and sentinel metadata.
Selected compilation, provenance, Python syntax, contract and whitespace passed.
No failures or skips. The initial passing fixture used default breakup state;
the final fixture strengthened this coverage and removed a large automatic
character buffer warning by allocating the buffer. These are revisions of one
test, not additional distinct evidence.

Artifacts: `build/spray-mpi-eb-checkpoint-checks/{build-final.log,build-selected-final.log,ctest-final.log}`;
the CTest script reports its retained archive/log directory.

This is payload restart evidence, not combined gas/particle continuation or an
application restart claim. Atomic namespace publication is not an fsync/power
loss guarantee. Whole-population root gathering is capped at one million and
has not been scaling-qualified. No selected runtime test, full matrix, external
parity or physical validation was run.
