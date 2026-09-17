# MPI AMR spray checkpoint and continuation — 2026-09-18

Implementation: `ce0d92738c4aae1d8d51ae1fe8482a4fdea747fe`, based on
`8e3758a`. The checkpoint library combines existing distributed indexed gas
I/O with a rank-zero, globally ID-sorted particle archive. It preserves parcel
composition, lineage/deformation/KH-RT fields, time/step count, next injection
event, regrid status, ID watermark and injected inventory. Model tables, gas
configuration, spray options/schedule and prescribed geometries bind the context.
Continuation limits/output and checkpoint controls may change.

An exclusive `.pending` file reserves each immutable prefix. Gas and particle
payloads carry a common generation; only closed successful payloads permit a
POSIX hard-link publication of the final manifest. Failed saves may leave
reserved prefixes/partial payloads but no published completion marker. Existing
prefixes are not replaced. This provides atomic visibility on the shared Linux
filesystem, not power-loss durability. Reads keep candidates private until gas,
particle and migration checks succeed. Gas is repartitioned round-robin and
particles migrate to their new finest owners. Global particle archive capacity
is one million; gas payloads are not gathered.

The application lands on one configured checkpoint time, processes coincident
regrid/injection events, saves and optionally stops. Restart checks the saved
event/topology position and run limits before continuing. Absolute lexical path
checks protect input/output/checkpoint paths from ordinary naming collisions.

## Completed focused verification

Reused GNU Fortran 13.3.0 Debug/Ninja fixed tests-enabled and selected H2/O2
tests-disabled builds, `/usr/bin/gfortran`, Open MPI `/usr/bin/mpifort` and
`/usr/bin/mpiexec`. Compiler/MPI/backend configuration did not change. No install
prefix was involved.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_restart$' --output-on-failure --output-log build/spray-amr-checkpoint-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_amr_spray_restart.py --exe build/spray-amr-reactive-selected-debug/pelef_mpi_amr_spray_3d_selected --mpiexec /usr/bin/mpiexec --source /home/gtt13/projects/pelec --work /home/gtt13/projects/pelec/build/spray-amr-checkpoint-checks/selected --selected
python3 -m py_compile tests/application/check_mpi_amr_spray_restart.py tests/application/check_mpi_amr_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final fixed restart CTest **1/1 passed**, no skips, **4.48 s** test / **4.52 s**
total. Selected restart checker **passed**. Each checker completes a two-rank
continuous run, two-rank save/stop and four-rank continuation, then verifies three
expected rejections: changed cone speed/context, truncated particle footer, and
overwrite of an existing prefix. The overwrite case checks every existing file
byte-for-byte. Rejections produce no gas CSV.

Three injections occur at 0, 10 and 15 ns. Saving at 10 ns occurs after the
second injection and coarsening; continuation applies only the remaining event.
At 20 ns, all gas CSV values and all exported particle values agree with the
continuous run at relative 2e-12 / absolute 1e-25; 12 particles and three applied
events confirm schedule continuation. Gas chemistry/transport are disabled in
this CLI fixture. Breakup fields are serialized, but active-breakup continuation
is not independently exercised here.

An earlier two-event restart fixture and the existing application regression
passed **2/2** in **9.86 s** (restart 3.98 s, application 5.84 s). Adding the
future third event passed restart 1/1 in 4.40 s; the final run above followed
path-collision hardening. No build/test attempt failed in this increment;
expected rejection subprocesses are part of passing tests. Provenance, Python
syntax, project contract and whitespace checks passed. Compiler warnings remain.
No crash/fault injection, direct in-memory read-rollback test, install audit,
external solver parity, convergence, scaling or physical validation was run.

Logs: `build/spray-amr-checkpoint-checks/` includes build attempts,
`ctest-initial.log`, final `ctest.log` and `test-selected.log`.
Artifact directories are printed by each checker: fixed runs live below
`build/spray-amr-reactive-debug/tests/spray_mpi_amr_restart/restart_*`; selected
runs below `build/spray-amr-checkpoint-checks/selected/restart_*`.
