# Serial EB spray checkpoint/restart — 2026-09-18

Implementation: `bec2b46166909ffd2f10699e102d0abc57966862`, based on `484d5c2`.

The EB application reuses the transactional spray archive with complete parcel
fields, component masses, last ID and event counter. An optional chemistry
ledger adds `_CHEM_V1` to the schema and a separate record; existing callers
without that argument retain their format. The source ledger stores wall
exchange, alongside exterior, injected and initial-gas inventories. Context
contains canonical gas settings, spray namelist, species/transport/reaction
tables and the original nodal level set. Output/continuation limits and save
controls are excluded. Context capacity remains 131072 characters.

Snapshots are published after accepted steps and all tied injections, as
`<checkpoint_file>.step<step>`. Existing destinations are rejected. Restart
validates archive/EOS, schedule position, retained count/IDs, capacity and fluid
locations before stepping. Final-time and maximum-step limits must admit the
saved state. Existing checkpoint mechanics stage reads and publish closed files
atomically on Linux; no power-loss durability claim is made.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; reused unchanged fixed
tests-enabled and selected H2/O2 tests-disabled configurations.
Compiler `/usr/bin/gfortran`, MPI wrapper `/usr/bin/mpifort`.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_(restart|schedule)$' --output-on-failure --output-log build/spray-eb-restart-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_integration --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_integration_checkpoint$' --output-on-failure --output-log build/spray-eb-restart-checks/ctest-legacy.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_restart$' --output-on-failure --output-log build/spray-eb-restart-checks/ctest-final.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_eb_spray.py --exe build/spray-amr-reactive-selected-debug/pelef_reactive_eb_3d_selected_spray --source . --work build/spray-eb-restart-checks/selected --selected --scheduled --restart
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/application/check_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Initial focused CTest **2/2 passed**, total **1.79 s**: restart **1.21 s**,
schedule **0.54 s**. Legacy checkpoint regression **1/1 passed**, **0.03 s**
(total **0.07 s**). After adding changed-geometry rejection, the final restart
test passed **1/1**, **1.36 s** (total **1.40 s**). These reruns are not additional
distinct cases. Selected scheduled restart passed independently. Provenance,
syntax, contract and whitespace checks passed. An initial build invocation
requested nonexistent target `test_spray_checkpoint`; the actual
`test_spray_integration` target was then built and checked successfully.
No runtime test failures or skips; no tolerance changes.

The schedule has tied events at 5 ns (including zero dose), followed by 10 ns
and 20 ns events. Save occurs after the 5 ns injections. Restart reaches 20 ns
with gas and particle CSVs byte-identical to uninterrupted output and all four
printed cumulative ledgers identical. Tests reject immutable overwrite,
changed injection temperature, modified nodal geometry and truncated payload.
Fixed/selected fixtures use their existing distinct gas compositions and do not
claim inter-mechanism parity. Primitive/breakup archive fields are saved, but
this application restart fixture does not activate breakup.

Artifacts: `build/spray-eb-restart-checks/` contains build/CTest/checker logs;
the final selected case is `selected/eb_spray_gzx4vj_n/`. Fixed cases remain under
`build/spray-amr-reactive-debug/tests/spray_eb_restart/`. Rejected/corrupt test
files are retained beside valid snapshots.

This is bounded serial restart, not distributed EB particle migration, a
cross-compiler guarantee, install audit or physical validation. Large contexts
are explicitly rejected; scalable geometry identity/I/O remains open.
