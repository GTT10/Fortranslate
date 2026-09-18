# Serial EB computational splitting — 2026-09-18

Implementation `602cde8`, base `9191dd7`.

The common splitter now accepts optional shared IDs and component masses. New
IDs exceed the supplied watermark even when higher-numbered parcels have retired.
It copies per-drop component masses and partitions multiplicity and KH inventories.
Capacity, invalid component inventory or ID overflow failures publish nothing.
The allocated serial EB reactive wrapper splits after both particle halves and
gas evolution, before committing any state. Fixed-population callers reject an
enabled cap. The application uses the allocated interface for both wall paths.

Checkpoint APIs optionally retain cumulative split IDs using a `_SCOUNT_V1`
format suffix and an integer record after the watermark. Read validates the
counter privately and publishes it with the other state. The application checks
injection IDs plus split IDs plus four-child splash groups. Cap-disabled callers
retain their existing checkpoint layout. The optional counter is forwarded via
an allocated scalar only when splitting is enabled.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja. Reused fixed tests-enabled and
selected H2O2 tests-disabled build directories; all executed tests were serial.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_impact_reactive test_spray_core --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_split|spray_eb_impact_reactive)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_eb_split_application|spray_eb_restart|spray_eb_app_splash)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
python3 -m py_compile tests/application/check_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Core/reactive tests **2/2 passed**, total **1.78 s**. The reactive case verifies
late capacity failure preserves gas, temperature, particles, component masses,
IDs and ledgers. Success with and without film arguments turns one parcel into
five, retaining ID 1 and allocating IDs 101–104 above retired watermark 100.
Composite inventories close within relative 3e-12, scale floor 1e-30. The core
test also rejects invalid components and integer overflow and checks descendant
component ancestry. Existing reactive deposition/splash checks remain active.

Application tests **3/3 passed**, total **4.73 s**: scheduled splitting/restart,
legacy splash/restart and legacy specular/restart. Their drivers performed 29
launches including expected negative cases. The split fixture has tied, zero-dose
and final-time injection events. It ends with 12 records: evolved parcels obey
the cap of 60; the final two freshly injected parcels remain unsplit. The step-two
checkpoint records four computational child IDs. Continuous and resumed gas,
particle, film CSVs and printed cumulative budgets are byte-identical. Altering
the saved split count by one rejects restart. Existing immutable-publication,
changed-context, truncated-payload, alias and schedule checks remain enabled.

Fixed and selected builds, provenance, syntax, project contract and whitespace
checks passed. No unexpected test failures or skips. Selected execution, broad
configuration matrix, installed audit, external solver parity and physical
validation were not run. Serial KH shed-mass child creation remains open.

Build/check logs: `build/spray-eb-split-checks/` (`provenance.log`,
`build-final.log`, `ctest.log`, `build-app.log`, `ctest-app.log`,
`build-selected-final.log`, `contract.log`). Application artifacts under
`build/spray-amr-reactive-debug/tests/`:

- `spray_eb_split_application/eb_spray__cxsm7b7`
- `spray_eb_app_splash/eb_spray_2dbf3yzm`
- `spray_eb_restart/eb_spray_q_3bxpea`
