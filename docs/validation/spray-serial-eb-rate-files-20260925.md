# Serial EB spray rate files — 2026-09-25

Implementation `0dc222a9b32dbb8e0a3973c280b51ce6b1e83da8`, based on `84f0efe`.

The shared fixed/selected serial EB application driver loads a common profile
file and per-port overrides through the existing transactional reader. The
optional `eb_spray_rate_files` group is separate from serialized `eb_spray`:
loaded values enter checkpoint identity, while paths/comments do not. The common
profile retains existing inheritance semantics. Inline/file conflicts and paths
beyond the active port count reject. Files load at startup; events and integration
intervals remain explicit inputs. No physical model or tolerance changed.

## Completed verification

GNU Fortran 13.3.0, Ninja Debug; existing fixed tests-enabled and selected H2O2
tests-disabled configurations. Both executable targets below built successfully;
the selected executable was not run. CTest **2/2 passed**, zero failures/skips,
2.73 seconds total:

- `spray_eb_rate_files_application`: 1.64 seconds.
- `spray_eb_port_profiles_application`: 1.05 seconds.

For both solid and hollow outlets, the new checker exercises a shared constant
profile inherited by port 1 and an independent linear port 2 file. Gas, particle
and film CSVs match inline input and checkpoint continuation byte-for-byte.
Renaming the common file/changing comments preserves continuation. Changed rate
values reject restart; a negative rate rejects startup before CSV output.
Existing analytic dose, discharge speed, property-fit, size and inventory checks
remain active. The second test checks the unchanged inline-only path.

Provenance preflight, Python syntax, project contract and whitespace checks
passed. No failed builds/tests occurred. No full matrix, external solver parity,
physical validation, scaling or installation audit was run.

## Commands and artifacts

Logs: `build/spray-serial-eb-rate-files-checks/` (`provenance.log`, `build.log`,
`build-selected.log`, `ctest.log`, `contract.log`). Per-run inputs, CSVs and
checkpoints: `build/spray-amr-reactive-debug/tests/spray_eb_rate_files_application/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_(rate_files|port_profiles)_application$' --output-on-failure
python3 -m py_compile tests/application/check_eb_spray_nozzle.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
