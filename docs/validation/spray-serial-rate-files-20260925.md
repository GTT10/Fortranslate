# Regular serial spray rate files — 2026-09-25

Implementation `1c50ed72bc0a37e30673477369ab5fee4515600d`, based on `1430a9e`.

The fixed/selected regular serial driver now loads external shared or independent
port profiles before existing input validation. The optional `spray_rate_files`
group is not serialized; loaded tables feed existing event integration and
checkpoint identity. Inline/file conflicts and inactive-port paths reject.
Existing shared/independent scheduling restrictions remain. No physical model,
event-generation rule or tolerance changed.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug; existing fixed tests-enabled and selected H2O2
tests-disabled configurations. Both application targets built; selected execution
was not run. CTest `spray_rate_files_application` **1/1 passed**, 3.28 seconds
(total 3.32), running two Python test methods. Zero failures/skips.

The new method checks shared and independent-port file input. Gas and parcel CSVs
match the inline case and resumed case byte-for-byte. Renamed files/comments
preserve continuation; changed numerical rates reject restart. Negative rates
and inline/file conflicts reject startup. Existing conservation checks remain
active. The second method reruns shared inline profile integration, analytic
injected mass, zero-dose events, restart and invalid-input cases.

Provenance, Python syntax, project contract and whitespace checks passed.
No failed builds/tests occurred. No full matrix, external solver parity,
physical validation, scaling or installation audit was run.

## Commands and artifacts

Logs: `build/spray-serial-rate-files-checks/` (`provenance.log`, `build.log`,
`build-selected.log`, `ctest.log`, `contract.log`). The existing unittest harness
uses temporary per-method directories and removes their inputs/CSVs/checkpoints
at teardown; CTest retains its test log.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_rate_files_application$' --output-on-failure
python3 -m py_compile tests/spray/test_spray_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
