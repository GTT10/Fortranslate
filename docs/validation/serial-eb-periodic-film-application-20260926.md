# Single-level EB boundary films in applications — 2026-09-26

Implementation `f0e88a91`, parent `b8da5392`, local main, version 0.246.0.
The serial particle interval forwards outer faces to the impact kernel and
periodic axes to film/gas sampling and deposition. Coupled reaction and application
input now accept EB films with periodic or domain-wall gas boundaries. Isothermal
gas walls map to specular particle reflection; no direct particle-wall heat model
or outer-domain film creation is added. Restart film validation uses the same
wrapped stencil. The gas boundary and spray context retain their existing identity.

## Focused checks

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug. Existing fixed tests-enabled
`build/spray-amr-reactive-debug` and selected tests-disabled
`build/spray-amr-reactive-selected-debug`, selected full H2/O2.

Three direct checker configurations passed, each with **3 successful executions**
(full, checkpoint stop, continuation):

- Fixed periodic deposition: 9 expected rejections.
- Selected periodic splash: 10 expected rejections, including the existing wet
  model-limit rejection.
- Fixed isothermal-domain reflection followed by EB deposition: 9 expected rejections.

Gas, particle and film CSVs are byte-identical after restart. Film positions
confirm wrapping before periodic deposition and reflection before wall deposition.
Existing finite/component checks, child IDs, boundary/context changes, immutable
checkpoint publication and truncated payload rejection execute unchanged.
There were no unexpected failures or skips. A fixed periodic-film CTest entry
was registered; direct runs are not counted again as CTest runs.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_eb_spray.py --exe build/spray-amr-reactive-debug/pelef_eb_spray_3d --source . --work build/serial-eb-periodic-film-app-checks/fixed --periodic --impacts deposit --restart
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_eb_spray.py --exe build/spray-amr-reactive-selected-debug/pelef_reactive_eb_3d_selected_spray --source . --work build/serial-eb-periodic-film-app-checks/selected --selected --periodic --impacts splash --restart
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_eb_spray.py --exe build/spray-amr-reactive-debug/pelef_eb_spray_3d --source . --work build/serial-eb-periodic-film-app-checks/wall --wall --impacts deposit --restart
python3 -m py_compile tests/application/check_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, contract and whitespace passed. Build warnings remain.
No tolerances or golden data changed. Artifacts under
`build/serial-eb-periodic-film-app-checks`: `build-fixed.log`, `build-selected.log`,
`fixed.log`, `selected.log`, `wall.log`, `fixed/eb_spray_789hxmew`,
`selected/eb_spray_1eftcmwa`, `wall/eb_spray_9r3f8bq2`.
Mechanism SHA: `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
PeleC: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
PelePhysics: `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.

No full matrix, installed audit, external solver comparison or physical impact/
film validation ran. Complete PeleC parity remains unproven.
