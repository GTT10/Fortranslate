# Serial EB KH boundary flight — 2026-09-26

Implementation `00169ffc`, parent `3415d0ed`, local main, version 0.246.0.
KH eligibility now previews the full private parent drift when explicit outer
faces are supplied. Periodic wrapping alone permits release; EB or reflecting
domain contact suppresses release for that source substep. Actual parent/child
flight and inventories publish through the existing transactional path.

## Focused checks

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug, existing fixed tests-enabled
`build/spray-amr-reactive-debug`. The fixed spray application rebuilt successfully.
Two direct checker configurations passed: periodic KH and isothermal-domain-wall
KH, both with restart. Each performed 3 successful executions (full, checkpoint
stop, continuation) and 10 expected rejections. No unexpected failures or skips.

Final populations were 34 parcels after periodic wrapping (y approximately
1.9999e-5 m) and 26 after reflection (y approximately 0.003980001 m).
Checks cover child IDs, component mass, velocity direction, positions and
byte-identical gas/particle/film CSVs after continuation. Restart rejection
checks include changed boundary/context, truncated payload and corrupted child
watermark. These fixtures do not independently exercise EB contact after wrapping;
they use the shared drift contact detection already checked at kernel level.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_eb_spray.py --exe build/spray-amr-reactive-debug/pelef_eb_spray_3d --source . --work build/serial-eb-boundary-khrt-checks/periodic --periodic --khrt --restart
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_eb_spray.py --exe build/spray-amr-reactive-debug/pelef_eb_spray_3d --source . --work build/serial-eb-boundary-khrt-checks/wall --wall --khrt --restart
python3 -m py_compile tests/application/check_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
ctest --test-dir build/spray-amr-reactive-debug -N -R '^spray_eb_(periodic|wall)_khrt_application$'
```

Provenance, syntax, project contract and whitespace passed. CTest lists both new
entries; application runs above are direct executions, not additional CTest runs.
No tolerances or golden data changed. Build warnings remain.
Artifacts: `build/serial-eb-boundary-khrt-checks/{build.log,configure.log,periodic.log,wall.log}`,
`periodic/eb_spray_f8xqbzwz`, and `wall/eb_spray_k9h7448k` beneath that directory.
Mechanism SHA: `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
PeleC: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
PelePhysics: `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.

No selected-mechanism run, full matrix, installed audit, external-solver parity
or physical breakup validation ran. Complete PeleC parity remains unproven.
