# Two-level EB application step retry — 2026-09-17

Implementation: `7debd28ffca8c9a483448d6ffc4eb78e3991e7ed`, based on
`94982e7` on local `main`.

The serial two-level application now halves a rejected trial interval and
retries the complete reaction/transport/hydro composition. The existing full
update publishes state and per-step budgets only on success. Time, accepted
step count, cumulative inventories and checkpoint cadence advance only after
acceptance. `maximum_step_retries` defaults to 12 and accepts 0–60; zero keeps
immediate rejection. Each new accepted step starts with the configured maximum
interval clipped by remaining time. A nonadvancing time or exhausted retry
limit stops before output publication. The printed rejected-trial count is
local to the invocation, including after restart.

The retry limit controls termination, not the accepted arithmetic, so it is
not a checkpoint compatibility field. The saved maximum interval and minimum
accepted interval retain their existing meanings. This adds no automatic CFL
estimate, temporal-error estimator, small-cell stabilization, chemistry solver
change or physical qualification.

## Independently run checks

GNU Fortran 13.3.0, Debug/Ninja, existing
`build/selected-tree-mpi-debug`; configured `/usr/bin/mpifort` Open MPI 4.1.6.
The exercised application is serial. Commands:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target pelef_reactive_eb_two_level_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^regression_reactive_eb_two_level_3d$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/check_reactive_eb_two_level_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 regression passed**, zero failures/skips, 2.98 s (3.01 s total).
Provenance, syntax, project contract and whitespace checks passed.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

The added hydro-only case requests a 1e-3 s maximum interval and final time.
It completes 9 accepted steps after 18 rejected trials, with normalized
inventory residual `3.0085911406930120e-16`. The first accepted interval is
1.25e-4 s. Its state and budget CSVs are byte-identical to direct execution
of that interval with retries disabled, testing that rejected trials did not
leak state or inventory. Restart after the first accepted step also reproduces
both complete final CSVs byte-for-byte. Retry exhaustion and invalid negative
limits preserve existing outputs. The original three-step reacting/transport
case and checkpoint corruption tests remain in this same checker.

Logs: `build/eb-two-level-retry-checks/{build,ctest,retry}.log`; the exploratory
`retry.log` is not counted as an additional regression. Public regression
artifacts: `build/selected-tree-mpi-debug/tests/reactive_eb_two_level_3d/`,
including `retry.chk`, `retry_direct.nml.log`, `retry_resume.nml.log` and
`exhausted.nml.log`.
No broader configuration matrix, MPI execution, installed audit, external
solver parity or physical validation was rerun.
