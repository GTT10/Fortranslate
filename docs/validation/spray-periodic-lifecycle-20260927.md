# Periodic EB spray lifecycle — 2026-09-27

Implementation `0ad2f38e`, based on `08ed07b118418f9d3ff33a5e94a8ea14b4baa017`,
local main, version 0.246.0.

## Change

EB coupled remapping and particle/film checkpoint migration now accept optional
periodic axes (omitted means nonperiodic). Load balancing forwards its existing
axes to EB particle and film migration instead of rejecting periodic EB input.
The application passes its paired boundary labels into all three operations.
The existing private-state transactions, component payloads, SDC reaction-source
history and archive format are retained. Checkpoint context remains caller-bound
to the gas/model; application gas restart checks boundary labels.

## Independently completed checks

Artifacts: `build/spray-periodic-lifecycle-checks/` (`provenance.log`, `build.log`,
`ctest.log`, `LastTest.log`). Existing fixed-mechanism Ninja Debug configuration,
GNU Fortran 13.3.0 and OpenMPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d test_spray_mpi_eb_balance --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_(periodic_balance_(2|4)|sdc_recursive_periodic_film_restart_regrid)$' --output-on-failure
python3 -m py_compile tests/application/check_recursive_eb_spray_impacts.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final result: **3/3 passed**, no skips. Periodic EB balance took 0.49/0.51 s on
two/four ranks, and the application restart/regrid check took 70.32 s; total
71.37 s. Provenance, build, Python syntax, project contract and whitespace passed.
An initial build caught an accidentally added periodic keyword on the generic
bin-exchange call; it was removed before testing. Existing warnings remain.

The balance fixture checks actual owner reassignment with particles, mixed films
and reaction history, plus rollback for rank-disagreed history presence, particle
capacity and film capacity. Periodic x/z faces are enabled in the new modes.

The application starts a parcel near the periodic z seam, wraps before EB splash,
forms a film, advances recursive gas/particle/film stages, then moves the refined
patch at time 1e-8. A two-rank checkpoint after that regrid is continued on four
ranks to 2e-8. Gas, particle and film CSVs **exactly match** the continuous
four-rank run. Both end with four particles and one film, normalized combined
budget residual **1.2425714608450861e-15** against the unchanged `3e-12` criterion.
Budget columns also match within that criterion; reading leaves all archive
bytes unchanged. The checkpoint explicitly contains film state and splash IDs.

Cases:
`build/spray-amr-reactive-debug/spray-sdc-recursive-periodic-film-restart-regrid/impact_0da5x_il/`.

The application fixture uses a planar EB; oblique-seam thermal evolution,
selected mechanisms, a full configuration matrix, physical qualification and
external-solver parity were not run. Library balance and application regrid/restart
are separate evidence; this does not claim measured application load-balancing
performance.
