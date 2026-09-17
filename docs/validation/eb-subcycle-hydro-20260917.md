# Two-level EB subcycled Euler hydro — 2026-09-17

Implementation: `03b608e815fc84250be3f4ea92b33026adad72ee`, based on
`4459d33` on local `main`.

`advance_reactive_eb_two_level_subcycled_euler_3d` advances the coarse grid once
and the fine patch `ratio` times over the same interval. Each fine update uses
PCM face donors from a linear interpolation of conserved coarse state between
the initial and unrefluxed coarse endpoint. Halo construction recovers donor
EOS temperature; interpolated temperature is only its solve guess. The final
fine interval consumes the remaining coarse time to avoid interval drift.

The routine accumulates accepted fine flux integrals and embedded-wall changes
at every fine step. It then restricts fluxes, refluxes unrefined coarse cells,
averages the final child state into covered parents and computes composite
boundary/wall ledgers. All four state arrays and both output ledgers publish
only after success. Existing common-step hydro and transport paths retain their
arithmetic. The new path shares their geometry consistency and state checks.

Scope: one static fine patch, serial, physical outflow, Euler hydro. This does
not increase temporal order or connect reaction/transport subcycling,
application timestep selection, restart mode, multiple patches or MPI. Callers
must bound the coarse interval and every fine interval appropriately; the
application still uses common-step coupling.

## Independently executed checks

GNU Fortran 13.3.0, Debug/Ninja, existing
`build/selected-tree-mpi-debug`. The exercised tests are serial.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target test_eb_two_level_3d test_eb_two_level_full_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^(unit_eb_two_level_3d|unit_eb_two_level_full_3d)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, zero failures/skips, 9.67 s total: two-level hydro 0.98 s,
existing full reaction/transport/hydro composition 8.65 s.
Provenance, project contract and whitespace checks passed; no Python changed.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

The expanded hydro test exercises ratios 2 and 3, an interior patch and one
meeting physical domain boundaries. Three subcycled coarse intervals of
1e-7 s each check all-component composite inventory balance against boundary
and wall ledgers, parent/child volume-weighted agreement, positive temperatures
and preserved covered storage. Fine evolution differs from a single fine
Euler step over the same interval. Large nonphysical intervals, inconsistent
geometry and invalid ratios reject without changing state or ledgers. The
original common-step cases remain in the same test. This does not specifically
inject a failure after a successful fine substep.

Logs: `build/eb-subcycle-hydro-checks/{build,ctest}.log`.
No full matrix, installed audit, application subcycling, external solver parity,
convergence study or physical validation was rerun.
