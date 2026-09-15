# Regular coupled-step wall impacts — 2026-09-16

Source: parent `73ccfef` plus this change; GNU Fortran 13.3 Debug,
`build/port-spray-debug`, fixed H2/O2 mechanism, MPI disabled.

## Implemented boundary

Optional wall settings and persistent static film arrays connect the impact
transaction to both spray coupling halves. Crossings use the first intersected
plane, sampled gas pressure, selected liquid properties and overlapping static
patch heights. Parents and new children advance for the remaining interval;
film inventories join system mass/species/energy and wall impulses/work enter
the boundary ledger. All candidate state commits only with the full coupled step.

This API requires stationary walls and has no film flow, heat transfer or
re-evaporation. Impact occurs at the end of the particle substep; it is not exact
hit-time integration. Corner handling, scalable film spatial indexing and physical
qualification remain open. The existing wet leftover-deposition closure remains
explicit. CLI switches, film output and checkpoint persistence are not connected.

## Focused evidence

Seven CTest entries passed, zero failures/skips, in 17.01 s: wall event/state/
correlations, integration rollback/materials/wall coupling, application CLI.
The CLI contains the preceding 19 tests (16.05 s CTest wall time).
The new integration case uses 10 and 300 m/s droplets crossing the low-z no-slip
wall. Low speed deposits completely; high speed creates four children and a
film patch. It verifies gas/drop/film mass conservation, energy and momentum
against the surface ledger, child lineage, and rollback of gas, parent, children,
films and ledger after a deliberately invalid later chemistry backend.

After adding finite surface-tension and film vapor-index guards, the wall-coupling
entry was rebuilt and rerun successfully in 0.70 s. The initial build rejected an
overlong Fortran declaration; wrapping it resolved the compile failure.
Provenance, project contract and whitespace passed. No Python source changed.
No broad matrix or physical wall-impact qualification was run.

Artifacts: `build/port-spray-debug/wall-coupling-build.log` and the latest focused
rerun in `build/port-spray-debug/Testing/Temporary/LastTest.log`. Earlier seven-entry
results are recorded in this report. Reference code identities and wet closure
are in `references/spray_wall_impact.json`.

```sh
cmake --build build/port-spray-debug --target test_spray_integration test_spray_core pelef_spray_3d --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(integration_(wall_coupling|materials|rollback)|wall_(event|state|impact)|application_cli)$' --output-on-failure
cmake --build build/port-spray-debug --target test_spray_integration --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_integration_wall_coupling$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```
