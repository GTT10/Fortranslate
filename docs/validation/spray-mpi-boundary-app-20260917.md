# MPI application physical boundaries and restart — 2026-09-17

Source: parent `9e8607a` plus this change. Source hashes and logs:
`build/spray-mpi-boundary-app-checks/`. Existing GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration. Fixed H2/O2 application exercised here; the selected
frontend shares the implementation but was not independently rerun in this change.

The `mpi_spray` namelist exposes six `boundary_faces`, `wall_temperature(6)` and
`wall_velocity(3,6)`. Inlets use the configured initial gas state. Physical labels
and wall data pass through the full coupled step. Injection wraps periodic axes
only, rejecting positions outside nonperiodic faces. Tangential gas wall velocity
requires viscosity; isothermal walls require conduction. Particle walls retain
stationary specular reflection, without mixed deposition/splash or film evolution.

Physical-boundary history adds cumulative `boundary_mass,boundary_energy`, positive
into gas plus particles. Error columns subtract injection and boundary increments.
The complete boundary inventory persists through the existing checkpoint
boundary/film schema with zero films. Face and wall settings extend checkpoint
context; changing them rejects restart. Default periodic history columns and
checkpoint context remain unchanged. A runnable example is
`cases/mpi_spray_3d/boundary_mixture.nml`.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_mpi_spray_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_application_[24]$' --output-on-failure
python3 -m py_compile tests/spray/test_mpi_spray_application.py
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Final application entries passed 2/2 in 85.82 s, with no failures or skips.
Python syntax, provenance, contract and whitespace checks passed. Detailed outputs:
`build/spray-mpi-migration-gnu-debug/tests/spray_mpi_application_2/run-ubgy49ie/`
and `build/spray-mpi-migration-gnu-debug/tests/spray_mpi_application_4/run-gczjfgx1/`.

New application cases use open x faces, no-slip/isothermal y faces, tangential
wall motion, SGS and four escaping parcels. The final population is zero and
boundary energy is nonzero. Same-rank restart reproduces every gas/particle
CSV value exactly; 2-to-4 and 4-to-2 restart agrees within relative 2e-11 and
absolute 1e-25. Cumulative boundary budgets use those same tolerances, and corrected
mass/energy errors remain below 2e-11. A changed wall temperature rejects
before history creation. Existing periodic, injection, breakup, SGS and restart
application checks remain in the same two entries.

The first compile caught an overlong history-header statement; splitting it fixed
compilation. No numerical tolerance or reference hash was relaxed. This is focused
input/output and restart evidence, not physical spray validation, production
scaling or a full configuration matrix.
