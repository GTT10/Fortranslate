# MPI AMR spray physical-boundary application/restart — 2026-09-18

Implementation: `c68aa5aacbba43aab3ffcfa3c893a97f558476e9`, based on `a636921`.

Fixed and selected frontends share the updated driver. Resolved gas faces now
reach the reacting spray wrapper; derived periodic axes reach injection,
regrid, adaptation tagging, balancing and checkpoint particle migration.
Reservoir composition and state reuse existing gas routines, with current-tree
finest-spacing hydro/transport bounds. Gas wall temperature and tangential
velocity reach transport; particle walls remain specular.

Accepted steps accumulate a separate boundary inventory, printed in conserved
mass/momentum/energy/species order. The checkpoint stores it in the existing gas
budget field, validates finite rank agreement on write, and publishes it with
all restored state on successful read. The optional boundary output cannot be
omitted when reading a nonzero budget. Nonperiodic axis flags additionally bind
the generation context; the gas configuration already binds face/reservoir/wall
settings. Particle payload V1 is unchanged. No blanket cross-build checkpoint
compatibility is asserted.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja, `/usr/bin/gfortran`,
`/usr/bin/mpifort`, `/usr/bin/mpiexec`. Reused unchanged fixed tests-enabled and
selected H2/O2 tests-disabled configurations.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_(physical_application|inflow_application|restart)$' --output-on-failure --output-log build/spray-amr-open-checks/ctest-final.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_amr_spray_restart.py --exe build/spray-amr-reactive-selected-debug/pelef_mpi_amr_spray_3d_selected --mpiexec /usr/bin/mpiexec --source . --work build/spray-amr-open-checks/selected --selected --multiport --physical
python3 -m py_compile tests/application/check_mpi_amr_spray_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **3/3 passed**, zero failures/skips: reservoir inflow **4.66 s**,
wall/outflow **4.62 s**, periodic restart regression **4.61 s**; total **13.93 s**.
The selected wall/outflow restart checker also passed, recorded separately from
CTest. Provenance, syntax, contract and whitespace checks passed. Builds passed
with existing warnings.

Physical fixtures use two ports, four events (including tied events), a
prescribed coarsening and a 10 ns checkpoint followed by continuation to 20 ns.
Wall reflection occurs before saving; subsequent high-face outflow removes eight
parcels. The reservoir variant keeps the first port inside with incoming gas.
Two-rank uninterrupted results match two-rank save/four-rank continuation in
all gas/particle CSV values and cumulative boundary inventory at relative
2e-12/absolute 1e-25. Eight surviving particles, outlet identities/compositions,
positive surviving x velocities and signed boundary mass are asserted. Changed
physical faces or spray properties, truncated particle payloads and attempts to
overwrite immutable prefixes are rejected; original files remain intact.

Initial attempts: the first fixed physical fixture rejected a zero radius with
a finite solid/hollow nozzle model. The fixture now retains its positive outlet
radii. Initial manual selected and inflow runs could not find topology because
`--source .` remained relative after subprocess directory changes; the checker
now resolves that root. These failures were setup issues; no tolerances were
changed. Subsequent final checks above passed.

Artifacts: `build/spray-amr-open-checks/{build,build-selected,ctest-final,selected-final}.log`;
initial `ctest.log`, `selected.log`, `inflow.log` retain failures. Additional
intermediate successful checks are `ctest-physical.log` and `inflow-final.log`.
Final selected cases: `build/spray-amr-open-checks/selected/restart_3xme0fu0/`.
Fixed generated cases are below `build/spray-amr-reactive-debug/tests/` in the
three named CTest work directories, each with inputs, logs, CSVs and checkpoints.

These are application/restart checks with synthetic liquids and chemistry/gas
transport disabled. Library coupled transport evidence remains separately
recorded. Thermal/moving-wall operation, physical-boundary automatic adaptation
and load balancing were connected but not newly exercised here. Wall films,
EB impacts, convergence, external solver parity and physical validation remain
open; no full matrix or production scaling claim is made.
