# Distributed EB spray application connection — 2026-09-18

Implementation `c1ac24a3e440ad5322784000d19f5e72458806e4`, base `6e08e81`.

The existing EB tree application's shared source now has an opt-in spray target,
`pelef_mpi_eb_spray_3d`. Selected MPI EB tree bundle targets also produce a
`_spray` executable. Both accept tree and spray input paths. Initial mixed-liquid
cone parcels are migrated to their finest fluid owners before the first step.
The existing retry loop invokes the transactional gas/spray wrapper. Composite
budgets include liquid mass, momentum, kinetic/internal energy and component
species; output includes owner-local particle CSVs with global IDs. Gas-only
builds retain their existing preprocessor path.

## Independently completed checks

GNU Fortran 13.3.0 `/usr/bin/gfortran`, Open MPI 4.1.6 `/usr/bin/mpifort` and
`/usr/bin/mpiexec`, Debug/Ninja. Unchanged existing fixed tests-enabled and
selected H2/O2 tests-disabled configurations were reused. The shell's unqualified
`mpiexec` is a different installation; all executions below explicitly use the
configured `/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_mpi_eb_application|regression_mpi_eb_tree_application_3d)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_eb_spray.py --exe build/spray-amr-reactive-selected-debug/pelef_mpi_reactive_eb_tree_3d_selected_spray --mpiexec /usr/bin/mpiexec --source . --work build/spray-mpi-eb-application-checks/selected --selected
python3 -m py_compile tests/application/check_mpi_eb_spray.py tools/check_project_contract.py tools/check_mechanism_provenance.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **2/2 passed**, **50.20 s** total: new fixed application **21.74 s**,
existing gas-only application regression **28.42 s**. The selected application
script also passed. Provenance, Python syntax, contract and whitespace passed;
no failures or skips. The initial direct fixed script run also passed and is
not counted as an additional distinct test. Selected linking passed the existing
MPI runtime dependency audit (ABI 40).

Each new application script runs two and four ranks and compares every gas and
particle CSV value at relative 3e-12, floor 1e-30. It verifies exactly two unique
particle IDs, per-drop component mass, finite fields, reflected wall-normal
velocity, fluid-side positions and changed particle temperatures. A requested
restart is rejected before output. Fixed gas uses the existing elementary
initialization; selected gas uses pure N2. The fixture has a 4 mm cube, ratio-two
left child, y=1.25 mm planar nodal geometry and two water parcels, with implicit
chemistry and all molecular transport enabled for 20 ns. The application retains
its existing combined-budget check; its scale floor is one, so this is not an
independent tight small-mass conservation qualification. Tight coupled conservation
is covered separately by the preceding library test.

Artifacts in `build/spray-mpi-eb-application-checks/`: `build-final.log`,
`build-selected.log`, `ctest.log`, `application.log`, `application-selected.log`.
Direct run outputs: `fixed/eb_spray_cjowc2u7/` and
`selected/eb_spray_mi4kssun/`; the CTest fixture directory is printed in its log.

This connection supports initial injection and a fixed hierarchy. It explicitly
rejects checkpoint/restart, regridding and rebalancing until particle lifecycle
state is integrated. Timed events, independent ports, film impacts and child
creation are not provided here. Gas subcycling/automatic timestep controls pass
through but were not independently exercised by this new fixture. Particle CSVs
are analysis outputs, not restart files. No full matrix, installed-binary audit,
scaling, convergence, external solver parity or physical validation was run.
