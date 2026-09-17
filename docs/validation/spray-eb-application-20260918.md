# Fixed/selected serial EB spray applications — 2026-09-18

Implementation: `1bd6c18059e807af6c5e69599dee468821484d10`, based on `1d85862`.

The shared driver reads the existing EB gas configuration and an `eb_spray`
namelist, initializes one multicomponent cone event, advances the coupled solver
under hydro/transport/user timestep bounds, retries failed steps by halving, and
writes gas/particle CSVs plus cumulative exterior/wall/chemistry inventories.
The configured geometry builder optionally returns its original nodal field
transactionally for level-set geometries, so gas and particles share the same
input. Fixed and selected frontends use their respective mechanism loaders.
Output alias guards cover both CSV paths against all three input files after
lexical absolute normalization. Checkpoint/restart requests are rejected.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; existing fixed tests-enabled and
selected H2/O2 tests-disabled configurations reused. `/usr/bin/gfortran` and
`/usr/bin/mpifort`; unchanged compiler and optional backends.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d test_spray_eb_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_(application|reactive)$' --output-on-failure --output-log build/spray-eb-app-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_eb_spray.py --exe build/spray-amr-reactive-selected-debug/pelef_reactive_eb_3d_selected_spray --source . --work build/spray-eb-app-checks/selected --selected
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/application/check_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **2/2 passed**, total **0.53 s**: application **0.34 s**, coupled library
regression **0.15 s**. The selected checker passed independently. Provenance,
syntax, contract and whitespace checks passed. No build/test failures or skips;
existing compiler warnings remain. A preliminary manual fixed run passed; its
particle CSV format was subsequently given explicit component column names and
the redundant trailing comma removed before automated checks.

The case has a planar nodal EB at x=1.25 mm in a 4 mm cube, two water parcels
starting 1 nm from the wall at 20 m/s toward it, total represented liquid mass
1e-10 kg, 10 micrometre drops and 350 K liquid. Fixed gas uses the existing
elementary default composition; selected gas uses named pure N2 with the H2/O2
bundle. Both enable chemistry and molecular transport and reach 20 ns. Checks
require two finite particle records, reflection into fluid, heating, evaporation,
component mass consistency, 64 gas records at final time and finite ledgers.
An absolute-output/relative-input alias is rejected without changing the spray
input bytes. The checker does not claim fixed/selected numerical parity.

Artifacts: `build/spray-eb-app-checks/` has `build-final.log`,
`build-selected-final.log`, `ctest.log`, `selected.log`, preliminary build/run
files and `selected/eb_spray_7mznu16s/`. Fixed CTest cases are retained under
`build/spray-amr-reactive-debug/tests/spray_eb_application/`.

This is an initial-injection, serial, single-grid application with a first-order
particle source/drift split. Source liquids are synthetic. CSVs are not restart
files. Timed injection, coupled checkpoint/restart, outer-domain particle events,
wall films and MPI/AMR ownership remain open. No install audit, full matrix,
convergence, external parity or physical validation was performed.
