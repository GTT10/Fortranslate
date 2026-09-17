# MPI AMR film and reacting-wrapper composition — 2026-09-18

Implementation: `f06eaa876b556b452bc16f404c75531a275dcf12`, based on `48917d4`.

The optional film arrays/settings in `advance_mpi_amr_spray_reactive` wrap the
existing particle/gas split in F/2–S/2–gas–S/2–F/2. Arguments must occur together;
film presence agrees collectively before any optional film collective is entered.
Both film halves use the existing AMR interval, adding their wall-energy ledgers
to gas/particle boundary exchange. All state, IDs, diagnostics and inventories
publish only on complete success. Separate film/particle counts return zero on
failure. Film-free calls retain the original operator sequence.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja, `/usr/bin/gfortran`,
`/usr/bin/mpifort`, `/usr/bin/mpiexec`. Reused unchanged fixed tests-enabled and
selected H2/O2 tests-disabled build configurations.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_film_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_film_reactive_(2|4)$' --output-on-failure --output-log build/spray-amr-film-reactive-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_restart$' --output-on-failure --output-log build/spray-amr-film-reactive-checks/ctest-app.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

New composition tests **2/2 passed**: two ranks **20.82 s**, four ranks **19.89 s**,
total **40.76 s**. Film-free application/restart regression **1/1 passed**,
**4.14 s** (CTest total **4.18 s**). No failures or skips. Fixed/application and
selected-library builds, provenance, Python syntax, contract and whitespace
checks passed. Existing warnings remain. No build/test attempt failed.

The three-level test retains two mixed films at opposite physical walls with
active heating and evaporation. Common and recursive gas updates enable
viscosity, conduction and species diffusion. Independent composite gas-plus-film
inventories, including wall heat and gas boundary flux, close at 2e-12 of
`max(1e-12,abs(initial))` per component. Films heat, lose mass and retain exact
wall positions. Exhausting the first film interval and rejecting the gas solver
after prior film exchange both preserve original gas/temperature, film
state/geometry/composition and sentinel ledger; diagnostics remain unchanged.

The new fixture has no airborne parcels and disables chemistry to isolate film
composition with hydro/transport. Existing chemistry/particle checks are separate
historical evidence, not newly rerun here. Selected evidence is compilation only.
The application regression uses airborne parcels without films, including changed
rank restart and checkpoint rejection checks.

Artifacts: `build/spray-amr-film-reactive-checks/{build,build-app,build-selected,ctest,ctest-app}.log`.
New executable: `build/spray-amr-reactive-debug/tests/test_spray_mpi_amr_film_reactive`.

Impact-created AMR films, whole gas/particle/film remap, application controls and
film checkpoint state remain open. This is not EB impact, lateral film transport,
physical validation, convergence, external parity or production scaling evidence.
