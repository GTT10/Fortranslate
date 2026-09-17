# Reactive EB gas / spray composition — 2026-09-18

Implementation: `55fe6493e195e4ae2d5373acd37cab1f8d4a7270`, based on `e26dcee`.

The serial wrapper applies a particle half-interval, the existing EB
reaction/transport/hydro operator, and a second particle half-interval to private
state. Publication includes gas, temperature, parcels, component masses,
particle-step count, transport limiter and three ledgers. Exterior gas exchange,
combined gas/particle EB wall exchange and gas chemistry remain separate.
All use the inventory-change sign. Failure preserves fields and ledgers and
returns zero particle steps, including when gas evolution rejects after the
first successful particle half. The selected spray tree library now explicitly
links the selected EB runtime that supplies the gas operator.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; existing fixed tests-enabled and
selected H2/O2 tests-disabled configurations reused. Compiler
`/usr/bin/gfortran`, MPI wrapper `/usr/bin/mpifort`; no backend or compiler change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_(reactive|step)$' --output-on-failure --output-log build/spray-eb-reactive-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_mechanism_provenance.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **2/2 passed**, total **0.20 s**: coupled reactive test **0.14 s**,
existing particle-interval test **0.02 s**. Selected library compilation,
provenance, syntax, contract and whitespace checks passed. Initial test compilation
failed on two overlong source lines; continuation lines fixed them. No runtime
failures or skips and no tolerance changes. Existing compiler warnings remain.

The coupled fixture uses a 4 mm cube of 1 mm cells, x=1.25 mm planar EB, 1000 K
N2 gas with nonzero velocity, and a 350 K water/peroxide parcel 1 nm from the
wall. It advances 10 ns with Rusanov, StateRedist, implicit chemistry, viscosity,
conduction, species diffusion and barodiffusion enabled. Assertions require
heating, evaporation, drag, reflection, positive wall impulse, a fluid endpoint,
nonzero chemical species change, valid transport limiter and covered-cell
preservation. Gas-plus-liquid inventories close against boundary+wall+chemistry
at relative 3e-12 with a 1e-30 scale floor. An invalid gas solver exercises
rollback after the first particle half. The liquid property tables are synthetic.

Artifacts: `build/spray-eb-reactive-checks/build.log` (initial compiler failure),
`build-final.log`, `build-selected.log`, `ctest.log`.

This is a serial single-grid composition. It inherits the first-order particle
source/drift split and root-contained trajectory restriction. The caller chooses
the gas timestep; this wrapper does not adapt or retry it. A symmetric outer
split does not prove second-order convergence. EB injection/restart application
integration, MPI/AMR particle evolution and liquid wall films remain open.
Selected execution, convergence, external solver parity and physical validation
were not run.
