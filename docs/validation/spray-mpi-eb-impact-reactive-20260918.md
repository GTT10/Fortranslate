# MPI EB reactive impacts and evolving films — 2026-09-18

Implementation `f375eee`, base `22f15b6`.

`advance_mpi_eb_spray_reactive` accepts an optional shared `last_id` together
with films, film component masses and wall options. This enables impact handling
in both particle source halves of the particle/film/reactive-gas/film/particle
composition. The watermark and all arrays remain private until full success.
A source model rejection propagates through optional `model_rejected` without
publication. Calls without a watermark retain the existing specular behavior.
The optional-state presence flags are agreed collectively before any source work.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`),
Debug/Ninja. Existing fixed tests-enabled and selected H2O2 tests-disabled builds
were reused. The new test enables runtime checks and invalid/zero/overflow traps.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_impact_reactive test_spray_mpi_eb_film_reactive test_spray_mpi_eb_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_(impact_reactive|film_reactive|reactive)_[24]$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **6/6 passed**, total **6.42 s**, no failures or skips. New impact tests:
two ranks 1.07 s; four ranks 1.06 s. Existing supplied-film and specular reactive
tests passed at both rank counts. Selected library compilation, provenance,
Python syntax, project contract and whitespace passed; selected execution was
not run. A source comment was clarified after the fixed run, with no executable
code changes.

The two-level planar fixture uses a mixed water/peroxide parcel, normal impact
speed 200 m/s, an initially empty film population, wall temperature 400 K and
interval 1e-8 s. Four children with IDs 2–5 and parent ID 1 finish on the fine
owner. The new deposited film heats above 350 K in the film source halves.
Implicit reaction, HLLC hydro and viscosity, conduction, species diffusion and
barodiffusion are enabled. Chemical species exchange is nonzero. Composite gas,
particle and film inventories close with wall, outer-boundary and chemical
ledgers to relative 3e-12 with scale floor 1e-30.

A deliberately invalid gas solver fails after the first private particle and
film halves. Original particle fields, gas/temperature, empty films, component
matrices, shared watermark and all diagnostic sentinels remain unchanged.
Rank-disagreed optional film presence also rejects collectively.

Artifacts: `build/spray-mpi-eb-impact-reactive-checks/{build.log,ctest.log,build-selected.log,provenance.log,contract.log}`.

Films created in the final particle half begin thermal evolution on a later
step. Existing empirical-impact, uniform-wall-temperature, bounded contact,
stationary-film and caller-identical property contracts remain in force.
MPI application inputs, output, regrid/rebalance and checkpoint integration of
these films remain open. No broad matrix, installed audit, external parity,
scaling or physical validation was run.
