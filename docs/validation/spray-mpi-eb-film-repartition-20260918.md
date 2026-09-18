# MPI EB film remap and rebalancing — 2026-09-18

Implementation `490c4a7`, base `993a373`.

`rebalance_mpi_amr_spray` now allows film state with `nodal_phi`, forwarding the
original nodal field to film migration. `remap_mpi_eb_spray` accepts optional
`wall_films`, `film_mass`, `maximum_films` together, agrees their presence across
ranks, and migrates private films after gas and particles. Neither operation
publishes any proposed state until all migrations succeed.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, `/usr/bin/mpifort` and `/usr/bin/mpiexec`,
Debug/Ninja. Existing fixed tests-enabled and selected H2O2 tests-disabled
configurations were reused.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_remap test_spray_mpi_eb_balance test_spray_mpi_amr_balance --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(eb_remap|eb_balance|amr_balance)_[24]$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The first fixed build failed because a test call exceeded Fortran's line limit.
A CTest invocation launched after that failure used stale binaries and is **not
counted**. After splitting the call, the fixed build succeeded and the final
CTest passed **6/6**, no skips, total **3.13 s**. EB balance: 0.46/0.47 s at
2/4 ranks; EB remap: 0.51/0.53 s; regular AMR balance: 0.54/0.57 s.
Selected library compilation, provenance, Python syntax, project contract and
whitespace passed. Selected execution was not run.

Both EB fixtures add a mixed water/peroxide stationary film on phi=y-0.00125,
initially on the final rank. Rebalancing moves it to the fine owner. Remapping
moves the refined patch to the opposite half, coarsens to the root, then refines
the original half while retaining the original finer nodal field. Film position,
normal, mass, temperature, diameter, height, material ID and component masses
are checked exactly after migration. Existing gas-remap conservation remains
at relative 3e-13, scale floor 1e-30; particle identity/composition are preserved.

Film capacity zero rejects after private gas/particle movement. Original owners,
particle/component populations, films, gas and temperature remain unchanged;
remap also retains the original patch placement. Existing particle-capacity
rejection remains in both tests. The regular AMR balance regression also passes.

Artifacts: `build/spray-mpi-eb-film-repartition-checks/` contains `build.log`
(the line-limit failure), `build-final.log`, `ctest.log` (excluded stale run),
`ctest-final.log`, `build-selected.log`, `provenance.log`, and `contract.log`.

These are library connections. MPI EB application film output/checkpoint and
routing through these APIs remain open. No broad matrix, installed audit,
external parity, scaling or physical validation was run.
