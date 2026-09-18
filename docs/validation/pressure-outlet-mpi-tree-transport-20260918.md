# MPI patch-tree pressure outlet transport — 2026-09-18

Implementation `606c685`, base `b7e8bff`.

Common-step and level-specific MPI tree molecular transport accept optional
`outlet_pressure(6)`. Normalized targets participate in collective control
agreement after finite/range/face validation. Owners apply the shared tree
pressure helper after ordinary viscous ghost sampling, before both RK2 transport
flux stages. The recursive wrapper forwards targets to every level. Existing
private candidate publication, flux limiting and reflux remain unchanged.
MPI reactive wrappers and application input remain pending.

## Focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6, `/usr/bin/mpiexec`; reused Debug/Ninja fixed
tests-enabled and selected H2O2 tests-disabled configurations.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_transport_3d test_mpi_amr_patch_tree_transport_subcycle_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_amr_patch_tree_transport(_subcycle)?_3d_[0-9]+$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_tree_reactive_3d_selected --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Both common and recursive tests add x-low/z-high 102325 Pa targets with moving
isothermal y walls and PCM/linear coarse donors. Three-node trees have ratios
2 then 3, with recursive counts 1/4/36. Two/three-rank checks compare every owned
state and temperature entry, boundary budget, reflux maximum and limiter minimum
to serial exactly. Empty ranks participate. A rank-local pressure mismatch rejects
without changing state, temperature, diagnostics or recursive counters; existing
owner-local invalid transport and timestep mismatch checks remain.

Initial fixed build failed because a test-edit replacement incorrectly added
pressure arguments to ALLOCATED/DEALLOCATE. Corrected before running tests;
`build-final.log` records the successful rebuild. No tolerances changed. This is
implementation parity, not physical validation or scaling evidence. The selected
MPI executable build is separate from fixed-mechanism runtime qualification.

Artifacts: `build/pressure-outlet-mpi-tree-transport-checks/` contains
`provenance.log`, initial `build.log`, `build-final.log`, `build-selected.log`
and `ctest.log`.

Final CTest **4/4 passed**, zero skips, 175.54 s total. Common transport
took 7.34/6.94 s at 2/3 ranks; recursive transport took 81.23/79.99 s.
Fixed/selected builds, provenance, Python syntax, project contract and diff
checks passed. No failed or interrupted runtime tests.
