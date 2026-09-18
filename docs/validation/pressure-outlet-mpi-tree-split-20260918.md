# MPI patch-tree pressure outlet coupled splitting — 2026-09-18

Implementation `980bab5`, base `45ba7b2`.

The MPI full split API forwards optional `outlet_pressure(6)` to both transport
halves and hydro, including recursive branches. Chemistry acts on private
candidates, and publication still follows completion of the full
R/2–T/2–H–T/2–R/2 interval. Pressure validation/agreement is performed by the
flow/transport operators before their updates. MPI application input and the
standalone reaction/hydro Strang API remain separate pending connections.

## Focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6; reused Debug/Ninja fixed tests-enabled and
selected H2O2 tests-disabled configurations, `/usr/bin/mpiexec` launcher.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_split_3d test_mpi_amr_patch_tree_subcycle_split_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_amr_patch_tree_(subcycle_)?split_3d_2$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_tree_reactive_3d_selected --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Common and recursive tests add a distributed three-node case with x-low/z-high
203650 Pa targets, moving isothermal y walls, characteristic PLM, implicit
chemistry and all molecular transport operators. Conserved states, temperatures,
boundary budgets, reflux maxima and limiter minima are compared exactly with
serial. Existing post-source hydro failure, enabled-operator mismatch and
recursive-mode mismatch retain rollback checks. Original zero-target cases
remain. No tolerances changed; this is implementation parity, not physical
validation or scaling evidence.

The initial fixture accidentally used 102325 Pa against initial pressures above
202650 Pa. The extrapolated ghost pressure became nonphysical, and the serial
oracle rejected the case. Targets were corrected to 203650 Pa; rejection gates
were preserved. Initial 2/3-rank results are in `ctest.log`; final reruns are
limited to two ranks for each coupled path in `ctest-final.log`.

Artifacts: `build/pressure-outlet-mpi-tree-split-checks/` contains provenance,
fixed/selected build logs and CTest logs. Selected executable compilation does
not constitute selected-mechanism runtime qualification.

Initial CTest **0/4 passed** (153.08 s), all rejected at the serial oracle
with the incorrect pressure fixture. Final two-rank CTest **2/2 passed**, zero
skips, 115.83 s total (common 10.17 s, recursive 105.62 s). Three-rank final
fixtures were not rerun. Fixed/selected builds, provenance, Python syntax,
project contract and diff checks passed. No interrupted runs.
