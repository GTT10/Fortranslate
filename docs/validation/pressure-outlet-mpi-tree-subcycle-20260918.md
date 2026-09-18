# MPI patch-tree pressure outlet recursive hydro — 2026-09-18

Implementation `813d2fe`, base `a57d04c`.

The MPI level hydro API accepts optional `outlet_pressure(6)`, validates its
finite/range/face constraints, and includes normalized values in collective
control agreement before updating private candidates. Each owned patch at the
active depth applies the shared pressure ghost helper after time-interpolated
coarse/peer sampling. Both RK2 stages use the pressure ghosts. The recursive
hydro wrapper forwards targets to every level call; state, temperature, boundary
budget, maximum correction and optional level counters publish only on success.

MPI transport, reactive wrappers and tree application inputs remain pending.
This is implementation parity, not physical validation or scaling evidence.

## Focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6, `/usr/bin/mpiexec`; reused Debug/Ninja fixed
tests-enabled and selected H2O2 tests-disabled configurations.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_subcycle_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_amr_patch_tree_subcycle_3d_[0-9]+$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_tree_reactive_3d_selected --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Two/three-rank recursive tests add PCM/characteristic-PLM pressure cases to the
existing periodic, outflow, slip-wall and inflow matrix. Three nodes use ratios
2 then 3 with 1/2/6 level steps. X-low/z-high targets are 102325 Pa, y faces are
slip walls, and a refined child touches x-low. Owned state/temperature entries,
boundary budgets, reflux maxima and counters equal the serial reference exactly.
Ranks without patches still participate. A rank-local pressure mismatch or
pressure on a wall rejects collectively, retaining state, temperature, budget,
maximum correction and counters. Existing solver/dt failures remain covered.
No tolerances changed. Selected MPI executable compilation is separate from
these fixed-mechanism runtime comparisons.

Artifacts: `build/pressure-outlet-mpi-tree-subcycle-checks/{provenance,build,build-selected,ctest}.log`.

Final CTest **2/2 passed**, zero skips, 31.03 s total (2 ranks 15.55 s;
3 ranks 15.43 s). Fixed/selected builds, provenance, Python syntax, project
contract and diff checks passed. No failed or interrupted test runs.
