# MPI patch-tree pressure outlet common hydro — 2026-09-18

Implementation `7429dce`, base `6a9e1bf`.

The common-step MPI tree hydro API accepts optional `outlet_pressure(6)`.
Absent targets normalize to zero. Finite/range/face validation and collective
value agreement precede state updates. After distributed ordinary ghost sampling,
each owner applies the shared serial tree pressure helper to its private patch
halos before both RK2 flux stages. Physical-face selection uses tree geometry;
interior peer/coarse-fine boundaries retain their existing sampling.
Failures are agreed collectively before reflux or state/diagnostic publication.

Level/subcycle hydro, transport, reactive wrappers and MPI tree application
input still require forwarding. This is implementation parity, not physical
validation, nonreflection or large-tree scaling evidence.

## Focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6; reused Debug/Ninja fixed tests-enabled and
selected H2O2 tests-disabled builds. MPI launcher `/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_amr_patch_tree_hydro_3d_[0-9]+$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_tree_reactive_3d_selected --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Two- and three-rank tests add PCM/characteristic-PLM cases with x-low/z-high
102325 Pa targets and slip y walls on a three-node, three-level tree. The x-low
boundary touches a refined child. Every owned conserved and temperature entry,
boundary budget and maximum reflux correction equals the serial reference
exactly. Ranks with no owned patches participate in collective checks.
Rank-local pressure mismatch and pressure assigned to a wall reject without
changing state, temperature or diagnostics. Existing periodic/inflow/outflow,
invalid-solver and timestep-mismatch cases remain; no tolerances changed.

Final CTest **2/2 passed**, zero skips, 9.99 s total (2 ranks 5.04 s;
3 ranks 4.90 s). Fixed/selected builds, provenance, Python syntax, project
contract and diff checks passed. No failed or interrupted test runs.
The selected executable was built but not runtime-qualified here.

Artifacts: `build/pressure-outlet-mpi-tree-hydro-checks/{provenance,build,build-selected,ctest}.log`.
