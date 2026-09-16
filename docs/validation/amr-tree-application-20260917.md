# Serial regular 3D patch-tree application — 2026-09-17

Source parent: `2f7384a64b1ad13e4e943db5e8a3d52da1596a36` on `main`, plus this change.
PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Frontend SHA256: `017925d3b179333f88850803de35a8cef881a3b1ab7aaf19cabb5363d343d46a`.
Driver SHA256: `5e97bd95574b7b6a5f0ea6361dc3102f7b7f61aad42c026b4609c46263192a2a`.
Regression script SHA256: `d07b5308ad8a5ef8b53efb696265500eec4cdebac6c8d9036c4f24d3551b58a2`.

## Scope

New `pelef_amr_tree_reactive_3d` frontend and generic driver connect explicit
static topology input, node-local analytic initial conditions, synchronized
parent temperatures, common-timestep R/2-T/2-H-T/2-R/2, physical wall/inflow
options, rejection retry and exposed-cell CSV plus composite budgets.
The frontend loads fixed elementary/full-H2/O2 data; the regression exercises
elementary H2/O2. The CMake target is registered for installation, but no installed
binary audit was performed. Dynamic regridding, MPI, checkpoint/restart,
selected-mechanism dispatch, subcycling and PLM coarse/fine donors remain open.

## Completed checks

Reused `build/spray-mpi-migration-gnu-debug`: GNU Fortran 13.3 Debug, bounds/FPE
checks, OpenMPI 4.1.6, MPI/spray enabled. No toolchain/backend change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_amr_tree_reactive_3d pelef_reactive_3d -j 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^regression_amr_patch_tree_application_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Final build and **1/1 focused regression passed, 0 failures/skips, 2.28 s**.
The regression comprises:

- A three-level/four-node reacting, diffusive, moving-wall case. Its 127 exposed
  CSV cells cover exactly the domain volume; integrated output agrees with
  driver budgets, Euler boundary balances hold, species close, and chemistry
  changes species beyond the boundary contribution.
- Root-only hydro evolved for two steps to `1e-7 s`, compared by every regular
  CSV column against the existing regular-grid executable.
- An invalid parent reference rejected without replacing existing output.

The first fixture used two root cells per direction, rejected by the existing
configuration's minimum-four requirement; the fixture was corrected. A budget
check initially normalized almost-cancelled y momentum by its near-zero global
sum. The observed difference was about `3.36e-27`; the check now uses the
absolute cell inventory as the cancellation-aware scale with the unchanged
`2e-10` relative threshold. No solver tolerance or golden hash was changed.

Provenance/project contract and whitespace checks passed. Python AST checks of
provenance, contract and the new regression script passed (3/3).
Build/CTest logs: `build/amr-tree-app-checks/`.
Executable logs and CSVs:
`build/spray-mpi-migration-gnu-debug/tests/amr_tree_application_3d/`.
This is focused implementation evidence, not physical or scaling validation.
