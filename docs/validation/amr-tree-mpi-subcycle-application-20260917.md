# MPI tree application subcycling — 2026-09-17

## Source and implementation

Base `56f293aab0e10442a08d1dc0a8a57b4a586a9748`, local `main`.
The shared fixed/selected MPI application accepts `amr_tree_timestepping` and
selects recursive composition plus scaled hydro/transport/wall timestep bounds.
Reservoir bounds use root spacing in recursive mode. Checkpoint context appends
the same mode token as serial; default false retains existing common-mode identity.
Unknown tree groups and malformed mode controls reject explicitly. Mechanism data,
serial algorithms, tolerances and golden hashes are unchanged.

SHA-256:

- `src/driver/mpi_amr_patch_tree_application_3d_mod.F90`:
  `763aac5c5f56d6793bce870b561c0dc46c75b0a0f5c68d43a9b817880590f83c`
- `tools/check_mpi_amr_tree_restart_3d.py`:
  `bb307468b48547376e9e7cfcf87e9cd878dfa01b52259fff74b88e312f878d58`

## Commands and scope

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, existing
`build/selected-tree-mpi-debug`; no compiler/backend changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target pelef_mpi_amr_tree_reactive_3d pelef_amr_tree_reactive_3d test_selected_elementary_mpi_amr_tree_3d test_selected_elementary_amr_tree_3d -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_mpi_amr_tree_(subcycle_application|selected_application)_3d$' --output-on-failure
# After correcting only the new fixture, rerun only its failed regression:
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_mpi_amr_tree_subcycle_application_3d$' --output-on-failure
python3 -m py_compile tools/check_mpi_amr_tree_restart_3d.py tools/check_mpi_amr_tree_application_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final fixed/selected builds passed. The existing selected common-mode regression
passed in 6.87 s. The corrected recursive application regression passed in
33.79 s (CTest total 33.82 s). Both targeted regressions passed, with no final
failures or skips. Provenance, Python syntax, project contract and whitespace
checks passed. Existing compiler warnings remain.

The new `--subcycling` restart fixture uses a 4-cubed root and one 2-cubed
ratio-two child, with the existing reacting/transport hotspot and moving no-slip
walls. CFL 0.002 yields two coarse steps to 1e-9 seconds. Both fixed and selected
elementary H2/O2 frontends compare serial and MPI continuous CSV/budget bytes
and step counts, stop at step one on two ranks and resume on three ranks and
serial. Serial and MPI checkpoint bytes must match, and a one-rank MPI restart
from the serial checkpoint must reproduce outputs. The three-rank run includes
an empty owner. Malformed subcycling and a changed restart mode reject without
altering outputs, along with existing physics/ledger/truncation/alias guards.

The initial new fixture used a 2-cubed root, rejected by the application's
minimum grid-size contract. A subsequent CFL 0.004 trial completed in one step
and failed the checker's continuation requirement. The final input preserves the
4-cubed grid and CFL 0.002; no application validation was weakened.

Logs: `build/amr-tree-mpi-subcycle-app-checks/{build,initial-ctest,single-step-ctest,ctest}.log`.
Case artifacts: `build/selected-tree-mpi-debug/tests/mpi_amr_tree_subcycle_application_3d/`.
The existing selected common-mode application regression also checks branched
and root-only serial/2/3-rank output equality and input rejection.

## Limits

The new recursive application regression is static and uses PCM hydro. Dynamic
scheduled/automatic recursive combinations, PLM at application level and inflow
bound dispatch need separate focused qualification. Lower-level recursive tests
cover PCM/PLM, inflow and moving/isothermal walls. No full matrix, installed audit,
external solver parity, physical validation or production scaling was run.
Root-gather output and whole-patch round-robin ownership remain unchanged.
