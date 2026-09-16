# Combined regular 3D tree reaction/transport/hydro update — 2026-09-17

Source parent: `a37f7857f4d2c575c0b46cec6c33dfb31d7a7e55` on `main`, plus this change.
PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Split module SHA256:
`af0caeb7c19870582c2ff71caf3c8dc65e17f39a1211c06af3fcbe3f9d928d61`.
Test SHA256:
`c8b759a83a8ac35e72f4348488494e62635d615382e1e80b5f5eb69245cda58c`.

## Scope

The new library entry point composes R/2-T/2-H-T/2-R/2 with the existing tree
operators. It sums transport/hydro boundary ledgers, takes the maximum stage
reflux correction and minimum transport limiter, and publishes all results only
after the complete composition succeeds. Reaction and transport can be disabled;
unused reaction/transport data are then not evaluated. Barodiffusion requires
species diffusion. The caller supplies the common timestep and handles retry.

This connects the operators in `pelef_core`; it does not provide the tree
application, MPI ownership, remapping, checkpoint/restart, PLM coarse/fine
interpolation or subcycle scheduler. It is not physical validation.

## Focused checks

Reused GNU Fortran 13.3 Debug/OpenMPI 4.1.6 build with bounds and FPE checks;
MPI/spray enabled. No toolchain change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_amr_patch_tree_split_3d -j 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^unit_amr_patch_tree_split_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

The test uses three levels/four nodes with adjacent root children and a perturbed
deep leaf, elementary H2/O2, explicit/implicit reaction integration and periodic/
outflow boundaries. It checks total mass/momentum/energy and H/O/N balance with
the combined boundary ledger, species closure, final restriction consistency,
and exact rollback of state/temperature/diagnostics after hydro rejection
following reaction and transport. It also checks active transport changes the
coupled solution, and disabling transport or both sources gives exact parity
with existing reaction/hydro or hydro entry points.

Initial compilation caught a long test declaration; it was wrapped before
execution. The first focused test passed. A second final run adds the active
transport comparison; no tolerance or golden-hash change was made.
Provenance, project contract, whitespace and Python AST parsing of both checker
scripts passed. Logs: `build/amr-tree-split-checks/build.log` and `ctest.log`.

Final result: **1/1 focused test passed, 0 failures/skips, 45.29 s**.
The final target build also passed. No full configuration matrix was run.
