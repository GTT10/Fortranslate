# Synchronous regular 3D tree hydro — 2026-09-17

## Source and scope

Parent: `919f618860f478257016883fc1eff960ac1023fc` on `main`, plus this change.
PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
New module SHA256:
`5b7f8b47f7c3d19b3431b48852c5131d00184f0246b116a4faa627dabc8104fe`.
New test SHA256:
`31e42d1c31e7a5c4bb034ce837e82a3c35928c4f50b334b93745a51ce5d88741`.

This adds transactional common-timestep SSPRK2 hydro, all-node CFL selection,
PCM/PLM face fluxes, reflux, restriction, temperature recovery and physical
boundary accounting. Coarse/fine spatial interpolation remains PCM. Subcycles,
chemistry/transport, tree MPI/application/restart and physical validation are
not established by these tests.

## Completed focused checks

Reused `build/spray-mpi-migration-gnu-debug`: GNU Fortran 13.3 Debug,
bounds checking and invalid/zero/overflow traps; OpenMPI 4.1.6,
MPI and experimental spray enabled. No compiler/backend change.

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_amr_patch_tree_hydro_3d test_amr_patch_tree_3d test_amr_patch_tree_ghost_3d test_amr_patch_tree_reflux_3d -j 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^unit_amr_patch_tree(_ghost|_reflux|_hydro)?_3d$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Final focused CTest: **4/4 passed, 0 failures, 0 skips, 8.02 s**.
Provenance, project contract, whitespace and Python AST parsing of both check
scripts passed. Logs are under `build/amr-tree-hydro-checks/`.

The new test checks root-only PCM/PLM against the existing serial SSPRK2 kernel,
with periodic and mixed outflow/slip-wall boundaries. A branched four-level,
six-node hierarchy includes adjacent children of different parents and refined
domain surfaces. It checks all-node CFL against independent local kernel calls,
composite inventory against the returned boundary ledger, positive temperature,
species closure, and exact state/temperature/diagnostic rollback after an
inadmissible large timestep.

An initial temperature-oracle failure exposed a different final EOS starting
estimate; final recovery now uses the first-stage estimate, matching serial.
An initial test incorrectly assumed that the fastest root cell was covered by
the finest patch; the CFL oracle now evaluates the actual node populations.
No tolerance was relaxed.

## Separate build failure

An initial unfiltered all-target build failed in the existing selected-AMR
checkpoint target: `amr_reactive_3d_checkpoint_mod.F90` imports
`h2o2_elementary_mechanism_mod`, which is unavailable on that target's module
search path. Log: `build/amr-tree-hydro-checks/all-target-build-failed.log`.
The focused targets subsequently built and passed. No full-build or selected
mechanism qualification is claimed; resolving that dependency is separate work.
