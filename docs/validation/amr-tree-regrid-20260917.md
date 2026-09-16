# Scheduled regular 3D tree regridding — 2026-09-17

Source parent: `177cfb4caa908fea2c0deb3a39676c7aced8c60b` on `main`, plus this change.
PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Driver SHA256: `95d602895b8d139d141bfdd7a4eed851abf165fd236de389f4373a1f071b1ee8`.
Frontend SHA256: `91a5f0123d60563a49652185f8957c6d4658b5cd6511b9ab280a67e45a832879`.
Regression SHA256: `0ca6a37be90aa0a0793f46ec879549d3e6e30a3315dce2991381f60c4e9865dd`.

## Scope

The optional third CLI argument supplies strictly increasing positive event
times and topology paths. Relative paths resolve against the schedule directory.
Topologies are preloaded; timestep clipping reaches each exact event time.
Conservative state remapping, EOS recovery and a component inventory check
precede replacement of the active tree/state. The inventory scale uses absolute
cell contributions to handle cancelling momentum; tolerance is `5e-12`.
Initial inventories and accumulated boundary fluxes remain unchanged.
Final-time events occur before output. Inflow timestep bounds cover the finest
spacing of every scheduled topology; all other bounds use the active tree.

This is prescribed-time regridding, not automatic solution-based tagging.
The existing common timestep and PCM remap/ghost donor boundaries remain.
MPI, checkpoint/restart, selected frontend and scaling/physical qualification
remain separate work.

## Completed focused checks

GNU Fortran 13.3 Debug with bounds and FPE traps, OpenMPI 4.1.6, reused
`build/spray-mpi-migration-gnu-debug`, MPI/spray enabled. No toolchain change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_amr_tree_reactive_3d -j 4
ctest --test-dir build/spray-mpi-migration-gnu-debug \
  -R '^regression_amr_patch_tree_(application|regrid)_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Build and **2/2 focused regressions passed, 0 failures/skips, 10.13 s**.
The dynamic elementary reacting/transport/moving-wall case executes events at
`2.5e-10`, `5e-10`, and `1e-9 s`: moved ratio-3 refinement, complete coarsening,
and restoration of the original branched tree at final time. Maximum relative
remap discrepancy reported by the application: **1.5696e-16**.

Checks include exact logged event times, final exposed-cell volume and species
closure, output integration against the composite ledger, Euler boundary
balance through all remaps, a final-time identity remap preserving the solution
and exact budget file, and out-of-order schedule rejection preserving existing
outputs. The unchanged static application regression also passes, including its
regular-grid reference comparison. No tolerance/golden changes or test failures
occurred in this focused run.

Provenance/project contract, whitespace and Python AST parsing of both checkers
plus the new regression script passed (3/3). Build/CTest logs:
`build/amr-tree-regrid-checks/`. Per-run logs and CSVs:
`build/spray-mpi-migration-gnu-debug/tests/amr_tree_regrid_3d/`.
