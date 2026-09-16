# Conservative regular 3D patch-tree reflux — 2026-09-17

Parent: `9c617e7`, local `main`. Contract and algorithm:
[patch-tree design](../design_decisions/amr-patch-tree-3d.md).

Added a transactional flux-synchronization/reflux API for multiple nodes and
arbitrary depth. Time-integrated child fluxes are area-averaged bottom-up.
Shared faces prefer finer replacement values and otherwise require consistency;
this includes cross-parent peers and periodic seams. Only uncovered cells get
flux-divergence corrections, and final restriction updates covered parents.
Physical-domain flux replacement is returned as a root-ledger correction.
The checked integer layout is now shared with the existing ghost builder.

This API consumes fluxes over one common physical interval. It does not produce
hydro fluxes, schedule subcycles, recover temperatures or assert EOS admissibility
of corrected states. General tree evolution, chemistry/transport, dynamic remap,
MPI/EB and application integration remain open. No scaling or physical validation
is implied.

## Reproduction and results

Reused GNU Fortran 13.3 Debug configuration with bounds/FPE checks:
`build/spray-mpi-migration-gnu-debug`. No compiler/MPI/backend change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_amr_patch_tree_reflux_3d test_amr_patch_tree_ghost_3d test_amr_patch_tree_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^unit_amr_patch_tree(_ghost|_reflux)?_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Final focused tests: **3/3 passed, 0 failures, 0 skips**, 0.07 s.
The new test builds four levels/six boxes, adjacent children on separate parent
branches and a deeper child touching their shared interface and a domain surface.
It applies analytically specified, level-dependent integrated flux divergences
to every node before reflux. Under physical boundaries, corrected composite
inventory equals an independent sum over exposed finest leaf faces; the returned
ledger correction equals that sum minus the original root ledger. Under periodic
boundaries the composite inventory is conserved. Restriction propagates the
corrected inventory to the root. Tests use two independent conserved components.
Conflicting shared peer flux and NaN flux reject without changing any level or
the output diagnostics. Existing tree inventory and temporal/physical ghost tests
also pass after sharing integer layout construction.

Build, provenance preflight, project contract, Python AST syntax for two
verification scripts, and whitespace checks passed. No full matrix, installed
binary audit, external solver comparison or experimental validation was run.
Logs and hashes: `build/amr-tree-reflux-3d-checks/`.

```text
20f92ec13a3196e9e78cfd005dee54a8764a2d998f6fdad679e8b7d89a9a0459  src/amr/amr_patch_tree_reflux_3d_mod.F90
0fb3feb569e1c2f09cea1079802c69407a3b8be78546d963025a807b9451823f  src/amr/amr_patch_tree_3d_mod.F90
8db8eae03f4f9f8b989e5e9860437ffe227ae7f1325ace911cf09e18f4ca7e95  src/amr/amr_patch_tree_ghost_3d_mod.F90
51ca490583f2ad19692f46b6a157a558a6b2ce428f7a0496457e9107dc358c78  tests/unit/test_amr_patch_tree_reflux_3d.F90
```
