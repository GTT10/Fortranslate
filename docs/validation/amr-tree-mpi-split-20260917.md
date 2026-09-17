# MPI transport bounds and full operator composition — 2026-09-17

## Source and scope

Base: `8cc5735d9977e071a0f791d15c653c53a6430fe7`, local `main`.
Adds global transport/isothermal-wall timestep selection and the transactional
common-step R/2-T/2-H-T/2-R/2 composition. Existing chemistry input cloning is
shared without changing its numerical behavior.

SHA-256 at verification:

- `src/parallel/mpi_amr_patch_tree_transport_3d_mod.F90`:
  `4f6e6b070c8b13a9c628afc3018a4c4132d5d21575617f6fb6b79b9b68cc58da`
- `src/parallel/mpi_amr_patch_tree_chemistry_3d_mod.F90`:
  `589c4ef49c657242ab53fac25f62c55362a31b8f68da20171534b879b3fb788c`
- `src/parallel/mpi_amr_patch_tree_split_3d_mod.F90`:
  `4b879760780f28f2a83bea7d2868fd5b16a1a61a893c69ec393b6a51a98742a9`
- `tests/unit/test_mpi_amr_patch_tree_transport_3d.F90`:
  `049861ec05f9de3b2ec9c738d0e6e0865bd2547132cc0f1d4715c879511ffe1c`
- `tests/unit/test_mpi_amr_patch_tree_split_3d.F90`:
  `6d5337c10f13bb829958e0386451ee0c36770a437b56b37cea847ca61523edee`

## Completed checks

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6. Incremental build in the existing clean
system-MPI configuration; no compiler/backend option changes.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/tree-distribution-mpi-debug --target test_mpi_amr_patch_tree_transport_3d test_mpi_amr_patch_tree_split_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/tree-distribution-mpi-debug -R '^unit_mpi_amr_patch_tree_(transport|split)_3d_[23]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Builds succeeded. Focused CTest **4/4 passed**, 21.33 s total, no failures/skips.
Transport at 2/3 ranks took 4.64/4.54 s; full split took 6.06/6.06 s. Provenance,
contract and whitespace checks passed. No repository Python changed; syntax
checks were not applicable. Exact-control comparisons retain real-comparison warnings.

Transport fixtures now compare global dt and maximum diffusivity exactly to the
serial helper, including moving isothermal physical-wall cases. Invalid owner
transport data preserve both timestep outputs, in addition to the existing step
rollback assertions.

Full-split fixtures use a small three-level tree (root 2 cubed, ratios 2/3), seeded
elementary H2/O2, explicit/implicit chemistry, PCM/PLM flow, all transport terms,
periodic or moving-isothermal boundaries, and a disabled chemistry/transport mode
with inflow. Owner-local states/temperatures, total boundary budget, maximum reflux
correction and minimum theta match the serial common-step split exactly. Empty
ranks and all-on-one ownership are included. A hydro failure after preceding
reaction/transport work and inconsistent enable flags preserve every output.

Logs: `build/amr-tree-mpi-split-checks/{build,ctest}.log`.
Executables: `build/tree-distribution-mpi-debug/tests/test_mpi_amr_patch_tree_{transport,split}_3d`.

## Remaining scope

The caller still combines timestep limits and supplies replicated coefficient
records. Recursive level subcycling/coarse-time endpoints, rebalancing, selected
MPI tree runtime, application I/O/restart and broader qualification remain open.
These are bounded operator comparisons, not external solver parity, production
scaling or physical validation.
