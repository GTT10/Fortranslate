# MPI tree ownership and transfers — 2026-09-17

## Source and scope

Base: `8e65ef948f9088aba90bd0d95199d36a20ac2fa8`, local `main`.
The new distribution module scatters/gathers tree fields between an arbitrary
root and explicit zero-based patch owners. Payload arrays exist only on owners;
metadata is replicated. Collective checks reject inconsistent topology, roots,
owner maps, components, shapes or allocations before payload transfer.

SHA-256 at verification:

- `src/parallel/mpi_amr_patch_tree_distribution_3d_mod.F90`:
  `269a29d485ded212d76b4118b346212a23cde4758558045cfd9b6a268e5dfe39`
- `tests/unit/test_mpi_amr_patch_tree_distribution_3d.F90`:
  `d7993685ac3180c0aad28a1b9dc7fa4025b24bb3fd15b4fae5c91485f8e071d6`

## Completed checks

Fresh GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6 build with MPI enabled.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/tree-distribution-mpi-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DPELEF_ENABLE_MPI=ON -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DMPI_mpi_LIBRARY=/usr/lib/x86_64-linux-gnu/libmpi.so -DMPI_hwloc_LIBRARY=/usr/lib/x86_64-linux-gnu/libhwloc.so
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/tree-distribution-mpi-debug --target test_mpi_amr_patch_tree_distribution_3d -j 6
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/tree-distribution-mpi-debug -R '^unit_mpi_amr_patch_tree_distribution_3d_[23]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

- Clean configure/build succeeded. Focused MPI CTest: **2/2 passed**, 0.92 s
  total (2 ranks 0.44 s, 3 ranks 0.45 s); no failures or skips.
- `ldd` under the same environment resolved MPI and hwloc to system libraries.
- Provenance, project contract and whitespace passed. No Python repository source
  changed; Python syntax verification was not applicable.
- Both tests use three levels with ratios 2/3, a nonzero root rank, spatially varying
  two-component fields, distributed and all-on-one-rank ownership, and empty ranks.
  They verify exact cell/component ordering and gathered owner modifications,
  allocation only on owners, unchanged nonroot gather destinations, collective
  rejection of different valid topologies/owner maps/root choices, and missing
  owner data with destination preservation.

Preliminary tests in the existing MPI build also passed. Inspection revealed its
inherited environment loaded Intel's hwloc; rerunning without LD_LIBRARY_PATH and
LIBRARY_PATH resolved system hwloc and passed. The final evidence above comes
from a fresh build with the cleaned runtime environment, not the mixed-library
invocation. These are repeated focused checks, not additional matrix coverage.

Logs: `build/amr-tree-mpi-distribution-checks/{configure,build,ctest}.log`.
Executable: `build/tree-distribution-mpi-debug/tests/test_mpi_amr_patch_tree_distribution_3d`.

## Remaining scope

This is ownership/transfer infrastructure, not distributed tree evolution. Full
fields are gathered only when explicitly requested; no root-gather time stepping
is implemented. Partial local fields do not meet serial full-tree validity.
Ownership planning, halos, restriction/reflux, stepping, rebalancing and MPI tree
application integration remain open. Whole-node ownership cannot divide a large
root patch, and no performance/scaling result is claimed. MPI communication/runtime
failures remain subject to the communicator error handler; rollback claims here
cover collectively detected input/allocation failures. No external solver parity,
physical validation or selected-MPI application qualification is claimed.
