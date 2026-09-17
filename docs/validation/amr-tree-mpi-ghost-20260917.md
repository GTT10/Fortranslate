# Synchronous MPI tree halos — 2026-09-17

## Source and scope

Base: `8366e8d235044ecb62b7587d38100b59cef7cf4e`, local `main`.
The shared reader-backed ghost sampler connects sparse owner-cell requests to
synchronous MPI PCM/PLM boundary filling. Each owner receives only its own padded
outputs. Serial validation and temporal endpoint behavior remain available.

SHA-256 at verification:

- `src/amr/amr_patch_tree_ghost_3d_mod.F90`:
  `c5a90a5253850f1e4637f045e9c6a36f128ee8836e3fd6013ec7720361ad2dd1`
- `src/parallel/mpi_amr_patch_tree_ghost_3d_mod.F90`:
  `9380755d7bc64e48f89c93dcac2768dc28b10ad97c1398b300fcc090b08e2083`
- `tests/unit/test_mpi_amr_patch_tree_ghost_3d.F90`:
  `6b7baa877a1e395c932d43493e6bb58ee6a419397d12938deb4b60cf8ed1c2ec`

## Completed checks

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6. Reused the clean system-MPI configuration
from [distribution checks](amr-tree-mpi-distribution-20260917.md); no toolchain or
optional-backend configuration changed.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/tree-distribution-mpi-debug --target test_mpi_amr_patch_tree_ghost_3d test_amr_patch_tree_ghost_3d test_amr_patch_tree_plm_ghost_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/tree-distribution-mpi-debug -R '^unit_(mpi_amr_patch_tree_ghost_3d_[23]|amr_patch_tree_(plm_)?ghost_3d)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
readelf -W -l build/tree-distribution-mpi-debug/tests/test_mpi_amr_patch_tree_ghost_3d
```

- Build succeeded. Final focused CTest **4/4 passed**, 1.46 s; no failures/skips.
  MPI 2/3 ranks took 0.67/0.71 s; serial ghost/PLM tests took 0.03/0.02 s.
- MPI fixtures use four nodes and three levels with ratios 2/3, distributed and
  all-on-one ownership, empty ranks, periodic and mixed physical faces, moving
  viscous walls, PCM/MC/minmod and two species. Padded results match serial
  exactly. Missing owner outputs remain unallocated. Nonfinite donors, invalid
  boundary options on empty ranks and post-exchange negative-density moving-wall
  failures reject collectively and preserve every prior output.
- Existing serial tests retain independent dense/affine oracles, coarse-time
  endpoints and rejection coverage after the sampler extraction.
- Provenance, project contract and whitespace checks passed. No repository Python
  changes; Python syntax checks were not applicable.
- Reader type-bound methods avoid internal-procedure trampolines; the final MPI
  executable has a read/write, non-executable GNU_STACK segment.
- Build warnings include real equality/conversion and existing optional-array
  descriptor warnings in dependent modules; the build is not warning-free.

Logs: `build/amr-tree-mpi-ghost-checks/{build,ctest}.log`.
Executable: `build/tree-distribution-mpi-debug/tests/test_mpi_amr_patch_tree_ghost_3d`.

## Limits

MPI ghost filling is synchronous. MPI coarse-time endpoints, distributed level
stepping, restriction/reflux, rebalancing and tree application integration remain.
Planning repeats geometry walks and retains duplicate stencil requests, including
local interiors; no communication/performance scaling claim is made. Metadata is
replicated and ownership is whole-node. No selected-MPI tree application, external
solver parity or physical validation is established by these focused tests.
