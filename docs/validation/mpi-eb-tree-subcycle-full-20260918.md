# MPI EB recursive transport and reacting coupling — 2026-09-18

Implementation: `8479f5d02ab5551e5c99dcd86e4144a7c49d9659`, based on
`6ed7b7a`. Adds `advance_mpi_eb_tree_subcycle_transport_3d` and optional
`subcycling` dispatch in `advance_mpi_eb_tree_full_3d`.

## Implemented behavior

Molecular transport recursively advances refinement-ratio-squared child
intervals, using time-interpolated ancestor halos, closed physical diffusion
faces and the existing positivity-limited patch kernel. Peer corrections define
provisional endpoints before children advance. Accepted effective fluxes
accumulate and synchronize over each parent interval. Ratio squaring is checked
for integer overflow. The final child step uses the remaining parent interval.
The effective root flux supplies the boundary ledger; MPI minimum reduction
supplies the minimum accepted limiter across all levels and substeps.

The coupled API selects recursive hydro and transport inside R/2 T/2 H T/2 R/2.
Absent or false `subcycling` retains common steps. The mode must agree across
ranks before branching. Chemistry/transport switches remain independent;
explicit/implicit reaction uses the existing exposed-cell operator. Owner-local
U/T, boundary/wall/chemical ledgers and minimum limiter publish only after all
stages succeed. Geometry remains replicated; no full state gather is introduced.

Application dispatch, checkpoint mode integration and scaled MPI automatic
intervals remain separate work.

## Completed focused verification

Fresh GNU Fortran 13.3.0 Debug/Ninja build, `/usr/bin/gfortran`, Open MPI wrapper
`/usr/bin/mpifort` and launcher `/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/mpi-eb-tree-subcycle-full-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-subcycle-full-debug --target test_mpi_eb_tree_subcycle_transport_3d test_mpi_eb_tree_subcycle_full_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-subcycle-full-debug -R '^unit_mpi_eb_tree_subcycle_(transport|full)_3d_[23]$' --output-on-failure --output-log build/mpi-eb-tree-subcycle-full-checks/ctest.log
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**3/3 tests passed**, no failures/skips, total **75.09 s**:

- Two/three-rank transport: **13.58/13.61 s**. Three cut-plane directions,
  four/five-node trees and distributed/single-owner layouts produce **24** MPI
  inventory comparisons, maximum normalized residual **6.798092530336388e-15**.
  U/T, boundary ledger and limiter match serial exactly. Negative interval,
  rank-disagreeing interval/flags and invalid state reject with rollback.
- Two-rank coupled step: **47.86 s**. Four explicit chemistry/transport switch
  combinations and one implicit fully coupled single-owner case (with an empty
  rank) match serial U/T, boundary and limiter exactly. Wall/chemical ledgers
  meet the unchanged 3e-14 inventory-scaled bound. The five MPI inventory checks
  have maximum normalized residual **2.4649337620366713e-15**. Invalid hydro
  after successful reaction/transport preserves all original fields and ledgers.
  Rank-disagreeing chemistry and subcycling switches are rejected collectively.

Provenance, configure/build, Python syntax, project contract and whitespace
checks pass. No failed build/test attempt occurred; compiler warnings remain.
No Python source changed. These focused API checks do not qualify application
subcycling, temporal convergence, selected MPI, installed binaries, scaling,
a full configuration matrix or physical validation.

Logs: `build/mpi-eb-tree-subcycle-full-checks/{configure,build,ctest}.log`.
