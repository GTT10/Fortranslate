# MPI EB tree restriction and inventory — 2026-09-18

Implementation: `49b6333adca377f71c3f1b75347836e5de9631a0`, based on
`cefc9d8`. Module `mpi_eb_tree_transfer_3d_mod` is part of `pelef_mpi_support`.

## Implemented contract

- `validate_mpi_eb_tree_geometry_3d`: collective agreement on topology, owners,
  component count and all replicated EB geometry metrics/classifications;
  verifies each child coarsens to its parent section.
- `average_down_mpi_eb_tree_3d`: owner-local `amr_tree_state_3d` containers hold
  U and one-component temperature. The child owner applies the existing
  volume-fraction-weighted restriction/EOS kernel, then sends only coarse U/T
  to the parent owner. Reverse traversal handles remote parents at successive
  depths. Solid parent payloads are preserved. All state remains private until
  every child succeeds collectively; failures preserve both caller containers.
- `integrate_mpi_eb_tree_3d`: owners integrate exposed fluid cells, excluding
  cells covered by immediate children, then sum component inventories across
  the borrowed communicator. Invalid exposed values preserve output.

Existing regular-tree scatter/gather can distribute these containers. Numerical
restriction does not gather full fields to a root. Geometry remains replicated
and is checked on every call; this is not a scaling qualification. Species
thermodynamic data are a caller-supplied common mechanism contract. MPI halo
exchange, reaction/hydro/transport, selected-target wiring and an EB tree MPI
application remain separate work.

## Completed focused verification

Fresh Debug/Ninja build, GNU Fortran 13.3.0 `/usr/bin/gfortran`, OpenMPI
`/usr/bin/mpifort` and `/usr/bin/mpiexec`:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/mpi-eb-tree-transfer-openmpi-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-transfer-openmpi-debug --target test_mpi_eb_tree_transfer_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-transfer-openmpi-debug -R '^unit_mpi_eb_tree_transfer_3d_[23]$' --output-on-failure --output-log build/mpi-eb-tree-transfer-checks/ctest.log
python3 -m py_compile tools/check_project_contract.py tools/check_mechanism_provenance.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 tests passed**, zero failures/skips, 0.51 s each, 1.05 s total. Both rank
counts exercise a five-node, three-level cut-cell tree with distributed owners
and with all patches on rank zero (empty ranks). Gathered U and T match serial
restriction exactly, including inactive solid payloads. Maximum normalized
composite-inventory difference is **1.9809508747161377e-15**, below 3e-14.
A late invalid child EOS preserves U/T even after another private restriction.
Locally valid but rank-inconsistent translated geometry is rejected collectively;
NaN in an exposed cell preserves sentinel inventory outputs.

An initial configure in `build/mpi-eb-tree-transfer-debug` selected the Intel
MPI wrapper from the environment and failed the existing mpi_f08 compatibility
check. A separate fresh directory with explicit OpenMPI wrappers succeeded.
No test failed. Compilation retains existing warning diagnostics; this is not
claimed as a warning-free build. Provenance, Python syntax, project-contract and
whitespace checks passed. No full matrix, installed audit or physical validation
was run.

Logs: `build/mpi-eb-tree-transfer-checks/{initial-configure,configure,build,ctest}.log`.
Detailed test output:
`build/mpi-eb-tree-transfer-openmpi-debug/Testing/Temporary/LastTest.log`.
