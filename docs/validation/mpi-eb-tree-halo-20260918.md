# MPI EB tree face halos — 2026-09-18

Implementation: `9c0d8a2cc4ecf8c5420081016d34bb15fba69637`, based on
`24ee433`. Module `mpi_eb_tree_halo_3d_mod` belongs to `pelef_mpi_support`.

`fill_mpi_eb_tree_halo_3d` builds spatial face halos for all patches owned by
the calling rank. U/T use the existing owner-local `amr_tree_state_3d` storage;
geometry remains replicated. It averages a private distributed snapshot, selects
same-level donors before successively coarser ancestors on the checked integer
lattice, and uses adjacent interior cells at physical outflow boundaries.
Batched sparse requests fetch only face-donor U/T through the existing cell
exchange. No full field gather is used.

Peer spacing, origins, apertures and tangential centroids are checked, including
closed target faces. The existing geometry `face_metrics` helper is now public
and shared with the serial overlay. Donor EOS recovery precedes a collective
commit of both halo containers. Inputs remain unchanged; any rank's failure
preserves all halo outputs. Owner-local output fields have spatial lower bounds
zero, with zero interiors, solid faces, edges and corners. Nonowners allocate
no patch values. Temporal interpolation and coupled MPI stepping are not yet
connected; selected-runtime target wiring and the MPI application remain open.

## Completed focused verification

Fresh Debug/Ninja build, GNU Fortran 13.3.0 `/usr/bin/gfortran`, OpenMPI
`/usr/bin/mpifort`, launcher `/usr/bin/mpiexec`:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/mpi-eb-tree-halo-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-halo-debug --target test_mpi_eb_tree_halo_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-halo-debug -R '^unit_mpi_eb_tree_halo_3d_[23]$' --output-on-failure --output-log build/mpi-eb-tree-halo-checks/ctest.log
python3 -m py_compile tools/check_mechanism_provenance.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 tests passed**, zero failures/skips, 0.55 s and 0.57 s, total 1.15 s.
A five-node, three-level cut-cell fixture includes same-level peers belonging
to different parents, immediate/remote ancestor donors and physical outflow.
Both two and three ranks exercise distributed owners and all patches on rank
zero (empty ranks). Every owned halo U/T array matches the serial API exactly,
including inactive zeros, and input fields remain unchanged. An invalid
unrefined root boundary cell passes averaging but fails EOS after exchange;
all output halos remain unchanged. Rank-inconsistent translated geometry is
also rejected without input/output mutation.

Configure/build completed without failures. Existing compiler warnings remain;
this is not a warning-free-build claim. Provenance, Python syntax, contract and
whitespace checks pass. This focused spatial library comparison is not a full
matrix, MPI coupled-application, performance or physical-validation result.

Logs: `build/mpi-eb-tree-halo-checks/{configure,build,ctest}.log`.
Detailed output: `build/mpi-eb-tree-halo-debug/Testing/Temporary/LastTest.log`.
