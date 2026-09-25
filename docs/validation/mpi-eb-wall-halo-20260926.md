# MPI EB tree wall halos — 2026-09-26

Implementation 021e360, parent 91c6210, local main, version 0.246.0.
MPI halo controls now collectively distinguish slip_wall and no_slip_wall.
Physical indices reflect before peer/ancestor request routing; returned
momentum receives accumulated reflection signs after EOS recovery and inlet
handling. Full-halo volume lookup retains original exterior coordinates for
explicit wall geometry. Candidate outputs still commit only after collective
success. Rank-disagreeing wall labels reject before donor exchange.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3, Ninja Debug, existing
fixed tests-enabled build/spray-amr-reactive-debug. No backend changes.
The two- and three-rank halo tests passed **2/2**, 7.02 and 6.30 seconds
(CTest reported total 13.37 seconds), no failures or skips.

The fixture now uses nonzero x/y/z velocity. Distributed and all-on-rank-zero
ownership cover empty ranks. Mixed x-slip/y-no-slip walls match serial state,
temperature and volume arrays exactly for face-only, complete one-layer and
six-layer halos, including corners and ancestor time interpolation. Existing
outflow, inflow, periodic and rollback cases remain in the same tests.
Rank-local label mismatch leaves prior halo outputs unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_halo_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_halo_3d' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, project contract and whitespace passed. Existing compiler
warnings remain. No tolerances or golden data changed.
Artifacts: build/mpi-eb-wall-halo-checks/build.log and LastTest.log.
Mechanism SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
Pinned PeleC: bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8.
Pinned PelePhysics: 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.

## Limits

This qualifies MPI adiabatic wall halo construction. Coupled tree wall stepping,
reflected extended geometry, application input and isothermal tree walls remain
separate work. No full matrix, installed audit, external-solver parity or physical
validation ran. Complete PeleC parity remains unproven.
