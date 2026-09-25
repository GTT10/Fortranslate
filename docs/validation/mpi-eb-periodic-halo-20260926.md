# MPI EB tree periodic halos — 2026-09-26

Source `292ded9`, parent `6f20de9`, local main. MPI halo construction now
encodes periodic faces separately from prescribed inflow in its collective
control agreement. Paired faces and matching root metrics are required before
donor exchange. Coordinates wrap before peer/ancestor lookup; periodic faces
do not trigger prescribed inlet EOS conversion. Existing time interpolation
and private output publication remain in use.

Serial and MPI complete one-layer halos use target-local indices along
periodic axes when accessing embedded halo metrics; wrapped donor coordinates
can lie outside the target patch's metric array. Common and recursive MPI
hydro/transport still reject periodic time updates collectively until seam
flux reconciliation is implemented.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Ninja Debug, Python 3.12.3; existing
tests-enabled `build/spray-amr-reactive-debug`, no compiler/backend changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_halo_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_halo_3d_' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, two ranks 5.42 s, three ranks 4.97 s, total 10.43 s.
No failures or skips. Each test exercises distributed owners and all patches
on rank zero, leaving other ranks idle. New periodic checks compare all
owner-local gas and temperature halo values exactly against serial for face
donors, complete one-layer halos, six-layer halos and interpolated ancestor
states; volume fractions also agree. Rank-disagreeing periodic labels and
incompatible root metrics reject without changing halo outputs. Existing
prescribed inflow, invalid donor, width mismatch and temporal rollback checks
remain in the same test.

Provenance, Python syntax, project contract and whitespace checks passed.
Artifacts: `build/mpi-eb-periodic-halo-checks/build.log` and `LastTest.log`.
Mechanism SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

This is halo exchange evidence, not periodic time-update conservation or full
PeleC parity. Fine periodic seam metric reconciliation, application dispatch,
installed audit, full matrix, external solver parity and physical validation
were not qualified by this check.
