# MPI EB geometry context agreement and restart — 2026-09-25

Code `f534f5d`, parent `627a184`. Replicated EB geometry validation now
compares the optional fine-resolution halo's presence and every value across
ranks. A fixed-size record prevents mismatched allocation presence from creating
unequal collective message lengths. Existing local validity checks run first.
This validation is shared by averaging, stepping, integration and checkpoints.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Python 3.12.3, reused Ninja Debug
`build/spray-amr-reactive-debug`. **2/2 distinct CTests passed**, no failures or
skips: three-rank positive-viscosity EB tree 1.82 s (1.87 s total), EB tree
application 23.97 s (24.01 s total). Compiler warnings remain. Mechanism
provenance, Python syntax, project contract and whitespace checks passed.

The unit test changes a valid halo value on rank zero and separately removes
that rank's optional halo. Both cases reject collectively and preserve owner
states, temperatures and ledgers under distributed/single-owner layouts,
four/five-node trees and x/y/z cuts. Existing serial/MPI positive-step parity
and other failure cases remain.

The application checks serial/2/3-rank state output byte equality, 2-to-3-rank
restart, both serial/MPI restart directions, exact same-rank budgets and identical
serial/MPI checkpoint patch records. All five saved patches contain
`GEOMETRY_HALO_V1`. Its conservation residual is 3.18097960520324641e-15.
Existing corrupted-record/context/output-preservation checks pass. This test
uses the standard fixed-step common-level fixture; recursive MPI, regridding,
selected mechanisms, physical validation and external solver parity were not
rerun here.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_viscosity_3d_3$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^regression_mpi_eb_tree_application_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Logs: `build/mpi-eb-geometry-agreement-checks/{build.log,application-build.log,LastTest-unit.log,LastTest-application.log}`.
Application outputs: `build/spray-amr-reactive-debug/tests/mpi_eb_tree_application_3d/`.
