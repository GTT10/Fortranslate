# MPI composite EB particle sampling — 2026-09-18

Implementation: `3c7a8a7936a2a51614ddb76683838719218e309a`, based on `cff088a`.

Paired optional EB geometry/nodal arguments select composite EB stencils in
`sample_mpi_amr_tree_spray_field`. The existing request routing fetches only
required cell values from gas owners; fields are not globally gathered.
Presence and nonperiodic flags are collective. All ranks validate replicated
nodal values, patch bounds, fluid fractions and centroids, including ranks with
zero queries. Samples publish after collective success; regular AMR calls keep
their existing stencil path.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; reused unchanged fixed
tests-enabled and selected H2/O2 tests-disabled configurations.
Compiler `/usr/bin/gfortran`, MPI wrapper `/usr/bin/mpifort`, launcher
`/usr/bin/mpiexec`.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_sample --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sample_[24]$' --output-on-failure --output-log build/spray-mpi-eb-sample-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_gas --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_gas_2$' --output-on-failure --output-log build/spray-mpi-eb-sample-checks/ctest-regression.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/check_project_contract.py tools/check_mechanism_provenance.py
python3 tools/check_project_contract.py --source .
git diff --check
```

New MPI checks **2/2 passed**, two ranks **0.48 s**, four ranks **0.49 s**, total
**1.01 s**. Existing two-rank regular AMR gas exchange regression **1/1 passed**,
**0.55 s**, total **0.59 s**. Selected compilation, provenance, syntax, contract
and whitespace checks passed. No failures or skips; existing compiler warnings
remain. Selected runtime execution was not tested.

The fixture uses a 4-cubed root, ratio-two child and y=1.25 cut. Root and fine
fields have distinct constant values 10 and 20 and different MPI owners. Rank
zero queries (1.9,1.3,1.5), requiring both local and remote support; other ranks
have no queries. The sample matches the serial composite-stencil weighted
reference within 1e-13, between 10 and 20. A rank-local nodal perturbation, solid
query and zero request capacity reject without changing old samples. Globally
empty queries also succeed. This checks routing/composition, not an independent
interpolation-accuracy oracle.

Artifacts: `build/spray-mpi-eb-sample-checks/` contains `build.log`,
`build-selected.log`, `build-regression.log`, `ctest.log`, `ctest-regression.log`.

EB deposition, fluid-volume parent restriction, masked EOS recovery and coupled
particle evolution remain open. Geometry is replicated and checked in full,
bounded to one million finest cells; no scalable geometry-distribution or
production performance claim is made. No convergence, external solver parity
or physical validation was performed.
