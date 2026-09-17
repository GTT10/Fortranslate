# Selected-mechanism static MPI tree application — 2026-09-17

## Scope and source

Base: `4c9b60c`, local `main`. Adds the configure-time selected MPI tree frontend,
isolated MPI tree runtime, bundle target/installation routing and elementary
selected regression fixture. Numerical MPI kernels are unchanged.

Verified source SHA-256:

- `app/pelef_mpi_amr_tree_reactive_3d_selected.F90.in`:
  `65ade30241cdebc54c6c548a0d2cadb074927134e7143de8e0c8443d5404df3d`
- `cmake/targets/SelectedRuntime.cmake`:
  `8fe48e18e59d254b923ac4a96a6dd45eed62ce73124869a01132babbfbbe6c7a`
- `tools/check_mpi_amr_tree_application_3d.py`:
  `c2fe6640b7161f34ec684681be235550ce6504e0ce045087212480adca6d6375`

## Completed checks

Clean build directory `build/selected-tree-mpi-debug`: GNU Fortran 13.3.0 Debug,
OpenMPI 4.1.6, Ninja, spray disabled. Configuration initially rejected a mixed
Intel launcher/system wrapper and the kinetics-only elementary JSON. Explicit
system launcher selection and the schema-1 full bundle corrected configuration;
the elementary test uses the existing pinned bundle composer. No source data or
provenance checks were bypassed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/selected-tree-mpi-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON -DPELEF_MECHANISM_BUNDLE=/home/gtt13/projects/pelec/mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 4 --target pelef_mpi_amr_tree_reactive_3d_selected test_selected_elementary_mpi_amr_tree_3d test_selected_elementary_amr_tree_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_mpi_amr_tree_selected_application_3d$' --output-on-failure
python3 -m py_compile tools/check_mpi_amr_tree_application_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Build and MPI ABI 40 dependency audits passed. CTest **1/1 passed**, 6.91 s total
(test 6.88 s), no failures/skips. The test comprises five successful application
runs and two expected rejections: branched serial/2-rank/3-rank, root-only
serial/3-rank (two empty ranks), unsupported checkpoint input and topology/output
alias. CSV and budget outputs match the selected serial application byte for
byte. Python syntax, provenance, project contract and whitespace checks passed.
Existing compiler warnings remain.

The full-H2O2 production target was built, not run. Its Ninja link command
contains only selected runtime/mechanism archives and MPI/system libraries;
`pelef_core` and fixed database archives are absent. An installed-prefix audit
was not performed. This check does not qualify every selected mechanism.

Logs: `build/amr-tree-selected-mpi-checks/{configure,build,ctest}.log`.
Run artifacts: `build/selected-tree-mpi-debug/tests/mpi_amr_tree_selected_application_3d/`.

## Boundaries

Static synchronous owner-local stepping only, with round-robin patch assignment
and root CSV gathers. MPI tree restart, regridding and subcycling remain open.
This is selected serial/MPI implementation parity, not external solver parity,
physical validation or scaling qualification.
