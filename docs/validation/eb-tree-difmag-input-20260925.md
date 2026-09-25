# EB tree gas viscosity input — 2026-09-25

Code `31abf46`, parent `dfedb61`. Shared fixed/selected serial/MPI EB tree gas
application source accepts finite nonnegative `difmag` in `&eb_tree`, forwards
it to full composition and stores it in patch checkpoint configuration.
Positive coefficient changes therefore reject restart through the existing
patch solver identity. Omitted/zero keeps legacy behavior and identity.
MPI EB spray rejects positive input until its coupling forwards the coefficient.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Python 3.12.3, Ninja Debug; reused
`build/spray-amr-reactive-debug`. Fixed serial and MPI gas executables built.
Focused serial CTest **1/1 passed**, 26.94 s (26.98 s total), no failures or
skips. Compiler warnings remain. Selected variants and MPI runtime were not
independently checked in this increment.

The five-node/three-level planar fixture runs recursive stepping, chemistry
and molecular transport. Coefficient 0.037 gives exact continuous/restarted
state CSV equality; the interrupted CSV differs from the final reference.
Coefficient 0.05 on restart, -1 and NaN reject before output publication.
Omitted and explicit zero agree, including legacy-zero restart. This fixture
does not compare budget CSVs or positive versus zero dynamics. This is bounded
application/restart evidence, not full reference geometry or physical validation.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^reactive_eb_tree_3d_difmag_input$' --output-on-failure
python3 -m py_compile tools/check_eb_difmag_input.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

All listed checks passed. Artifacts: `build/eb-tree-difmag-input-checks/`
contains `build.log` and copied `LastTest.log`; the latter identifies the unique
fixture directory under `build/spray-amr-reactive-debug/tests/reactive_eb_tree_3d_difmag_input/`.
