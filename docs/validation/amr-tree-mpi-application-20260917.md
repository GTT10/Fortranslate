# Static MPI tree application — 2026-09-17

## Source and scope

Base: `b66b516fa4881672ca4ca2f96afca0b756591fb1`, local `main`.
Adds fixed-mechanism static MPI tree frontend/driver, owner-local analytic
initialization and stepping, common timestep selection/retry, and root final
CSV/budget output. Installs the executable when MPI is enabled.

SHA-256 at verification:

- `src/driver/mpi_amr_patch_tree_application_3d_mod.F90`:
  `49d4c3ad1580178908fc7ad79b4c9ce9b00cbd4a1b9f4581116f98801a578c62`
- `app/pelef_mpi_amr_tree_reactive_3d.F90`:
  `d1d961688a52b9c6a45d0154ec14991b7377cb55f36093c24c11a01b609c5469`
- `tools/check_mpi_amr_tree_application_3d.py`:
  `c86fab479146b2d0b497074495bb1ad909fb7c0cd46ec2215fabba515f118531`

## Completed checks

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6. Incremental build in the existing clean
system-MPI configuration; compiler/backend options unchanged.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/tree-distribution-mpi-debug --target pelef_mpi_amr_tree_reactive_3d pelef_amr_tree_reactive_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/tree-distribution-mpi-debug -R '^regression_mpi_amr_tree_application_3d$' --output-on-failure
python3 -m py_compile tools/check_mpi_amr_tree_application_3d.py
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Build succeeded, including the executable's MPI ABI 40 dependency audit. Final
focused CTest **1/1 passed**, 6.87 s total (test 6.84 s), no failures/skips. Python
syntax, provenance, project contract and whitespace checks passed. Integer-real
conversion and exact time-comparison compiler warnings remain.

The checker runs serial and 2/3-rank MPI on the existing elementary hotspot/branched
tree, then serial and 3-rank MPI on the root-only topology (two empty ranks).
Final exposed-cell CSV and component budget CSV are byte-identical. It also
checks unsupported checkpoint-group rejection and topology/output alias rejection,
with existing output/input contents preserved. These are five successful application
runs plus two expected rejections inside one CTest, not a configuration matrix.

The first checker run found that lexical-only comparison missed an absolute input
versus relative output alias. The driver now anchors both paths to the actual POSIX
working directory before comparison. Rebuilding and rerunning the same checker
passed; the failed run is not counted as successful qualification.

Logs: `build/amr-tree-mpi-app-checks/{build,ctest}.log`.
Per-run inputs, logs and CSVs:
`build/tree-distribution-mpi-debug/tests/mpi_amr_tree_application_3d/`.
Executable: `build/tree-distribution-mpi-debug/pelef_mpi_amr_tree_reactive_3d`.
The install rule was added; an installed-prefix audit was not run.

## Remaining scope

MPI restart/checkpoint, scheduled/automatic regridding, selected-mechanism runtime
and subcycling are not yet connected; the static frontend rejects extension groups.
Whole-node round-robin ownership and root initialization/output gathers have no
scaling qualification. Full-H2O2 loading is wired but was not run in this focused
application check. This is serial/MPI implementation parity, not external PeleC
parity or physical validation.
