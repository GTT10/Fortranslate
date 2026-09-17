# Two-level EB application — 2026-09-17

Implementation `528df042ffa4c3555196d77d0d0a0afa71e35e1b`, based on
`acd8be1` on local `main`.

## Delivered behavior

`pelef_reactive_eb_two_level_3d` reads a dedicated `eb_two_level` namelist,
constructs finest geometry, restricts it to the coarse domain, extracts one
fine patch and advances the existing combined reaction/transport/hydro library.
The supplied case completes three explicit 1e-8 s intervals. The executable is
listed for installation; an installed-binary audit was not performed here.

Outputs contain only composite active cells, with level/local indices, fluid
centroids, volume fractions, temperature and every conserved component. A separate
budget CSV records initial/final inventories and accumulated outer-boundary,
wall and chemical changes. Output occurs after the inventory gate passes.
Invalid refinement and lexical input/output aliases reject before output writing.
The fixed timestep is user-controlled; no automatic CFL bound or retry is claimed.

The frontend supports existing plane/extruded/nodal geometry builders, but the
new application regression exercises a plane. It uses fixed elementary H2/O2
initialization and chemistry. Checkpoint/restart, selected-mechanism frontend,
MPI and multi-patch/regridding remain open. The executable does not accept those
unsupported controls.

## Completed checks

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6, Ninja; existing matching build
`build/selected-tree-mpi-debug`, MPI on, spray off. No MPI execution or full
configuration/install matrix was run.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 2 --target pelef_reactive_eb_two_level_3d
python3 -m py_compile tools/check_reactive_eb_two_level_3d.py
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_reactive_eb_two_level_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **1/1 PASS, 0 failures, 0 skips**, test 0.69 s, total 0.72 s.
The application completes three steps to 3e-8 s, reports inventory error
`4.6057075715980637e-16` and minimum species-limiter factor one.
The checker verifies 104 composite cells, unique level/cell keys, no duplicated
refined parents, finite positive density/temperature, species closure and
initial mass from geometry/default densities. It independently integrates CSV
conserved fields and matches the budget's final inventories and residuals.
Invalid refinement preserves existing outputs; an output/input alias preserves
the input file.

The first regression found a trailing comma in composite rows. A terminating
colon in the unlimited Fortran format removed the extra field; the final run
passed without relaxing checks. A diagnostic-message initialization and comment
clarification followed the run; numerical and output formatting code was unchanged.
Provenance, Python syntax and source contract passed. Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

Logs: `build/eb-two-level-app-checks/{build,final-build,csv-build,ctest,final-ctest}.log`.
Artifacts: `build/selected-tree-mpi-debug/tests/reactive_eb_two_level_3d/`.
This is a focused application test, not external PeleC parity, convergence or
physical validation.
