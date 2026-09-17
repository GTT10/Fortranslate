# MPI EB temperature adaptation — 2026-09-18

Implementation: `f5c45f3f52cc6f6564535cfe2db6e5c7e5b98c51`, based on
`a28c03f`. The fixed-mechanism MPI EB tree application now connects temperature
tagging to distributed remapping at scheduled adaptation events.

## Implemented behavior

`tag_mpi_eb_tree_temperature_3d` reads owner-local temperatures and replicated
immutable geometry. Only exposed fluid cells vote; covered parents and solid
cells are skipped before temperature validation. Fine hot spots are retained
without first averaging them into parents. Hot-cell integer boxes and temperatures
are gathered across ranks and passed to the existing deterministic tile builder.
Solid-only nonroot tiles are pruned, and topology is published only after
collective success. The existing capacity bound applies before solid pruning.

The API checks shared refinement ratios, old geometry against the coarsened
finest geometry, controls across ranks and final topology agreement. It preserves
the previous output topology on rejection. The application reuses owner-local
remapping and the serial V4 adaptive checkpoint envelope. CSV/checkpoint output
still gathers state on the root; initial setup/restart temporarily replicates it.

## Completed focused verification

Fresh GNU Fortran 13.3.0 Debug/Ninja build, `/usr/bin/gfortran`,
`/usr/bin/mpifort` and `/usr/bin/mpiexec` (Open MPI).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/mpi-eb-tree-adaptive-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-adaptive-debug --target test_mpi_eb_tree_tagging_3d pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-adaptive-debug -R '^(unit_mpi_eb_tree_tagging_3d_[23]|regression_mpi_eb_tree_adaptive_3d)$' --output-on-failure --output-log build/mpi-eb-tree-adaptive-checks/ctest.log
python3 -m py_compile tools/check_mpi_eb_tree_application_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**3/3 tests passed**, no failures/skips, total **36.27 s**. Two/three-rank
unit tests took 0.91/0.97 s; the application regression took 34.36 s.
The unit fixture checks exact serial topology for distributed and single-owner
layouts (including empty ranks), fine hot spots, depth-specific thresholds,
root-only coarsening, solid-tile pruning, capacity rejection and exposed NaN
rejection. Inactive NaNs are ignored; rejected calls preserve the prior topology.

The fixed explicit application adapts at 1.5e-9 and 2.5e-9 s. Serial/two/three-rank
runs produce identical final state CSVs. Changed-rank restart before/after each
event and both serial/MPI restart directions also match. Same-rank restart
matches budget bytes; changed-rank budgets meet the unchanged 3e-14 scaled
inventory bound. Independent volume integrals and ledger checks pass with
maximum normalized residual **3.31928290111239962e-15**. A cold-threshold run
coarsens to the root and matches serial exactly; capacity exhaustion is rejected.
Changed tagging policy and existing corrupt-record/interval-policy failures
preserve sentinel outputs.

Configure/build, provenance, Python syntax, project contract and whitespace
checks pass. Compiler warnings remain. This focused check does not qualify
automatic-plus-adaptive integration independently, MPI subcycling, selected MPI,
installed binaries, scaling, a full configuration matrix or physical validation.

Logs: `build/mpi-eb-tree-adaptive-checks/{configure,build,ctest}.log`.
Runtime artifacts: `build/mpi-eb-tree-adaptive-debug/tests/mpi_eb_tree_adaptive_3d/`.
