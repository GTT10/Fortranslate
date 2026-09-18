# MPI EB external acceleration — 2026-09-18

Implementation `07c41d8`, base `ba57bea`.

The MPI EB reactive wrapper accepts uniform acceleration and a separate source
inventory. Collective finite/range and exact rank-agreement checks precede force
work. Existing geometry and owner-local state validators guard the half kicks.
Gas momentum/energy changes use the same midpoint-work convention as serial;
local source contributions use EB fluid volume and exclude cells covered by
children. Covered coarse fluid state still receives the kick for consistent
donors. Live parcels contribute represented-mass impulse/work; attached films
remain stationary. A collective success check precedes MPI_SUM of the local
source inventory. Both kicks and the source output remain inside full rollback.

## Completed focused checks

GNU Fortran 13.3.0/Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Selected library compiled only.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_impact_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_impact_reactive_(2|4)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **2/2 passed**, 0 failures/skips, 6.75 s total, configured
`/usr/bin/mpiexec`. Builds, provenance, syntax, contract and diff checks passed.

The added nested cut-cell fixture applies (1000,-2000,3000) m/s² to gas and one
mixed parcel. Composite inventory closes against boundary, wall, chemistry and
force at existing relative 3e-12 tolerance with 1e-30 floor. Force has nonzero
momentum and zero mass/species entries. A rank-0 acceleration mismatch rejects
collectively and preserves the source sentinel. Invalid hydro after the first
kick preserves owner-local gas/temperature, parcel velocity, watermark and source
sentinel. Existing impact/film, KH, splitting, outflow and rollback cases remain.

Artifacts: `build/spray-mpi-eb-force-checks/{provenance.log,build.log,ctest.log,build-selected.log}`.
MPI application/checkpoint exposure is not yet connected. No physical or
external-solver validation was run.
