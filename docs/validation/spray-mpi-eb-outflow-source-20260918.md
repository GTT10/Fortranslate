# MPI EB source coupling with exterior inventories — 2026-09-18

Implementation `10d5b81`, base `8cbcf56`.

`advance_mpi_amr_mixture_coupling` accepts optional `eb_outflow_change(nvar)`
with nodal EB geometry and impact films. The optional presence participates in
collective mode agreement. Each trial forwards a private component-order
inventory to the MPI impact interval. Accepted substeps map components into
gas vapor-species slots and accumulate already-global values without another
MPI sum. Existing `boundary_change` remains the EB system wall change.

KH eligibility clips its candidate trajectory to the first root crossing before
querying the level set. Final particle and film migrations remain inside the
whole-interval transaction; both ledgers publish only on complete success.
The API requires impact films for this mode. Specular EB exterior flight,
reactive wrapper and application connection remain future work.

## Completed focused checks

GNU Fortran 13.3.0/Open MPI 4.1.6, Debug/Ninja, reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Selected library compiled only.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_impact_step --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_impact_step_(2|4)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **2/2 passed**, 0 failures/skips, 1.71 s total with configured
`/usr/bin/mpiexec`. Builds, provenance, syntax, contract and diff checks passed.
Existing deposition, transfer, KH release, ownership and rollback checks remain.

The new fixture starts one mixed-liquid parcel on the last rank with velocity
(20,-200,5) m/s and 0.25/0.75 component masses. Drag, heat, evaporation and KH
are disabled for this flight fixture. In a 1 ms source interval it splashes and
all four children escape, leaving watermark 5. Composite gas/particle/film
inventory closes against wall plus exterior change at existing relative 3e-12
tolerance and 1e-30 floor. Each mapped exterior species follows its original
mass fraction to 1e-23 kg. A one-substep cap rejects the interval, preserving
particle/component/film state, watermark and both ledger sentinels.

Artifacts: `build/spray-mpi-eb-outflow-source-checks/{provenance.log,build.log,ctest.log,build-selected.log}`.
No physical validation, external-solver parity, or MPI application run is claimed.
