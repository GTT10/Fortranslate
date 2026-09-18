# MPI reacting EB spray with particle outflow — 2026-09-18

Implementation `ff548fa`, base `ff1a2e9`.

The MPI reactive EB spray wrapper accepts optional `particle_outflow`. True
requires films, component film masses, film settings and the shared ID watermark.
The effective flag participates in collective mode agreement. Each particle half
forwards private optional exterior storage to the MPI source interval; their
inventories are accumulated separately from gas boundary updates and added only
at whole-step publication. Wall and chemical ledgers keep their existing meaning.
The MPI application is not connected yet; specular EB outflow is also pending.

## Completed focused checks

GNU Fortran 13.3.0/Open MPI 4.1.6, Debug/Ninja, reused fixed tests-enabled and
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

CTest **2/2 passed**, 0 failures/skips, 5.94 s total using configured
`/usr/bin/mpiexec`. Builds, provenance, syntax, contract and diff checks passed.
Existing impact, film, KH, splitting and rollback cases remain.

The new two-component fixture starts one parcel on the last rank near the
embedded y wall and exterior z face, velocity (20,-200,5) m/s. A 10 ns reactive
step creates four splash children with watermark 5; some escape and some remain.
Reaction, molecular transport, hydro and film coupling are enabled. Particle drag,
heat and evaporation are disabled in this added fixture. Composite gas, particle
and film inventory closes against wall, exterior and chemical change within the
existing relative 3e-12 tolerance and 1e-30 floor.

With an invalid gas solver, the first particle half remains private: parcel and
component arrays, empty films, watermark, owner-local gas and temperatures, all
ledger sentinels and limiter diagnostic are preserved. No external-solver or
physical validation was run.

Artifacts: `build/spray-mpi-eb-outflow-reactive-checks/{provenance.log,build.log,ctest.log,build-selected.log}`.
