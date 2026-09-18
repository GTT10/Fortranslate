# MPI specular EB outflow application and restart — 2026-09-18

Implementation `c8e52fe`, base `4607f6b`.

The MPI reactive wrapper accepts outflow without films, forwards optional
exterior storage in both specular source halves, and accumulates it into the
final boundary inventory. Collective mode checks still enforce agreement and
valid film argument groups. The application no longer requires impact films for
`particle_outflow`. Existing checkpoint ID-watermark rules already accommodate
removed particles, so no checkpoint validation or format changes were needed.

## Completed focused checks

GNU Fortran 13.3.0/Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Selected executable compiled only.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d test_spray_mpi_eb_impact_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_mpi_eb_specular_outflow_application|spray_mpi_eb_impact_reactive_(2|4))$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
python3 -m py_compile tests/application/check_mpi_eb_impact.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **3/3 passed**, 0 failures/skips, 26.15 s total. Reactive regression tests
passed at 2/4 ranks (2.87/2.85 s); the new application test passed in 20.38 s.
Builds, provenance, syntax, contract and diff checks passed.

The application fixture injects two water parcels near x=0.004 m moving outward.
Both escape, leaving no particles or films. Four launches compare 4-rank full
execution against 2-rank checkpoint and 4-rank resume, then reject changed wall
settings through context matching. Three regrids and rebalancing remain enabled.
Gas, particle, film and budget rows match within the existing relative 3e-12
comparison and 1e-30 floor. This is restart-comparison evidence; source-level
specular reflection/escape conservation is recorded separately in
`spray-mpi-eb-specular-outflow-20260918.md`. No all-component application residual
qualification, external-solver parity or physical validation was run.

Artifacts: `build/spray-mpi-eb-specular-app-checks/{provenance.log,build.log,ctest.log,build-selected.log}`;
fixture `build/spray-amr-reactive-debug/tests/spray_mpi_eb_specular_outflow_application/impact_buwaeh7w/`.
