# MPI EB spray outflow application/restart — 2026-09-18

Implementation `ef1d947`, base `880e08d`.

The MPI fixed/selected EB spray application accepts `particle_outflow` in
`eb_spray`, false by default. It requires impact films and forwards the flag
to reactive stepping. Namelist serialization binds it into restart context.
Existing gas and particle/film checkpoint formats retain exterior history,
sparse surviving IDs and the watermark; no format change was needed.

## Completed focused checks

GNU Fortran 13.3.0/Open MPI 4.1.6, Debug/Ninja, reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. Selected executable compiled only.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_outflow_application$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
python3 -m py_compile tests/application/check_mpi_eb_impact.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **1/1 passed**, 0 failures/skips, 20.23 s total. Four launches cover
4-rank full, 2-rank checkpoint, 4-rank resume and changed-context rejection.
Three scheduled regrids and rebalancing remain enabled. Gas, particle, film and
budget records match full versus resumed runs using the existing relative
3e-12 comparison with 1e-30 floor. The full run retains two films and four
particles (IDs 3,4,7,8) after half the splash children escape. The changed wall
temperature case rejects before output. Builds and static checks passed.

A separate ad-hoc assertion that every budget residual be below
3e-12*max(abs(initial),abs(final),1e-30) **failed** for near-zero momentum.
This was not a CTest check and was not relaxed or counted as passed. For example,
x momentum final is 6.0047277677e-16 with residual -2.0471168274e-24;
the y momentum balance contains cancelling wall/boundary terms around 4.32e-8.
Mass boundary change is -2.7426979062e-17 kg and mass residual 1.7520689723e-23 kg.
Thus this record establishes application wiring and restart comparison, not a
new all-component application conservation qualification. Focused library
conservation evidence remains in `spray-mpi-eb-outflow-reactive-20260918.md`.

Artifacts: `build/spray-mpi-eb-outflow-app-checks/{provenance.log,build.log,ctest.log,build-selected.log}`;
fixture `build/spray-amr-reactive-debug/tests/spray_mpi_eb_outflow_application/impact_2dxnjuwh/`.
Specular MPI EB exterior modes and broader physical qualification remain open.
