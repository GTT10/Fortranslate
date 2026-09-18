# MPI EB impact/film application integration — 2026-09-18

Implementation `7b3412d`, base `9e5f8ae`.

Fixed and selected MPI EB spray executables expose wall impacts and wall options,
connect the reactive impact/film composition, include films in coupled budgets,
write protected per-rank film CSVs, and pass films through remap/rebalance and
checkpoint restore. The normalized namelist binds wall settings to the generation
context. Impact restart admits retired parents and four-child ID reservations.
The application exits immediately on an empirical-model rejection.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpifort`, `/usr/bin/mpiexec`),
Debug/Ninja. Existing fixed tests-enabled and selected H2O2 tests-disabled
configurations were reused.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_eb_spray_3d --parallel 2
python3 tests/application/check_mpi_eb_impact.py --exe build/spray-amr-reactive-debug/pelef_mpi_eb_spray_3d --mpiexec /usr/bin/mpiexec --source . --work build/spray-mpi-eb-impact-app-checks/fixed
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_eb_impact.py --exe build/spray-amr-reactive-selected-debug/pelef_mpi_reactive_eb_tree_3d_selected_spray --mpiexec /usr/bin/mpiexec --source . --work build/spray-mpi-eb-impact-app-checks/selected --selected
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_eb_spray.py --exe build/spray-amr-reactive-debug/pelef_mpi_eb_spray_3d --mpiexec /usr/bin/mpiexec --source . --work build/spray-mpi-eb-impact-app-checks/legacy
python3 -m py_compile tests/application/check_mpi_eb_impact.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

All **three application drivers passed**, no unexpected failures or skips.
Each impact driver performs eight MPI launches: full four-rank, stopped two-rank,
resumed four-rank and deliberate wall-context mismatch, for both deposition and
splash. The legacy driver performs three launches for specular rank parity and
missing-checkpoint rejection. Both executable builds, provenance, Python syntax,
contract and whitespace checks passed. The fixed driver is registered as
`spray_mpi_eb_impact_application`; this run invoked it directly, not through CTest.

The fixture uses two initial water parcels, a planar nodal wall, 400 K wall
heat transfer, implicit chemistry, all molecular transport and HLLC hydro.
Deposition at 20 m/s leaves two films and no particles. Splash at 200 m/s with
1e-16 kg total injection leaves two films and eight outgoing particles (IDs 3–10).
Both film temperatures exceed 350 K; component mass equals film mass. Prescribed
regrids move the fine patch, coarsen to the root and restore refinement. A
checkpoint after step two contains V2 particles/films; continuation at four
ranks matches gas, particle, film and cumulative-budget columns within relative
3e-12, scale floor 1e-30. A changed wall temperature rejects restart before output.

Rebalancing was enabled every step, but these application logs show no accepted
owner change. Actual film ownership redistribution has separate library evidence
in `spray-mpi-eb-film-repartition-20260918.md`; this application run qualifies
routing and the no-change path, not an observed load-balancing move.

Artifacts under `build/spray-mpi-eb-impact-app-checks/`: `build.log`,
`build-selected.log`, `fixed.log`, `selected.log`, `legacy.log`, `provenance.log`,
`contract.log`; inputs, outputs, manifests, payloads and run logs are under
`fixed/impact_dllxqubf`, `selected/impact_fnfp15jn`, `legacy/eb_spray_knxvco84`.

The test retains the empirical thin-film boundary; it does not qualify arbitrary
wet-wall impacts, film lateral flow, adaptive tagging, multiport/profile impacts,
large-scale runs or experimental physics. No broad matrix, installed audit or
external solver parity was run.
