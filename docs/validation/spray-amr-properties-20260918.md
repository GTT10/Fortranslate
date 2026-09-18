# MPI AMR spray component properties — 2026-09-18

Implementation `e5203434bff68ca282417e0b841d10d454cafd95`, based on `015085a`.

Fixed/selected MPI AMR spray applications accept `&amr_spray_properties` for
component density/heat-capacity polynomial fits, Antoine pressure fits and
the existing Clausius surface cap. The shared liquid validator checks model
controls, and unused component entries must retain defaults. These initialized
liquids feed existing injection, particle/film thermodynamics and coupling.
The supplied group serializes into checkpoint context; omission preserves the
previous context. No constitutive equation or kernel tolerance changed.

## Completed verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing configurations unchanged.
Fixed `pelef_mpi_amr_spray_3d` and selected H2O2
`pelef_mpi_amr_spray_3d_selected` built. Selected executable was not run.
CTest **2/2 passed**, zero failures/skips, 9.85 seconds:
`spray_mpi_amr_properties_application` (4.83 s),
`spray_mpi_amr_restart` (4.98 s).

The mixed water/peroxide fixture configures polynomial density and heat capacity
for both components and Antoine pressure for water, over 300–450 K. Continuous
2-rank gas/parcel output matches checkpoint continuation on 4 ranks through
scheduled coarsening/injection at the existing `2e-12` relative and `1e-25`
absolute tolerances. A changed density coefficient rejects restart; a negative
density fit rejects the run before CSV output. Existing context, truncated
payload and immutable-checkpoint checks pass. The default-input restart test
also passes. Surface-cap behavior, films and alternative physical states were
not separately rerun. This checks application wiring and continuation, not
experimental accuracy of these synthetic fits.

Mechanism provenance, Python syntax, project contract and whitespace checks
passed. No failed builds/tests, changed tolerances, full matrix, external solver
parity, scaling or physical validation.

## Commands and artifacts

Logs: `build/spray-amr-properties-checks/{provenance,build,build-selected,ctest}.log`.
Run artifacts: `build/spray-amr-reactive-debug/tests/spray_mpi_amr_properties_application/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_(properties_application|restart)$' --output-on-failure
python3 -m py_compile tests/application/check_mpi_amr_spray_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
