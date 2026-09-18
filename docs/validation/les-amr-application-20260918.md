# MPI AMR spray LES application — 2026-09-18

Implementation `fbaa855fba45af2e6a0d38b7d44aa40bbbc6cb4e`, based on `34f8a1c`.

The shared fixed/selected MPI AMR spray application reads optional
`&amr_spray_les` controls for static model, eddy coefficient, CI, PrT, ScT and
diffusion CFL. All are validated. Nondefault controls serialize into checkpoint
identity; omitted/default controls retain the previous serialization. Both
film-enabled and ordinary branches forward LES controls to the transactional
coupled integrator. The distributed LES interval estimate limits each proposed
application step, after regridding/injection and before retry. Later stage
checks remain active, including after particle feedback. LES uses common half
steps even when the gas operator subcycles.

## Completed verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing configurations unchanged.
Fixed `pelef_mpi_amr_spray_3d` and selected H2O2
`pelef_mpi_amr_spray_3d_selected` built. Selected executable was not run.
CTest **2/2 passed**, zero failures/skips, 9.99 seconds:
`spray_mpi_amr_les_application` (4.93 s) and `spray_mpi_amr_restart` (5.03 s).

The LES case uses Vreman 0.17 and CI 0.1. Continuous 2-rank output matches
checkpoint continuation on 4 ranks within the existing `2e-12` relative and
`1e-25` absolute tolerances for gas and parcels, through scheduled coarsening
and injection. A changed CI rejects restart before CSV publication. Existing
spray-context, truncated payload and immutable checkpoint rejection checks
also pass, as does the default-disabled restart fixture. These short event-bound
steps exercise selector dispatch but do not prove operation at the LES stability
limit. Film-plus-LES, adaptive/rebalance-plus-LES, physical boundaries and the
other two closures were not separately run at application level in this change.

Provenance, Python syntax, project contract and whitespace checks passed.
No failed builds/tests occurred and no tolerances changed. No full matrix,
external solver parity, scaling or physical validation was performed.

## Commands and artifacts

Logs: `build/les-amr-application-checks/{provenance,build,build-selected,ctest}.log`.
Per-run inputs, logs and CSV/checkpoint artifacts are under
`build/spray-amr-reactive-debug/tests/spray_mpi_amr_{les_application,restart}/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_(les_application|restart)$' --output-on-failure
python3 -m py_compile tests/application/check_mpi_amr_spray_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
