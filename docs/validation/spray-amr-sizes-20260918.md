# MPI AMR spray size inputs — 2026-09-18

Implementation `bbccf6eb04d607238fec7f235b7ff075ad160283`, based on `3c4b362`.

The shared fixed/selected MPI AMR spray application reads `&amr_spray_sizes`
with `port_sizes(16)`, passes all six existing distribution controls to the
MPI cone injector and binds a supplied group to checkpoint identity. Omission
preserves previous input/context behavior. The sampler, component-mass and
global-ID logic are reused unchanged; event validation remains in the injector.

## Completed verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing configurations unchanged.
Fixed `pelef_mpi_amr_spray_3d` and selected H2O2
`pelef_mpi_amr_spray_3d_selected` built; selected executable was not run.
CTest **2/2 passed**, 9.09 seconds, zero failures/skips:
`spray_mpi_amr_sizes_application` (4.98 s), `spray_mpi_amr_restart` (4.07 s).

The two-port fixture uses a truncated normal mass-weighted distribution on one
port and a truncated Weibull number-weighted distribution on the other. Each
port produces varying individual drop masses (maximum exceeds 1.1 times minimum).
Continuous 2-rank gas/parcel output matches 2-to-4 rank continuation through
scheduled coarsening and injection at the existing `2e-12` relative / `1e-25`
absolute tolerances. A changed standard deviation rejects restart; a negative
standard deviation rejects injection before CSV output. Existing context,
truncated-payload and immutable-checkpoint checks pass, as does the separate
default-input restart test. Lognormal and chi-squared were not independently
rerun at AMR application level; this is not distribution convergence or physical
injector qualification.

Provenance, Python syntax, project contract and whitespace checks passed.
No failed builds/tests, changed tolerances, full matrix, external solver parity,
scaling or physical validation.

## Commands and artifacts

Logs: `build/spray-amr-sizes-checks/{provenance,build,build-selected,ctest}.log`.
Run artifacts: `build/spray-amr-reactive-debug/tests/spray_mpi_amr_sizes_application/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_(sizes_application|restart)$' --output-on-failure
python3 -m py_compile tests/application/check_mpi_amr_spray_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
