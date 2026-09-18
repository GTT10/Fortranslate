# MPI AMR spray rate-derived speed — 2026-09-18

Implementation `035c7d2d463fb0efad08d726951a6cd20185cec9`, based on `b2e1c2d`.

Fixed/selected MPI AMR spray applications accept per-port `rate_velocity` and
`discharge_coefficients` through `&amr_spray_velocity`. Enabled ports require
rate profiles and finite bounded radii/coefficients. Each event's speed uses
its integrated mass divided by interval length, mixture density at the port
temperature, circular area and discharge coefficient. Zero-mass events use
zero speed; positive events reject nonpositive denominators or speed above
`1e7 m/s`. The existing mixture property API supplies density, including fits.
The supplied group binds checkpoint identity. Defaults retain prescribed speed
and the previous checkpoint context. No pressure-drop model is introduced.

## Completed verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing configurations unchanged.
Fixed `pelef_mpi_amr_spray_3d` and selected H2O2
`pelef_mpi_amr_spray_3d_selected` built; selected executable was not run.
CTest **2/2 passed**, zero failures/skips, 9.73 seconds:
`spray_mpi_amr_discharge_application` (5.68 s), `spray_mpi_amr_restart` (4.01 s).

Two opposite-axis ports use different component fractions, radii, rate profiles
and discharge coefficients. With zero cone angle/swirl, the accumulated axial
injection momentum matches an independent calculation of the four event masses
times their rate-derived speeds within `2e-12` relative. The rate fixture also
checks total mass `8.125e-11 kg`. Continuous 2-rank gas/parcel output matches
2-to-4 rank continuation through scheduled coarsening/injection within existing
`2e-12` relative / `1e-25` absolute tolerances. Changed discharge coefficient
rejects restart; zero coefficient rejects input. Existing rate/context,
truncated-payload and checkpoint overwrite checks pass. Default restart passes
separately. Temperature-dependent density fits were not separately rerun in
this fixture; this is not physical nozzle qualification.

Provenance, Python syntax, project contract and whitespace checks passed.
No failed builds/tests, changed tolerances, full matrix, external solver parity,
scaling or physical validation.

## Commands and artifacts

Logs: `build/spray-amr-discharge-checks/{provenance,build,build-selected,ctest}.log`.
Run artifacts: `build/spray-amr-reactive-debug/tests/spray_mpi_amr_discharge_application/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_(discharge_application|restart)$' --output-on-failure
python3 -m py_compile tests/application/check_mpi_amr_spray_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
