# MPI AMR spray rate profiles — 2026-09-18

Implementation `2eb747877c095763d9938387675e21a20b0e1de9`, based on `f489682`.

Fixed/selected MPI AMR spray applications read optional per-port rate tables and
event end times through `&amr_spray_rates`. The existing piecewise-linear
integrator computes event masses before schedule validation. Profile-driven
events require zero explicit mass, avoiding ambiguous input precedence. Knot
counts, unused entries, intervals and final mass/count consistency are checked.
The supplied group and resolved event masses bind restart identity. Omission
retains the previous behavior. The whole interval mass is injected at its start;
this does not introduce continuous source deposition, automatic event refinement
or rate-driven velocity.

## Completed verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing configurations unchanged.
Fixed `pelef_mpi_amr_spray_3d` and selected H2O2
`pelef_mpi_amr_spray_3d_selected` built; selected executable was not run.
CTest **2/2 passed**, 9.91 seconds, zero failures/skips:
`spray_mpi_amr_rates_application` (5.87 s), `spray_mpi_amr_restart` (4.00 s).

Two independent port profiles feed four events, including two simultaneous
events on different ports. Total injected mass matches the analytic trapezoidal
sum `8.125e-11 kg` within `2e-12` relative. Continuous 2-rank gas/parcel outputs
match 2-to-4 rank continuation through scheduled coarsening and injection within
the existing `2e-12` relative / `1e-25` absolute tolerances. A changed rate rejects
restart and a negative rate rejects the run before CSV publication. Existing
context/truncation/overwrite checks and the default-input restart case pass.
No temporal-convergence or physical injector validation is claimed.

Provenance, Python syntax, project contract and whitespace checks passed.
No failed builds/tests, changed tolerances, full matrix, external solver parity
or scaling qualification.

## Commands and artifacts

Logs: `build/spray-amr-rates-checks/{provenance,build,build-selected,ctest}.log`.
Run artifacts: `build/spray-amr-reactive-debug/tests/spray_mpi_amr_rates_application/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_(rates_application|restart)$' --output-on-failure
python3 -m py_compile tests/application/check_mpi_amr_spray_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
