# MPI AMR wall-film application integration — 2026-09-18

Implementation: `61ccaee11343b8fa2a1b8bff07629ad884af7ac6`, based on `df63496`.

The shared fixed/selected driver exposes `walls` in `amr_spray`. Any enabled
face activates film arrays and routes them through the coupled step, prescribed
or automatic remap, gas-work rebalance and checkpoint APIs. Defaults keep the
film-free path. Wall options and each port's impact viscosity are checked;
surface tension is required. Film settings bind restart context. Per-rank film
CSVs contain all film fields/components and are analysis output, not restart.
Output/checkpoint alias guards include film CSV paths.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja, `/usr/bin/gfortran`,
`/usr/bin/mpifort`, `/usr/bin/mpiexec`; unchanged fixed tests-enabled and selected
H2/O2 tests-disabled configurations reused.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_(film_application|restart)$' --output-on-failure --output-log build/spray-amr-film-app-checks/ctest-final.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_amr_spray_restart.py --exe build/spray-amr-reactive-debug/pelef_mpi_amr_spray_3d --mpiexec /usr/bin/mpiexec --source . --work build/spray-amr-film-app-checks/adaptive --multiport --physical --films --adaptive --rebalance
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_mpi_amr_spray_restart.py --exe build/spray-amr-reactive-selected-debug/pelef_mpi_amr_spray_3d_selected --mpiexec /usr/bin/mpiexec --source . --work build/spray-amr-film-app-checks/selected --selected --multiport --physical --films
python3 -m py_compile tests/application/check_mpi_amr_spray_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **2/2 passed**: film application **4.58 s**, existing film-free restart
**4.15 s**, total **8.77 s**. Fixed adaptive/rebalance film continuation and
selected film continuation each passed separately as manual checker invocations.
No build/test failures or skips. Provenance, syntax, contract and whitespace
checks passed. Existing compiler warnings remain. An initial edit helper stopped
on an assertion before writing files; replacing its multiline-call parser resolved
that preparation error.

Two ports and tied scheduled events deposit at low x and exit high x. The low
port starts 1 nm inside and moves toward the wall at 1 m/s; synthetic surface
tension/viscosity are 0.072 N/m and 0.001 Pa s. Four films exist at the 10 ns
checkpoint, eight at the 20 ns endpoint. Tests compare every gas/particle/film
CSV value and cumulative boundary ledger against uninterrupted two-rank results
after two-rank save/four-rank continuation, at relative 2e-12/absolute 1e-25.
They assert exact wall location/normal, positive film geometry/mass, and continued
loss of the saved film population's mass after restart. Changed face or port
properties, truncated payload and immutable overwrite are rejected.

The base film case performs prescribed coarsening. The adaptive variant instead
forms nine nodes; logged balancing actually migrates at steps 2 and 3 in both
full and resumed runs. The selected run uses the configured H2/O2 bundle.
Gas chemistry and molecular transport are disabled in these application fixtures;
film/particle transfer is active. Library coupled-transport evidence is separate.

Artifacts: `build/spray-amr-film-app-checks/` contains build logs, `ctest-final.log`,
`adaptive.log`, `selected-final.log` and earlier successful checks.
Adaptive cases: `adaptive/restart_37du73z5/`; selected cases:
`selected/restart_y8axq9qp/`, under that directory. Fixed CTest cases are below
`build/spray-amr-reactive-debug/tests/spray_mpi_amr_film_application/`.

Newly serialized wall settings may reject older application checkpoint contexts;
no blanket cross-build restart compatibility is claimed. These are bounded
implementation cases with synthetic liquids, not calibrated wall-impact/film
validation. EB particles, lateral film transport, convergence, external parity,
production scaling and broader completion qualification remain open.
