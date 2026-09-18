# LES model selection in coupled spray updates — 2026-09-18

Implementation `4d5d125`, base `6748a0a`.

`advance_spray_reactive_3d` and `advance_mpi_mixture_reaction` accept optional
trailing `les_model` (`smagorinsky`, `wale`, `vreman`; omitted means Smagorinsky).
Both LES half steps receive the selection. Serial validates the name before
particle/gas work; MPI includes its normalized code in collective coupled-control
agreement, even when the coefficient disables LES. Existing candidates retain
whole-step publication of gas, parcels, mixture/film data and diagnostics.
Application input, timestep-driver dispatch and checkpoint identity are pending.

## Focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6; reused Debug/Ninja fixed tests-enabled and
selected H2O2 tests-disabled configurations. No MPI runtime test in this increment.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_integration --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_integration_coupled(_wale|_vreman)?$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The existing 5^3 reacting/evaporating coupled fixture retains its Smagorinsky
case and adds WALE/Vreman with crossed sinusoidal gas velocities. Tests retain
mass/energy relative 1e-12 and momentum absolute 1e-18 conservation gates,
require active evaporation and chemistry, and verify invalid model rejection
preserves gas/temperature and parcel masses/positions with zero LES coefficient.
MPI coupled forwarding is compiled in fixed/selected runtimes; its runtime
qualification is separate from prior isolated MPI LES parity.

CTest **3/3 passed**, zero skips, 14.38 s total. Fixed/selected builds,
provenance, Python syntax, project contract and diff checks passed. No failed
or interrupted runs, no existing tolerance changes, no physical-validation claim.

Logs: `build/les-coupled-model-checks/{provenance,build,build-final,build-selected,ctest}.log`.
