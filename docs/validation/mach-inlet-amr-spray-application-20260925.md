# MPI AMR spray application Mach inlet integration — 2026-09-25

Implementation `d89c8d0`, based on `ae4ae60`.

Fixed/selected inputs enable Mach controls. Active faces bypass reservoir
construction; per-step distributed ghosts provide hydro/molecular bounds with
the application's existing finest-level policy. The LES selector and both
wall-film/non-film coupled calls receive Mach and masked inlet pressure. The
existing full-configuration checkpoint context binds these controls without a
format change.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug. Fixed and
selected H2O2 applications built; selected execution was not run.
`spray_mpi_amr_dynamic_mach` passed: **1/1 tests**, no failures or skips,
17.80 s (CTest total 17.84 s).

The fixture uses Dynamic LES, Mach 0.2/1000 K with extrapolated pressure,
three injection events, scheduled coarsening and 2-to-4-rank checkpoint
continuation. It retains the existing state/particle/budget continuation checks,
context and truncated-checkpoint/overwrite guards. Changed Mach, inlet pressure
and temperature reject restart. Chemistry and molecular transport are disabled
in this application fixture; their coupled-operator checks are separate. Wall
films and selected execution were not exercised. Provenance, Python syntax,
project contract and whitespace passed. No tolerance changed. These bounded
checks do not establish nozzle, diesel ignition or general turbulence validation.

## Commands and artifacts

Logs: `build/mach-inlet-amr-spray-app-checks/`, containing `build-fixed.log`,
`build-selected.log` and copied `ctest.log`. Application artifacts:
`build/spray-amr-reactive-debug/tests/spray_mpi_amr_dynamic_mach/restart_*/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_dynamic_mach$' --output-on-failure
python3 -m py_compile tests/application/check_mpi_amr_spray_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Environment overrides are host-local workarounds, not application defaults.
