# MPI AMR gas application Mach inlet integration — 2026-09-25

Implementation `0e3702d`, based on `0ddcbd6`.

## Change

Fixed and selected MPI tree frontends enable Mach input. The application skips
reservoir construction for active Mach faces and forwards temperature, Mach and
masked inlet pressure to coupled stepping. Each step packs only owned nodes,
exchanges the required ghosts and computes local hydro/molecular timestep bounds.
An MPI minimum supplies the shared bound; empty ranks participate without owning
payloads. Recursive mode uses refinement and squared-refinement scaling matching
the serial application. The existing full-configuration checkpoint context binds
the inlet controls. MPI AMR spray application integration remains pending.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug. Fixed and
selected elementary serial/MPI targets built and executed.
`regression_mpi_amr_tree_mach_restart_3d` passed: **1/1 tests**, no failures or
skips, 59.33 s (CTest total 59.37 s).

Mach 0.2/900 K with extrapolated pressure, chemistry and molecular transport
completes three common steps and three scheduled regrids to 1e-9 s. Fixed and
selected variants preserve exact state/budget CSV bytes across continuous serial
and MPI runs, 2-to-3-rank continuation, MPI-to-serial continuation and serial-to-
one-rank continuation. Serial/MPI checkpoint bytes match. A root-only restart
with two empty ranks passes. Changed Mach, inlet pressure and temperature reject
restart, alongside existing schedule/physics/checkpoint guards. Provenance,
Python syntax, project contract and whitespace passed. No tolerance changed.
This fixture checks common-step integration; recursive operator evidence is
separate. Its short duration does not stress the timestep cap or validate nozzle
physics.

## Commands and artifacts

Build log: `build/mach-inlet-mpi-amr-app-checks/build.log`.
Application logs and checkpoints:
`build/spray-amr-reactive-debug/tests/mpi_amr_tree_mach_restart_3d/{fixed,selected}/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_tree_reactive_3d test_selected_elementary_mpi_amr_tree_3d pelef_amr_tree_reactive_3d test_selected_elementary_amr_tree_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^regression_mpi_amr_tree_mach_restart_3d$' --output-on-failure
python3 -m py_compile tools/check_mpi_amr_tree_restart_3d.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Environment overrides are host-local workarounds, not application defaults.
