# MPI spray Mach input and changed-rank restart — 2026-09-25

Implementation `5af0e83`, based on `472458c`.
Regular fixed/selected MPI spray namelists accept `inflow_mach(6)` (default -1),
`inflow_temperature(6)` (1200 K) and `inflow_pressure(6)` (101325 Pa). Active Mach
0 through 100 requires an inlet face, positive temperature and pressure from
zero to 1e15 Pa. Zero pressure extrapolates the interior. Mach faces replace the
initial-reservoir velocity/composition using current interior donors; inactive
inlets retain their original initial-gas reservoir. Controls reach the complete
coupled MPI step. Distributed current inlet ghosts contribute to hydro and
molecular diffusion bounds before the global timestep reduction. Active Mach,
temperature and pressure bind checkpoint context; default inactive inputs retain
the previous encoding.

## Completed checks and limitation

GNU Fortran 13.3.0, OpenMPI 4.1.6 via cached `/usr/bin/mpiexec`, Ninja Debug.
Fixed and selected H2O2 MPI spray applications built; selected execution was not
run. Final `spray_mpi_dynamic_mach_restart`: **1/1 passed**, 23.19 seconds, total
23.22, no failures/skips. The fixture has 8x4x4 cells, four mixed parcels,
molecular transport, implicit chemistry and Dynamic LES (Gaussian test ratio 4,
grid ratio 2). It advances two steps to 2e-8 s. Same-rank continuation matches
serialized gas/particle values exactly; 2-to-4-rank continuation meets existing
2e-11 relative/1e-25 absolute bounds. Relative mass/energy errors are
2.644972160187063e-16 and -6.322278091200028e-16. This is nearly balanced
throughflow: net boundary mass is only 9.69e-26 kg. Changed/removed Mach, changed
inlet temperature or pressure, and existing model/filter changes reject restart
before history creation. An active Mach on a non-inlet face rejects input.

**Unresolved startup limitation:** the initial fixture used zero initial gas
velocity, 800 K gas and a Mach 0.2/900 K inlet. It repeatedly reduced its timestep
(the first accepted step was 1.52587890625e-13 s) and exceeded the checker's
90-second timeout. No completed qualification is claimed for that startup.
The final fixture keeps Mach 0.2 but sets both gas/inlet to 800 K and initial
normal velocity to 132.615839434881 m/s, derived from the pinned NASA7 sound speed
for H2/O2/N2 mass fractions .03/.23/.74. This removes the initial inlet jump.
The abrupt-start behavior still requires investigation; no limiter, tolerance,
retry bound or timeout was relaxed. The timed-out MPI processes were confirmed
absent before retrying. Initial artifacts are retained alongside the final run.

Mechanism provenance, Python syntax of the application checker and contract
checker, project contract and whitespace passed. This is bounded MPI application
and restart evidence, not AMR integration, abrupt-start robustness, nozzle parity
or physical validation.

## Commands and artifacts

Logs: `build/mach-inlet-mpi-app-checks/`, including initial/final CTest output.
Final artifacts:
`build/spray-amr-reactive-debug/tests/spray_mpi_dynamic_mach_restart/run-_u3mlxfp/`.
Timed-out fixture: sibling `run-xmz3b2f3/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_dynamic_mach_restart$' --output-on-failure
python3 -m py_compile tests/spray/test_mpi_spray_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

MPI environment overrides are host-local workarounds, not application defaults.
