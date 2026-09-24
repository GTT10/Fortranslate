# Dynamic LES timestep estimate on active faces — 2026-09-25

Implementation `b4e3e65`, based on `05fa47e`.
This resolves the bounded abrupt-Mach-start failure recorded in
[mach-inlet-mpi-application-20260925.md](mach-inlet-mpi-application-20260925.md).

## Cause and change

Temporary instrumentation located the rejection in the MPI Dynamic diffusion
gate. On the first half step (5e-9 s), one rank's padded maximum diffusivity was
47846.14066210483 m2/s, while the maximum on faces contributing to that rank's
active-cell divergence was only 1.5805789901747857e-9 m2/s. The padded maximum
produced a stability quantity of 717.69 against the existing 0.2 threshold.
Halving the step initially increased the padded estimate. The artificial periodic
workspace and unused outer donors were constraining cells that were not advanced.

The ghosted Dynamic wrapper now passes its interior bounds to the full-workspace
RHS routine through optional `diffusion_bounds(2,3)`. For direction d, the estimate
includes cell-face indices from lower(d)-1 through upper(d), and only interior
indices in the two transverse directions. Thus both physical boundary faces and
all interior divergence faces contribute; unused padding and the periodic seam
do not. The same mask applies to the doubled isothermal-wall diffusion estimate.
The full workspace still computes exactly the same coefficients, fluxes and RHS;
only the maximum used for timestep acceptance is restricted. Periodic callers
without bounds retain their preceding estimate. Existing thresholds, EOS checks,
retry limits and tolerances are unchanged. All padded states must still be
thermodynamically admissible.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug.
Fixed application and serial/MPI boundary tests built; selected H2O2 MPI spray
application built (selected execution not run). **4/4 CTests passed**, no failures
or skips, total 45.97 s:

- `spray_mpi_dynamic_mach_restart`: 23.51 s. Restored the original abrupt startup:
  stationary 800 K gas and Mach 0.2/900 K inlet, interior pressure. Completes two
  steps to 2e-8 s. Same-rank continuation is exact; 2-to-4-rank continuation meets
  existing 2e-11 relative/1e-25 absolute bounds. Relative mass/energy errors are
  -1.639549784590493e-17 and -1.6950467145357245e-16. Net inlet mass/energy are
  2.0361653148288622e-11 kg and 2.14804907981321e-5 J. Context/input guards pass.
- `spray_les_mpi_boundary_2`: 8.39 s; `_4`: 8.29 s. Exact serial donor/RHS/state
  comparisons, physical budgets and rollback retain their existing checks.
- `spray_les_dynamic_boundaries`: 5.74 s. Periodic behavior, serial/AMR-root
  comparisons, walls, four-model Mach checks and rollback pass.

Two diagnostic executions intentionally stopped at `max_steps=1` after reproducing
the reduced first step; their exit status is not counted as a passed test.
Instrumentation was removed and production targets rebuilt before the final
checks. Mechanism provenance, Python syntax of the application checker and
contract checker, project contract and whitespace passed. This establishes the
specific startup/timestep correction, not general turbulence or nozzle validation.

## Commands and artifacts

Logs, original diagnostic sources and instrumentation output:
`build/mach-startup-diagnosis/`. Final application artifacts:
`build/spray-amr-reactive-debug/tests/spray_mpi_dynamic_mach_restart/run-7mry6zvo/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_spray_3d test_les_dynamic_boundaries test_les_mpi_boundary --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_mpi_dynamic_mach_restart|spray_les_dynamic_boundaries|spray_les_mpi_boundary_(2|4))$' --output-on-failure
python3 -m py_compile tests/spray/test_mpi_spray_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

MPI environment overrides are host-local workarounds, not application defaults.
