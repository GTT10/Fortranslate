# MPI EB spray SDC particle acceleration — 2026-09-27

Implementation `769bd1d3`, based on `a108a792` on local `main`.
PeleC pin `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`;
PelePhysics pin `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.
H2O2 source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation and boundary

`advance_mpi_eb_spray_sdc` accepts optional `particle_acceleration(3)`.
It defaults to zero, is independent of the gas `force` argument, and requires
finite components bounded by 1e6 m/s² with exact agreement across ranks.
For common gravity, callers supply the same acceleration to both arguments.

A particle half kick precedes the old physical kick/drift; a second half kick
follows each provisional new physical kick. New evaluations always start from
the saved midpoint. Only the final accepted particle state and work publish.
Each external kick retains mass and internal energy and records represented
impulse `m*dv` and kinetic work `m*(v+dv/2).dv`. Evaporation and exits between
kicks naturally change the mass receiving the second kick. The global particle
work is added to the gas external source inventory; it is never deposited as
an internal gas/particle exchange. All fields and ledgers remain private until
completion, including final particle splitting and finite-ledger checks.

Pinned PelePhysics `Source/Spray/Drag.H:427` adds `body_force` to particle
acceleration while gas feedback uses drag. This implementation preserves that
external/internal distinction but brackets the existing physical integrator
with external half kicks; it is not the reference's simultaneous drag/body-force
update. No finite-step trajectory parity or order claim is made.
Application/restart wiring, film integration and recursive particle scheduling
remain open.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed-mechanism tests-enabled
`build/spray-amr-reactive-debug`, OpenMPI 4.1.6 through `/usr/bin/mpifort`
and `/usr/bin/mpiexec`. Existing real-comparison and executable-stack warnings
remain. Builds and tests had no failures.

CTest **2/2 passed**, 16.01 s total:

- `spray_mpi_eb_sdc_2`: 5.92 s, zero external acceleration regression,
  combined relative inventory residual **1.1968665763966046e-15**.
- `spray_mpi_eb_sdc_force_2`: 10.04 s, gas and particle acceleration
  `[1000,-2000,3000]` m/s², chemistry, transport, mixed-liquid heat/evaporation/
  drag, EB reflection and computational splitting. Combined budget after
  external work has relative residual **2.1293427042534180e-15** (bound 3e-12).
  Existing invalid-hydro and late split-capacity rollback checks also run
  with nonzero acceleration.

The force test additionally disables gas forcing, drag, heat and evaporation,
places the parcel inside fluid, and checks one and three SDC iterations against
constant-acceleration motion: velocity error below 1e-14 m/s, position error
below 1e-18 m, external impulse/work relative error at most 1e-13.
Internal gas-particle exchange is exactly zero. Thus the work and motion checks
independently detect omitted or repeated particle acceleration. Rank-disagreeing
acceleration fails without publishing source outputs.

Provenance, Python syntax, project contract and diff whitespace passed.
No selected-mechanism rerun, full matrix, installed audit, upstream solver
parity, order study or physical spray validation was performed.

## Commands and artifacts

Artifacts: `build/spray-sdc-force-checks/` (`provenance.log`, `build.log`,
`build-final.log`, `ctest.log`, `detail.log`, `contract.log`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc(_force)?_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
