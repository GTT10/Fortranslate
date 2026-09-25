# AMR tree hydro artificial-viscosity controls — 2026-09-25

## Source and scope

Code commit `2b4b327`, parent `19dd9f6`. GNU Fortran 13.3.0, OpenMPI,
Ninja Debug, Python 3.12.3; existing fixed-mechanism tests-enabled build
`build/spray-amr-reactive-debug`.

Serial and MPI patch-tree common-step and single-level hydro now accept a
trailing optional `difmag`. Omission resolves to 0.1 for PPM/pelec and zero
otherwise. Explicit values must be finite and nonnegative; positive values
require PPM. The resolved coefficient reaches final halo CTU flux correction.
MPI checks it as part of collective control agreement. Serial/MPI recursive
subcycle entry points forward the optional coefficient to every level call.
Failures preserve the public state, temperatures, ledgers and step counters.

This is the hydro API connection. Reaction/transport composition, application
namelist enabling, AMR spray forwarding and restart binding are still pending.
The older two-level AMR stepping family was not extended in this increment.

## Completed checks

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_hydro_3d test_mpi_amr_patch_tree_subcycle_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(unit_mpi_amr_patch_tree_ppm_3d_2|unit_mpi_amr_patch_tree_ppm_subcycle_3d_2)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**: common-step 4.27 s; recursive 13.27 s; total elapsed 17.59 s.
Build, provenance, Python syntax, contract and whitespace checks passed.
No failed checks in this increment. Compiler warnings remain.

Existing PPM fixtures now use coefficient 0.037 for both serial and MPI
references. They cover periodic and physical boundaries, prescribed inflow,
Mach inflow with extrapolated/prescribed pressure, pressure outlets, alternate
owner maps and ranks without owned patches. Conserved state, temperature,
boundary budgets and reflux diagnostics remain exactly equal to serial.
Recursive updates retain exact level counts 1/2/6. A rank-dependent coefficient
and a negative value on one rank both reject collectively; the existing
unchanged-state helper verifies gas, temperature and diagnostics, plus recursive
counters. Existing bad solver, timestep and boundary-control rejection checks
also pass. No tolerance was changed.

Artifacts: `build/difmag-amr-tree-hydro-checks/build.log`, `test.log`, and
`LastTest.log` (copied before any later CTest run).

## Evidence limits

These are bounded fixed-mechanism serial/MPI hydro checks. No selected build,
application/restart run, full configuration matrix, installed audit, external
PeleC solver parity or physical validation was performed. The coefficient kernel
and default dispatch have separate earlier evidence; this record does not count
those checks as rerun.
