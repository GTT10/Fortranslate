# Dynamic LES AMR pressure outlets — 2026-09-25

Implementation `3c1d0c8`, based on `dc71c61`. The shared acoustic pressure-ghost correction now accepts
halo widths 1 through 128 (default 2) and explicit physical-domain donor indices.
The AMR wrapper locates these indices from global integer coordinates. Thus
even an interior patch whose wide halo crosses a domain boundary receives the
correct physical pressure extrapolation. In-domain donor cells remain unchanged.
Face ordering at corners and the underlying acoustic model are preserved.

Serial/MPI Dynamic AMR selectors and SSPRK2 now apply pressure correction to
their private wide halos. Coupled spray and application dispatch accept the
same nonzero pressure targets on outflow faces. Existing admissibility checks,
collective failure, conservative reflux and rollback remain in force.

Large pressure differences can produce an inadmissible outer ghost because
the acoustic extrapolation scales with distance from the physical boundary.
This remains a rejected state, not a clamped or relaxed boundary model.

## Checks and scope

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug. Fixed tests-enabled targets and
selected H2O2 tests-disabled application built. Selected pressure-outlet
execution and a full configuration matrix were not performed.

Final focused CTest: **3/3 passed**, no failures or skips, total 123.93 seconds:
pressure-halo unit 0.06 seconds, pressure-outlet application 16.55 seconds and
MPI LES comparison 107.28 seconds. The latter verifies composite boundary
conservation, serial/MPI timestep and state agreement, empty/distributed owners
and rollback. Application checks include 2-to-4-rank restart and rejection of
changed filter context, truncated checkpoints and overwrite attempts.
Provenance preflight, Python syntax, project contract and whitespace passed.
The halo unit also passed directly before the final CTest run.

The new `unit_amr_pressure_halo_3d` places a four-cell child strictly inside an
eight-cell domain. Its seven-cell halo crosses both ends. Every padded cell
matches the acoustic state at its known physical distance, including unchanged
in-domain cells; invalid width rejects without changing the buffer.
The MPI LES test adds a 110000 Pa target to its nonuniform mixed-inflow/outflow
case, retaining periodic, extrapolated-outflow, moving-wall and static cases.
The application checker uses a 100000 Pa outlet with a uniform initial pressure
of 101325 Pa, scheduled regridding and changed-rank continuation.

The first CTest pair had one pass (application, 16.83 seconds) and one failure
(MPI LES, 7.92 seconds). In the original nonuniform test fixture the maximum
interior pressure was 109431 Pa, so seven-layer extrapolation toward 100000 Pa
would give `109431 + 14*(100000-109431) = -22603 Pa`. The algorithm rejected
the state. The fixture target was changed to 110000 Pa, preserving the solver's
admissibility checks and all comparison tolerances.

These are bounded implementation checks, not external PeleC parity, general
nonreflecting-boundary validation, physical qualification or selected-fuel
pressure-boundary evidence. Dynamic EB and LES subcycling remain open.

## Commands and artifacts

Build and initial/final test logs are under
`build/les-dynamic-amr-pressure-checks/`.
Final application artifacts:
`build/spray-amr-reactive-debug/tests/spray_mpi_amr_dynamic_pressure/restart_hy1x0a4b/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_amr_pressure_halo_3d test_les_mpi_amr_tree pelef_mpi_amr_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(unit_amr_pressure_halo_3d|spray_mpi_amr_dynamic_pressure|spray_les_mpi_amr_tree_2)$' --output-on-failure
python3 -m py_compile tests/application/check_mpi_amr_spray_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
