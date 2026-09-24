# Distributed regular LES physical halos — 2026-09-25

Implementation `450864f`, based on `7d262bf`.
`fill_mpi_les_boundary_ghosts_3d` builds 1..128-layer physical halos on uneven
regular x slabs. Interior state remains local; only the two halo sets of raw
x planes are requested from their owners. Shared coordinate mapping then applies
periodicity, reservoir inlet precedence, outflow clamping and repeated wall
reflections in the same axis/order convention as serial. Moving-wall momentum
changes retain their kinetic-energy correction. Temperature is copied from the
donor before optional pressure-target correction.

The serial mapping was extracted into `map_reactive_boundary_donor_3d`, used
by both serial and MPI fillers. The pressure helper now accepts a one-cell
local dimension: its physical donor indices refer to the global boundary in
the padded buffer, not a local slab edge. Corrections run on every slab whose
halo crosses the physical boundary. Global dimensions must still be at least
two in this MPI API. Empty spatial slabs remain unsupported.

Grid, width, faces, wall velocities, pressure targets and active inlet data
agree collectively. Local state shape/finiteness, temperature and allocations
are checked. Ghost state and temperature publish together only after all ranks
complete filling and pressure correction. This is the donor dependency for
regular MPI Dynamic physical fluxes; those fluxes and application dispatch
remain to be connected. No global field gather is introduced.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug.
Fixed unit targets and selected H2O2 spray runtime built. Selected execution
was not run. **3/3 tests passed**, total 1.58 seconds:

- `unit_amr_patch_tree_ghost_3d`: 0.25 seconds; retains default-width and
  analytic wide-wall coverage after extracting shared mapping.
- `spray_les_mpi_boundary_2`: 0.65 seconds.
- `spray_les_mpi_boundary_4`: 0.63 seconds.

MPI tests use a 5x3x2 global domain and seven halo layers, including one-cell
local x slabs. Periodic, reservoir-inflow/outflow, slip, adiabatic no-slip,
isothermal no-slip and pressure-outlet/moving-wall combinations match the
serial global-halo slice exactly for state and temperature. Opposing wall
velocities differ, exercising multiple reflections and corner composition.
Rank-disagreeing width, wall velocity and pressure reject without changing
previous ghost values. No failures, skips or tolerance changes occurred.
Mechanism provenance, Python syntax of the contract checker, project contract
and whitespace passed. This is distributed boundary-donor evidence, not
Dynamic flux/application execution, external solver parity or physical validation.

## Commands and artifacts

Logs: `build/mpi-les-boundary-checks/`. Launcher overrides are host-local.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_boundary test_amr_patch_tree_ghost_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(spray_les_mpi_boundary_(2|4)|unit_amr_patch_tree_ghost_3d)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
