# Boundary-touching AMR molecular transport — 2026-09-17

Source: parent `1a55084` plus this change. Logs and source hashes:
`build/amr-boundary-transport-checks/`. Reused GNU Fortran 13.3.0 Debug build
with bounds and floating-point checks.

Serial AMR transport now accepts patches touching physical boundaries. Fine
subcycles construct two physical/coarse-fine ghost layers and evaluate the
existing physical transport kernel at fine spacing. Only actual domain faces
receive wall heat flux and tangential traction work; internal patch faces keep
their coarse-interpolated ghost data. Fine subcycling still uses ratio squared.
Interface positivity correction applies only where an uncovered coarse neighbor
exists. Physical surface fluxes replace covered coarse contributions through
the boundary-aware reflux API, followed by average-down and temperature recovery.
The enclosing SSPRK2 and full reaction/transport/hydro split retain private
candidates, cumulative boundary transfers and rollback.

The fine physical ghost builder is shared with hydro. Its viscous mode reflects
all momentum components at no-slip walls, composes wall velocity shifts at
corners, and adjusts total energy to preserve internal energy. Hydro retains
normal reflection. Reservoir and outflow handling keep the established ordering.
The state-based transport timestep helper now accepts valid boundary-touching
geometry; callers still supply appropriate wall/inlet timestep restrictions.

The transport fixture adds a boundary-face patch, a three-face corner at ratio
three, and full-domain refinement. It exercises cold moving walls, transverse
moving adiabatic walls, inflow, transport-only stepping and the full split with
chemistry disabled/enabled. Checks cover composite species/Euler or H/O/N
balance after subtracting boundary transfer, energy removal at a cold wall,
limiter bounds and exact rollback of states and diagnostics for invalid wall
temperature or inlet density. Full-domain viscous ghost arrays are compared
exactly with the existing regular-grid boundary fill. The hydro/chemistry unit
is rerun because its ghost construction now uses the shared builder; it retains
the exact full-refinement/regular-grid hydro comparison.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_amr_reactive_transport_3d test_amr_reactive_chemistry_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^unit_amr_reactive_(transport|chemistry)_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 2/2 tests passed, no skips (transport 76.74 s, hydro/chemistry 57.02 s;
total 133.79 s). Provenance, project contract and diff checks passed. No Python
source changed.

MPI, regridding, checkpoint and application activation for boundary-touching
patches remain open. Periodic seams are still rejected by stepping. This is
bounded numerical evidence, not physical or scaling validation.
