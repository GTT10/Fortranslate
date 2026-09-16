# Sparse MPI AMR physical molecular transport and full split — 2026-09-17

Source: parent `2a7bb0d` plus this change. Logs and source hashes:
`build/amr-mpi-physical-transport-checks/`. Existing GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration; launcher `/usr/bin/mpiexec.openmpi`.

Sparse MPI transport now accepts physical faces, reservoir data, wall temperature
and wall velocity. The coarse Euler stages compute physical fluxes from the
rank-local padded state; internal rank faces remain neighbor interfaces.
Domain surface transfer includes only actual x-end faces and the slabs' y/z
faces, then is reduced once across ranks. Fine transport uses physically padded
coarse endpoints with the shared physical interpolation path. Existing interface
limiting, reflux and average-down retain sparse ownership. SSPRK2 returns the
mean stage transfer only after committing both levels.

The full reaction / transport / hydro / transport / reaction split forwards
boundary data and sums the two transport transfers plus hydro transfer privately.
Failure leaves the hierarchy unchanged and resets boundary/reflux diagnostics;
transport-disabled calls forward to the physical reaction/hydro split. Boundary
option presence is part of the collective algorithm contract. Wall-temperature
presence and values also agree collectively before local flux evaluation.

The shared physical ghosted transport kernel now permits one-cell dimensions
and visits lower and upper faces independently. Previously its cell loop selected
only the upper face when a dimension had one cell. Explicit face indices now
apply wall conduction and traction work at both ends; ordinary multi-cell
arithmetic is unchanged. The serial boundary wrapper still supplies its usual
valid ghost layers. This enables one-column MPI slabs without a full-grid gather.

The MPI near-wall fixture adds transport-only and active full-split cases with
viscosity, conduction, species diffusion and barodiffusion. Opposing x walls are
600/900 K with tangential velocities +20/-10 m/s, alongside transverse inflow and
no-slip/outflow faces. A ratio-two patch covers coarse cells 2:3 in an 8-cubed
coarse grid. Coarse/fine states, temperatures and reflux are compared exactly to
serial APIs. Boundary transfers, composite Euler/species or elemental balances,
chemical activity and rank-local invalid-inlet rollback are checked. Four/eight
ranks exercise ranks without fine cells and one-column coarse slabs. Existing
periodic tests and the serial physical transport test are also selected.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_mpi_amr_reactive_3d test_amr_reactive_transport_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^(mpi_amr_reactive_3d_np(4|8)|unit_amr_reactive_transport_3d)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 3/3 tests passed with no skips (MPI four ranks 62.60 s, eight ranks
53.68 s, serial transport 17.02 s; total 133.32 s). Provenance, project contract
and diff checks passed.

No existing tolerance or Python source changed. The full hierarchy APIs are
connected, but MPI application inputs, boundary-aware timestep selection,
diagnostics and checkpoint persistence still need activation. This remains a
strictly interior two-level numerical check, not physical or scaling validation.
