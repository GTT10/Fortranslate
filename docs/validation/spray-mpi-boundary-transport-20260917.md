# MPI molecular transport with physical boundaries — 2026-09-17

Source: parent `fbc0e06` plus this change. Logs and source hashes:
`build/spray-mpi-boundary-transport-checks/`. Existing GNU Fortran 13.3.0
Debug / Open MPI 4.1.6 configuration.

The serial physical transport routine now delegates flux construction and wall
energy correction to `physical_transport_ghosted_fluxes_3d`. The new helper accepts
two supplied ghost layers; MPI fills them from neighboring slabs and physical
boundary data. This includes neighbor species donor budgets in the shared serial
limiter calculation. The minimum-theta diagnostic includes padded cells, as in
the serial physical-boundary kernel.

`advance_mpi_spray_transport` accepts optional faces, inlet state/temperature,
wall temperature/velocity and a boundary ledger. It uses the shared flux helper
at both SSPRK2 stages. Only global x faces and local transverse faces enter the
surface budget. The returned change is positive into gas and replicated across
ranks. Ghost filling checks geometry and wall velocity collectively; invalid wall
thermal data rejects collectively without committing gas, temperature or budget.
The existing periodic path is retained when face labels are omitted. Supplied
physical boundary data remain caller inputs; model consistency across ranks is
required. Physical-boundary slabs require at least two cells in every dimension.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi test_reactive_boundary_transport_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^(spray_mpi_migration_[24]|unit_reactive_boundary_transport_3d)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Final tests: 3/3 passed, no failures or skips. The initial new test declaration
exceeded Fortran's line limit and was split. Two early multicomponent fixtures
failed initial thermodynamic conversion; the final H2/O2/N2 fixture sets pressure
from mixture gas constant at 800 K. Those failures are retained in separate logs;
no comparison tolerance was relaxed.

Uneven 2/4-rank partitions cover periodic, wall/outflow and inlet/wall cases,
nonuniform composition and velocity, 600 K isothermal walls and tangential wall
motion. Results match serial conserved fields within 2e-12 scaled error,
temperature within 1e-8 K and surface budgets within 1e-10 relative error with a
1e-12 floor. Independent volume-inventory changes agree with boundary fluxes
within 2e-12 of initial absolute inventories. Invalid thermal data on one rank
preserves gas and temperature exactly and returns zero boundary change.
The focused serial boundary-transport regression also passes. Provenance,
contract and whitespace checks pass; no Python source changed.

Nonperiodic particle/SGS stepping and application integration remain open.
This is numerical implementation evidence, not physical spray/LES validation
or a full toolchain configuration matrix.
