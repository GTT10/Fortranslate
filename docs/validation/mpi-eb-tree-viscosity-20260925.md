# MPI common-step EB tree viscosity — 2026-09-25

## Source and scope

Code `e236810`, parent `7edaf81`. The MPI EB halo API optionally returns
owner-local volume-fraction fields and fetches complete fluid state/temperature
donors, including interior, edges and corners. Covered donors stay zero.
Geometry remains replicated; state requests use the existing owner exchange.
Allocation/count bounds account for the complete halo, and the request mode
is collectively agreed before exchange. Temporal ancestor donor selection
uses the existing interpolation path. Physical coordinates clamp to the domain;
multiple active inlet corners choose the first axis, matching serial policy.

MPI common-step hydro accepts finite nonnegative explicit coefficients,
collectively checks exact equality, and supplies complete halos to positive
coefficient Euler updates. Zero/omitted retains the previous face-only path.
Existing peer reconciliation, restriction/reflux and transaction boundaries
remain. Recursive/coupled MPI and application/default/restart wiring remain
outstanding.

The serial policy's PCM ancestor volume-fraction approximation is retained.
These tests do not establish true fine-scale ghost geometry or external
PeleC corner-policy agreement. Reference geometry fidelity remains open.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Ninja Debug, Python 3.12.3; reused
`build/spray-amr-reactive-debug`. Focused CTest **2/2 passed**: existing zero
coefficient three-rank hydro 1.68 s; coefficient 0.037 variant 1.83 s;
total 3.55 s. No failures or skips. Compiler warnings remain.

The 4/5-node nested-tree fixture covers x/y/z cut planes, distributed ownership
and all patches on the last rank. Positive-coefficient conserved states,
temperatures and boundary ledgers match serial exactly. Wall ledger and total
inventory retain the existing scaled 3e-14 bounds. A rank-zero coefficient
change to 0.05 rejects collectively with state/temperature/ledger rollback;
existing invalid solver, timestep mismatch and bad density checks also pass.
No tolerances changed. Positive-coefficient temporal halos and Mach corners
were compiled but not independently exercised in this increment.

This is focused MPI/serial hierarchy evidence, not a full matrix, external
solver parity or physical validation. Mechanism provenance, Python syntax,
project contract and whitespace checks passed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_(hydro|viscosity)_3d_3$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/mpi-eb-tree-viscosity-checks/build.log` and copied `LastTest.log`.
