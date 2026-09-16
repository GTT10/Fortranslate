# MPI physical-boundary ghost cells — 2026-09-17

Source: parent `5012e5c50d6c99aedbe86406c56d5c493922ade2` plus this change.
Hashes and logs: `build/spray-mpi-boundary-ghost-checks/`.
Existing GNU Fortran 13.3.0 Debug / Open MPI 4.1.6 configuration.

`fill_mpi_spray_ghosts` fills two ghost layers for an x-slab decomposition.
Physical faces use the existing serial periodic, inflow, outflow, slip-wall,
no-slip-wall and isothermal-wall rules, including optional viscous wall velocity.
Neighbor exchange includes transverse ghost rows, preserving corner composition.
Configuration mismatch or rejected boundary data leaves prior output arrays intact.
Each local dimension must contain at least two cells. Physical model values and
inlet/wall data remain caller inputs; face labels, dimensions and viscous mode
are checked collectively. Hydro, transport, SGS and application stepping are not
yet connected to this API.

Completed commands:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target test_spray_mpi --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^spray_mpi_migration_[24]$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

2/2 focused MPI entries passed in 10.05 s, with no failures or skips.
New cases compare every conserved field and temperature ghost cell exactly with
serial filling on an uneven x partition: fully periodic, mixed wall/outflow with
tangential wall velocity, and inflow with transverse walls. A face mismatch on
one rank rejects collectively and preserves both output arrays exactly.
Provenance, contract and whitespace checks passed; no Python files changed.
This is boundary-data evidence, not full nonperiodic flow or physical validation.
