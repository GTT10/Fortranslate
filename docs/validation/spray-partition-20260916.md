# Cartesian particle ownership and send ordering — 2026-09-16

Tested parent: `99d4299` plus this change. Exact source/test hashes and build/test
logs are in `build/spray-partition-checks/`.

The new partition module maps positions to Cartesian cell blocks, giving remainder
cells to low coordinates and using x-fastest rank numbering. Periodic positions are
wrapped only for ownership. Stable grouping produces counts and zero-based offsets
for send buffers, preserving entire parcel records and component columns. Inputs
must have increasing IDs and valid in-domain physical positions. Invalid input
preserves all caller arrays. This is MPI preparation; no MPI communication, remote
gas feedback or distributed run is claimed.

GNU Fortran 13.3.0 Debug, MPI disabled, bounds/FPE checks enabled:

```sh
cmake --build build/port-spray-debug --target test_spray_core --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_partition$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

One focused entry passed, no failures/skips. It checks unequal cell blocks, periodic
ownership, four-rank send counts/offsets with empty ranks, stable within-rank IDs,
exact KH inventory and component-column preservation, and invalid-owner rollback.
Provenance, project contract and whitespace checks passed. No Python source changed,
MPI executable, application, complete matrix or physical benchmark was run.
