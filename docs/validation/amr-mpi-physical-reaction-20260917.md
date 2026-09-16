# Sparse MPI AMR physical-boundary reaction split — 2026-09-17

Source: parent `3fdcde7` plus this change. Logs and source hashes:
`build/amr-mpi-physical-reaction-checks/`. Existing GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration, launcher `/usr/bin/mpiexec.openmpi`.

`advance_mpi_amr_sparse_strang_3d` forwards optional physical faces and reservoir
states through the hierarchy hydro step and returns its conserved boundary
transfer. The ledger stays private until both chemistry halves and hydro succeed.
Face-option presence joins the collective algorithm contract; a supplied ledger
must have the conserved-state size. Both hierarchy levels, temperatures and
boundary/reflux diagnostics retain transactional failure behavior.

The MPI physical hierarchy fixture now also exercises active reaction splitting
for PCM and characteristic PLM. Its ratio-two patch covers coarse cells 2:3 in
all directions of an 8-cubed domain with slip/no-slip, inflow and outflow faces.
For four/eight ranks, coarse/fine state, temperatures and reflux correction are
compared exactly to the serial reaction/hydro split. Boundary transfer matches
within relative 2e-11 with scale floor 1e-20. Euler inventory and H/O/N balance
include the boundary transfer; nonzero chemical conversion is required. Bad
inlet density on the last rank rejects the hydro middle phase and restores the
whole hierarchy, including the already-computed first chemistry half.
Existing physical hydro, ghost and periodic AMR tests remain in the executable.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_mpi_amr_reactive_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^mpi_amr_reactive_3d_np(4|8)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 2/2 tests passed (four ranks 53.65 s, eight ranks 44.95 s; total
98.63 s), with no skips. Provenance, project contract and diff checks passed.

No Python source or existing numerical tolerance changed. This is bounded
numerical split/recovery evidence. Physical MPI molecular transport, the complete
transport split and application activation remain open; no experimental or
scaling qualification is claimed.
