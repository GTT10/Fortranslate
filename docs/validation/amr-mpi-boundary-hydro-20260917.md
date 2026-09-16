# Sparse MPI AMR coarse physical hydro — 2026-09-17

Source: parent `aa7e5ae` plus this change. Logs and source hashes:
`build/amr-mpi-boundary-hydro-checks/`. Reused GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration; launcher `/usr/bin/mpiexec.openmpi`.

`advance_mpi_amr_sparse_boundary_hydro_3d` advances rank-local coarse slabs with
PCM or characteristic PLM and SSPRK2. The physical ghost API provides two layers
on every axis, including one-column slabs. Each stage computes fluxes from the
padded state and updates physical cells only. Returned mean flux arrays retain
both ends of each axis (zero-based in that axis), ready for subsequent reflux
and boundary accounting. State and temperature commit only after collective
success at both stages. The existing periodic API is unchanged.

The boundary ledger integrates only the global x-end faces plus each slab's
y/z-end faces. An MPI sum produces the same global conserved transfer on every
rank. Internal x interfaces are excluded. The ledger is zero on failure.
Species, geometry, timestep, algorithm and physical-boundary metadata are checked
collectively before candidate states are committed.

The MPI AMR unit test compares both reconstructions with the serial physical
SSPRK2 solver for slip-wall/outflow, transverse reservoir inflow and no-slip-wall
faces. Coarse state, temperature and every returned upper-cell mean face flux
must match exactly. Global surface budgets agree within relative 2e-11 with
absolute scale floor 1e-20, and a nonzero budget is required. A bad inlet density
on rank zero rejects collectively, retaining all local states/temperatures and
returning zero budget. Four and eight ranks cover two and one coarse columns
per rank on the 8-cubed grid. Existing periodic AMR tests also run.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_mpi_amr_reactive_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^mpi_amr_reactive_3d_np(4|8)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 2/2 tests passed (four ranks 44.83 s, eight ranks 37.29 s; total
82.15 s), with no skips. Mechanism provenance, project contract and diff checks
passed.

An initial build caught two misplaced continuation markers in new calls; these
were corrected. No numerical tolerance or Python source changed.

This is the coarse-level distributed solver. Physical coarse-to-fine
interpolation/subcycling, full hierarchy reflux, MPI transport and application
activation remain to be connected. It is not a claim of complete MPI AMR
physical-boundary stepping, performance scaling or physical validation.
