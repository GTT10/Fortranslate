# Sparse MPI AMR physical ghost layers — 2026-09-17

Source: parent `42323f9` plus this change. Source hashes and logs:
`build/amr-mpi-boundary-halo-checks/`. Existing GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration, launched through `/usr/bin/mpiexec.openmpi`.

`build_mpi_amr_sparse_boundary_halos_3d` adds two ghost layers on all axes for
rank-local coarse slabs. The existing two-depth neighbor exchange first provides
real x-neighbor columns, including one-column slabs. Transverse ghost filling
uses the serial boundary rules. At a physical x boundary, a two-column physical
window produces the correct reflected/outflow/reservoir extensions. These replace
all out-of-domain columns, including the second ghost column on the next inward
rank. This avoids a whole-level gather and supports one coarse column per rank.

Boundary kinds, optional-model presence, viscous-wall mode and numeric reservoir /
wall data must agree across ranks. Local dimensions and wall velocities are
validated collectively. Candidate output arrays are moved into the caller only
after every rank succeeds; failures preserve previously allocated ghost arrays.
The existing periodic API and its output shape are unchanged.

The MPI AMR unit program now compares every physical ghost entry, including
corners, exactly against the global serial fill for slip-wall/outflow,
tangentially moving no-slip walls with transverse inflow, and x-inflow with
transverse walls. One-rank nonfinite reservoir metadata and inconsistent face
names reject collectively and retain both ghost arrays exactly. The 8-cubed
coarse grid is exercised on four and eight ranks; the latter has one column per
rank. Existing periodic hydro, transport, sparse ownership and rollback checks
in the same executable also run.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_mpi_amr_reactive_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^mpi_amr_reactive_3d_np(4|8)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 2/2 tests passed (four ranks 41.95 s, eight ranks 34.89 s; total
76.87 s), with no skips. Mechanism provenance, project contract and diff checks
passed. No Python source changed. This is a prerequisite communication API: MPI AMR
physical hydro/transport stepping, coarse-fine interpolation and application
activation remain to be connected. No MPI physical-flow or scaling qualification
is claimed from ghost-layer equality.
