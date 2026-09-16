# Sparse MPI two-level physical hydro — 2026-09-17

Source: parent `cc5ba37` plus this change. Source hashes and logs:
`build/amr-mpi-physical-hierarchy-checks/`. Reused GNU Fortran 13.3.0 Debug /
Open MPI 4.1.6 configuration, launcher `/usr/bin/mpiexec.openmpi`.

`advance_mpi_amr_sparse_hydro_3d` now accepts optional physical face definitions,
reservoir states and a conserved boundary-change vector. Physical calls use the
coarse boundary solver and retain its mean upper-face fluxes for the existing
reflux operation. Both coarse time endpoints receive full physical ghost layers.
Fine subcycling then interpolates through actual global indices, including
physical ghost indices, instead of periodic modulo indices. The private
interpolation routines distinguish padded physical arrays from the original
transversely unpadded periodic arrays; default callers retain their old path.

Fine-interface flux accumulation, reflux and average-down retain sparse ownership
and collective checks. Boundary transfer is committed only with both hierarchy
levels after successful synchronization. Optional physical-path presence is part
of the common algorithm contract so ranks cannot enter different communication
branches. Invalid input leaves both levels and temperatures unchanged, with zero
boundary ledger and reflux diagnostic.

The new test places a ratio-two patch over coarse cells 2:3 in all directions on
an 8-cubed grid. This exercises interpolation one coarse cell from x/y/z physical
faces. Slip-wall/outflow, transverse reservoir inlet and no-slip faces are used
with both PCM and characteristic PLM. It checks exact serial equality of coarse
and fine conserved states, both temperature arrays and maximum reflux correction;
composite conservation against the boundary ledger; and collective rollback for
an invalid inlet on the last rank. Four/eight ranks include ranks without fine
cells; eight ranks also gives one coarse column per rank.

Boundary-ledger comparison uses relative 2e-11 with scale floor 1e-20; composite
balance uses 2e-11 times max(1, absolute initial inventory). Existing periodic
MPI tests, prior physical ghost tests and coarse-solver comparisons also run.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_mpi_amr_reactive_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^mpi_amr_reactive_3d_np(4|8)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 2/2 tests passed (four ranks 47.40 s, eight ranks 39.74 s; total
87.17 s), with no skips. Provenance, project contract and diff checks passed.

The first test build caught an omitted species argument in the new scatter call;
this was corrected. No existing numerical tolerance or Python source changed.
MPI physical transport, chemistry/full-split forwarding and application activation
remain open. This is a strictly interior two-level numerical comparison, not
physical validation or a general multipatch/multilevel AMR claim.
