# Boundary-touching AMR flux correction — 2026-09-17

Source: parent `d06978e` plus this change. Logs/source hashes:
`build/amr-boundary-reflux-checks/`. Reused GNU Fortran 13.3.0 Debug build
with bounds and floating-point checks.

`reflux_physical_coarse_3d` accepts the complete coarse face flux arrays,
including each lower domain face, plus fine interface fluxes averaged over
child faces and the coarse interval. Interior interfaces correct the adjacent
uncovered coarse cells. Interfaces on the physical domain surface instead
return the fine-minus-coarse surface transfer, positive into the domain.
Callers add that correction to the coarse boundary ledger. Covered cells are
left for average-down. The maximum correction describes actual coarse-cell
updates, excluding external surface transfers.

The operation uses private candidates and commits state, surface transfer and
maximum correction only after finite checks pass. Invalid input preserves the
state and zeroes output diagnostics. A patch touching a periodic seam is
rejected: correcting its opposite-side coarse neighbors requires a separate
periodic ownership rule. Interior patches retain the original reflux behavior.

The hierarchy unit test uses nonzero coarse fluxes and a manufactured linear
fine flux field, four conserved components and anisotropic spacing. Refinement
ratios two and three cover interior patches, each of six physical faces, lower
and upper corners, and the full domain. It verifies composite conservation,
the independent surface integral, untouched covered cells, exact equality with
the existing interior reflux routine, zero coarse updates for full-domain
refinement, and rollback for NaN fluxes and unsupported periodic seams.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_amr_hierarchy_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^unit_amr_hierarchy_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 1/1 test passed, no skips (0.01 s; CTest total 0.04 s). Provenance,
project contract and diff checks passed. No Python source changed.

This implements the boundary-touching flux-correction operator. Fine physical
ghost construction, composite boundary-ledger integration in stepping and
application activation remain open; application patch restrictions are retained.
This is a discrete conservation check, not physical validation.
