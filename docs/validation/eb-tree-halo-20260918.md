# Multilevel EB face halos — 2026-09-18

Implementation: `9f451f76be6697565e72b36dc9bc0872a9110432`, based on
`2762eff`. Adds `fill_reactive_eb_tree_halo_3d` to the core EB tree module.

## Contract

The API fills one layer of face ghost cells from a private, bottom-up-averaged
snapshot of the tree. Exact integer tree origins identify donors. Same-level
patches take precedence, including peers with different parents. Where no
same-level patch exists, successively coarser levels supply PCM values. This
supports a fine patch touching its parent's boundary without mistaking that
internal edge for a physical boundary.

At a physical domain edge, outflow uses the target's adjacent interior cell.
Closed EB faces and unused edge/corner entries remain zero. Open faces require
a fluid donor with finite positive temperature and valid EOS recovery.
Same-level peer overlays also check matching apertures and face centroids.
Source tree fields remain unchanged; output ghosts publish only after the
entire operation succeeds.

This is spatial snapshot filling with outflow boundaries, not temporal donor
interpolation, general physical boundary conditions or complete multilevel
hydro/transport evolution. Each call currently copies and averages the tree;
performance of larger trees is unqualified.

## Completed focused verification

GNU Fortran 13.3.0, `/usr/bin/gfortran`, Debug/Ninja, existing
`build/selected-tree-mpi-debug`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_eb_tree_halo_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_eb_tree_halo_3d$' --output-on-failure --output-log build/eb-tree-halo-checks/ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **1/1 passed**, 0.04 s (total 0.07 s), zero failures/skips. The fixture
has a 4-cubed root, two adjacent children and adjacent grandchildren belonging
to different parents, two ratio-2 refinement levels and a planar wall z=0.31.
Distinct density/pressure multipliers make each donor source identifiable.

Direct same-level, parent and grandparent donor copies match expected conserved
vectors exactly. Physical outflow, zero solid-face/corner ghosts, unchanged
source snapshot, invalid donor rejection and invalid target rollback pass.
A root boundary cell covered by children uses their restricted inventory.
Its comparison to the analytic uniform value initially used inappropriate
bitwise equality: eight-term restriction produced relative roundoff
1.2731531333948682e-16. The new analytic assertion now uses an eight-epsilon
summation bound; direct donor-copy assertions remain exact. No solver tolerance
or numerical update was changed. The diagnostic failure is retained in
`build/eb-tree-halo-checks/roundoff.log`.

Provenance, project contract and whitespace checks passed. No Python files
changed. Build/final test logs: `build/eb-tree-halo-checks/build.log` and
`ctest.log`. Multilevel flux synchronization, hydro/transport stepping,
application/restart wiring, MPI, external parity and physical validation remain
open; no full matrix or installed audit was rerun.
