# Two-level multiple-patch EB hydro — 2026-09-17

Implementation: `4630d5ae5623963709a1243a6efad0c48c519624`, based on
`5e19769` on local `main`.

`eb_reactive_patch_3d` groups geometry and allocated conserved/temperature fields.
`advance_reactive_eb_multipatch_euler_3d` advances a coarse domain and one or more
nonoverlapping fine patches at a common Euler interval. All patches use the
same spatial refinement ratio and inclusive coarse-grid lower indices.

The routine validates each coarsened patch against its coarse region, masks
the union of refined parents, and averages all initial fine states privately.
Coarse and fine updates use that common snapshot. Coarse-donor face halos are
overlaid by same-level peers; shared fluxes are reconciled using lexicographic
box-origin ownership. Restricted fine fluxes replace the corresponding coarse
faces before one coarse reflux. Final fine states are averaged into all refined
parents. Composite physical-boundary and embedded-wall changes are accumulated
in canonical patch order. All coarse/fine state and both ledgers publish only
after complete success. Patch array order does not define ownership.

Scope: serial, two levels, static boxes, common-step Euler hydro and physical
outflow. This does not connect multiple-patch reaction/transport, subcycling,
application/restart, adaptive regridding or MPI. It is not a convergence or
physical-validation result.

## Independently executed checks

GNU Fortran 13.3.0, Debug/Ninja, existing
`build/selected-tree-mpi-debug`; the new unit is serial.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target test_eb_multipatch_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^unit_eb_multipatch_3d$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **1/1 passed**, zero failures/skips, 0.84 s (0.87 s total).
Provenance, contract and whitespace checks passed; no Python changed.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

Ratios 2 and 3 each exercise two patches sharing a complete face, a partial
face, and separated interiors meeting domain boundaries. Three 1e-8 s updates
check every composite conserved component against boundary/wall changes,
nonzero fine evolution, and each refined parent's volume-weighted agreement
with its children. Reversing patch-array order reproduces coarse/fine states,
temperatures and ledgers exactly. A rectangular fine region split into two
patches agrees with the existing one-patch update within `3e-12 * max(1,abs(value))`
for conserved states and ledgers (also coarse temperatures). Overlapping boxes
and a nonphysical large interval reject without changing any patch or ledger.
No late-fine-step failure is specifically injected.

Logs: `build/eb-multipatch-hydro-checks/{build,ctest,final-build,final-ctest}.log`.
The final log owns the counts; the earlier successful run preceded the merged
patch comparison and explicit nonzero-evolution assertion. No broader matrix,
installed audit, external solver parity, convergence study or physical
validation was rerun.
