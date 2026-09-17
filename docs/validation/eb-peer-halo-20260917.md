# Same-level EB peer face halos — 2026-09-17

Implementation: `abff9ff11dbe7c414a8b42da55a9225b1562544e`, based on
`7c56096` on local `main`.

`overlay_reactive_eb_peer_halo_3d` overlays one same-resolution peer onto an
existing face halo. `peer_offset` is its origin minus the target origin in
cell indices. The pair must have matching spacing, aligned physical origins
and nonoverlapping interiors. The operation compares shared-face open fractions
and tangential centroids, recovers EOS temperature for active peer donors,
and replaces only shared open face slots. Other halos, interiors and
edges/corners retain their previous values. The entire pair operation stages
both conserved and temperature halos before publication.

This API is a building block for multiple EB patches. It does not itself
validate an entire patch collection, choose donor precedence, evolve multiple
patches, reconcile shared fluxes or exchange MPI messages. Callers still own
collection-wide topology validation and a shared-face flux policy, particularly
for viscous transverse gradients. The existing coarse-donor halo API is
unchanged.

## Independently executed checks

GNU Fortran 13.3.0, Debug/Ninja, existing
`build/selected-tree-mpi-debug`; the unit is serial.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target test_eb_patch_halo_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^unit_eb_patch_halo_3d$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, zero failures/skips, 0.07 s (0.10 s total).
Provenance, contract and whitespace checks passed; no Python changed.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

The expanded halo test retains existing coarse/fine cases and exercises all
six orientations with partial face overlap, position-dependent peer state,
open/cut/covered face slots and sentinel values elsewhere. It verifies exact
conserved donor copies, EOS-consistent temperatures and preservation of all
unrelated slots. A bad active donor, overlapping boxes and individually valid
but inconsistent shared geometry reject without modifying either halo array.
The geometry mismatch cases restore a valid peer state before rejection, so
they do not rely on the preceding invalid-density fixture.

Logs: `build/eb-peer-halo-checks/{build,ctest}.log`.
No multi-patch time integration, broader configuration matrix, installed audit,
external solver parity, convergence study or physical validation was rerun.
