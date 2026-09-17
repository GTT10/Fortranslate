# Multilevel EB ancestor time interpolation — 2026-09-18

Implementation: `90de36a387c9943f89db02060b03917da1b979e5`, based on
`d2bd37b`. Extends `fill_reactive_eb_tree_halo_3d` with optional
`ancestor_start`, `ancestor_end` and `ancestor_fraction` arguments.

All three optional arguments must be supplied together. Endpoint trees must
match the current topology and geometry; one finite fraction in [0,1] is
required for each level below the target (root first). When an ancestor is
the selected spatial donor, its conserved state and temperature guess are
interpolated at that ancestor level's fraction, followed by EOS temperature
recovery. Endpoint fractions copy the corresponding state exactly. Same-level
peers and physical outflow continue using the current snapshot.

Endpoint states represent provisional coarse updates and are not averaged
down from descendants at another time. Only donor cells actually used require
EOS acceptance. No optional arguments retains the previous spatial behavior.
Output ghosts/temperatures publish only after complete success; input trees
are unchanged. This is a temporal donor primitive for recursive stepping,
not an implemented subcycling integrator or a temporal accuracy claim.

## Completed focused verification

GNU Fortran 13.3.0, `/usr/bin/gfortran`, existing Debug/Ninja
`build/selected-tree-mpi-debug`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_eb_tree_halo_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_eb_tree_halo_3d$' --output-on-failure --output-log build/eb-tree-time-halo-checks/ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, zero failures/skips, 0.08 s (total 0.11 s). Existing spatial
halo checks remain. Added checks use the five-node three-level cut-cell tree
with doubled conserved endpoint states: fractions 0, 0.5, 1 give the analytic
root and immediate-parent donors within a four-epsilon arithmetic bound;
independent root/parent fractions 0/1 match exact endpoint vectors. The
same-level peer remains an exact copy of the current fine state. Fractions
outside [0,1], incomplete optional context and invalid selected endpoint
density reject while preserving sentinel ghosts/temperatures.

Provenance, project contract and whitespace checks pass. No Python source
changed. No failed build/test attempt occurred. Logs:
`build/eb-tree-time-halo-checks/{build,ctest}.log`. Recursive advancement,
accumulated substep synchronization and application restart mode remain open.
