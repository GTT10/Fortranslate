# Multilevel EB common-step hydro — 2026-09-18

Implementation: `9558304285c967e39bac60b2123f82729afb10eb`, based on
`c978a28`. Core API: `advance_reactive_eb_tree_hydro_3d`.

## Contract and scope

The serial routine averages a private initial tree, fills face halos from that
same snapshot, advances every node with plain EB Euler hydro and accumulated
Riemann fluxes, then invokes multilevel peer/coarse-fine synchronization and
state averaging. Spatial PCM donors and physical outflow faces follow the
existing tree halo contract. All nodes use one caller-supplied positive interval.

Embedded-wall impulse is integrated only on exposed cells, excluding immediate
child-covered regions. Physical boundary exchange uses synchronized effective
root fluxes. Tree state, temperatures, boundary and wall vectors publish together
only after every operation succeeds. Reported vectors represent this step,
not cumulative history. Finite input interval, output shapes and tree geometry
are checked before work; failed halo, hydro or synchronization leaves outputs
unchanged.

This is a library common-step hydro path. It does not add multilevel transport,
combined reaction splitting, subcycling, application/restart integration, MPI,
StateRedist, interval selection or broader physical boundary models. Plain
Euler retains the small-cell interval restriction. Per-node halo preparation
currently copies and averages the tree; large-tree performance is unqualified.

## Completed verification

GNU Fortran 13.3.0, `/usr/bin/gfortran`, existing Debug/Ninja build
`build/selected-tree-mpi-debug`; no compiler/MPI configuration change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_eb_tree_hydro_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_eb_tree_hydro_3d$' --output-on-failure --output-log build/eb-tree-hydro-checks/ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, zero failures/skips, 0.58 s (total 0.61 s). Six fixtures
combine planar x/y/z cuts at 0.31 with two adjacent children and either one
boundary-touching grandchild or adjacent grandchildren belonging to different
parents. A 4-cubed root and two ratio-2 refinements use node-dependent density
and pressure multipliers. An HLLC step of 1e-8 changes deep states and satisfies
`after - before - boundary - wall`, normalized by `max(1,abs(before))`, to at
most **5.6005960242831442e-15** (asserted bound 3e-14). Invalid solver and
invalid deep density reject with exact state, temperature and diagnostic
preservation. This is conservation/transaction evidence, not convergence or
external solver/physical validation.

Provenance, project contract and whitespace checks pass. No Python source
changed. No failed build/test attempt occurred. Logs:
`build/eb-tree-hydro-checks/{build,ctest}.log`; detailed fixture residuals:
`build/selected-tree-mpi-debug/Testing/Temporary/LastTest.log`.
