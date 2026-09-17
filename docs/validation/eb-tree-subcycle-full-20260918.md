# Coupled subcycled EB tree application — 2026-09-18

Implementation: `8ce55558c91edfa21998db0d4545f2568bf18b2b`, based on
`8e9bf25`. The full tree step accepts optional `subcycling`; the shared
application exposes the corresponding namelist switch, default false.

R/2 T/2 H T/2 R/2 now optionally uses recursive transport (ratio-squared fine
steps) and hydro (ratio fine steps). Reaction halves retain exposed composite
masking and explicit/implicit integration. Boundary/wall/chemical exchanges
and the minimum transport limiter publish with the complete candidate tree.
Failure after earlier stages preserves all caller inputs.

Subcycled checkpoints add `PELEF_EB_TREE_SUBCYCLE_V1` ahead of the existing
static/scheduled/adaptive policy envelope. Readers reject either direction
of stepping-mode mismatch. Common-step records retain their existing layout.
The automatic estimator currently uses its conservative common all-node
bound in either mode; refinement-scaled coarse estimates remain future work.
No higher-order temporal accuracy or MPI support is claimed.

## Completed focused verification

GNU Fortran 13.3.0, `/usr/bin/gfortran`, existing Debug/Ninja
`build/selected-tree-mpi-debug`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/check_eb_tree_application_3d.py
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_eb_tree_subcycle_full_3d pelef_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(unit_eb_tree_subcycle_full_3d|regression_eb_tree_subcycle_3d)$' --output-on-failure --output-log build/eb-tree-subcycle-full-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_eb_tree_application_3d$' --output-on-failure --output-log build/eb-tree-subcycle-full-checks/common-ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

Subcycled **2/2 passed**, zero failures/skips, total 74.78 s (unit 46.44 s,
application 28.31 s). Common application **1/1 passed**, 2.50 s (total 2.53 s).
These are three distinct focused tests, not a full configuration matrix.

The coupled unit covers reaction/transport on/off with explicit/implicit
chemistry on a five-node, three-level tree. Maximum normalized inventory
residual is **2.8441567677384985e-15**, below the unchanged 3e-14 bound.
The fully enabled explicit case matches separately invoked operators exactly
in all states, temperatures and diagnostics. A hydro rejection after successful
reaction/transport halves preserves state and sentinel diagnostics.

The scheduled application changes 5 nodes -> root only -> 4 nodes at 1.5e-9
and 2.5e-9. Output has 216 exposed cells; maximum normalized budget residual is
**5.94704866236712754e-15**. Checkpoints before and after either event continue
to byte-identical final state and budget CSVs. Independent volume integrals
and policy/stage/corruption rejection checks pass, including a subcycling-mode
mismatch. Selected-mechanism subcycling and adaptive subcycling were not rerun.

Initial compilation rejected an overlong application call; line continuation
fixed it. No test failed. Provenance, Python syntax, project contract and
whitespace checks pass. Logs are in `build/eb-tree-subcycle-full-checks/`,
including `initial-build.log`. Runtime artifacts:
`build/selected-tree-mpi-debug/tests/eb_tree_subcycle_3d/`.
