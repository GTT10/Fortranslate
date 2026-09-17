# Serial multilevel EB application — 2026-09-18

Implementation: `d00dc61fa4e6ec9b41bc4fe205e18feb10a04da8`, based on
`6d40e29`. Executable: `pelef_reactive_eb_tree_3d`.
[Case and input contract](../../cases/reactive_eb_tree_3d/README.md).

## Implemented integration

The application builds static parent-before-child topology and coarsenable
finest-lattice geometry, initializes fixed H2/O2 states, advances coupled
common-interval reaction/transport/hydro with bounded retry halving, checks
composite budgets and exports exposed-cell state and per-component inventory
CSVs. Accepted steps accumulate boundary, wall and chemical exchanges.

Checkpoint files wrap the tree stream record with time-step/integrator policy.
Restoration checks policy, existing node physics/geometry records, topology,
parent averages, composite inventory and EOF before any external output.
Saving stages the complete serialization in a scratch stream before replacing
the destination. Copy-time I/O failure is not filesystem-atomic. Input/output
literal-path alias checks follow the existing two-level application convention;
these checks do not resolve symlink/path aliases.

The executable is included in the install target, but no fresh installed-prefix
audit was run. This is a serial static-tree fixed-mechanism application;
selected frontend, regridding, automatic interval estimation, recursive
subcycling, MPI, scaling, convergence and physical qualification remain open.

## Completed focused verification

GNU Fortran 13.3.0, `/usr/bin/gfortran`, existing Debug/Ninja build
`build/selected-tree-mpi-debug`; no toolchain configuration change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/check_eb_tree_application_3d.py
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target pelef_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_eb_tree_application_3d$' --output-on-failure --output-log build/eb-tree-app-checks/ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, zero failures/skips, 2.35 s (total 2.39 s). Five-node,
three-level ratio-2 geometry with z=0.31 wall runs two 1e-9 HLLC steps with
explicit reaction and all molecular transport options enabled. Output has
**272 exposed fluid cells**, including all three levels. Independent CSV
fluid-volume sums match final inventory entries within 3e-14 scaled tolerance;
maximum normalized balance residual is **3.18097960520324641e-15** against
the existing 5e-10 application gate.

Stopping after step one and restarting yields byte-identical final state and
budget CSVs to uninterrupted execution. A changed time-step policy and a
corrupt final record marker reject while sentinel state CSV, budget CSV and
checkpoint destination remain unchanged. These are focused integration and
restart checks, not broader physics or installed-binary qualification.

Provenance, Python syntax, project contract and whitespace checks pass. No
failed build/test attempt occurred. Logs:
`build/eb-tree-app-checks/{build,ctest}.log`. Cases, outputs, checkpoint files
and per-invocation logs:
`build/selected-tree-mpi-debug/tests/eb_tree_application_3d/`.
