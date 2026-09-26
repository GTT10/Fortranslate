# Composite-tree weighted state redistribution — 2026-09-27

Implementation `cfb21f99`, based on `c5457b81` on local `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

The serial library now exposes `redistribute_eb_tree_state_3d`. It applies
zeroth-order weighted state redistribution to composite fluid leaves, using
physical fluid volumes. Hidden parent values are retained without participating
in neighborhoods; solid outputs are zero. The target for a small leaf is the
specified fraction (default 0.5) of that leaf's Cartesian cell volume.

An open outward receiver is preferred according to the embedded-boundary normal,
as in the single-patch implementation. A receiver must meet both its own target
and the source's target physical volume, so it does not need its own merge.
Overlapping receivers reserve the same membership-weighted shares as the existing
single-grid operator. If an outward receiver cannot supply the needed volume,
deterministic connected neighborhoods grow breadth-first through open faces.
Membership accounting partitions each member's volume between neighborhoods in
both averaging passes. Insufficient connected fluid volume rejects the operation.

Integer geometry searches cross peer/coarse-fine boundaries and implement paired
periodic coordinate wrapping. Boundary labels and root periodic metrics are
validated. A candidate is published only after successful finite calculation.
No gas EOS recovery, SDC dispatch or MPI wrapper is connected in this increment.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`; existing compiler warnings remain.

- Mechanism provenance and the new test target built successfully.
- `unit_eb_tree_state_redistribution_3d`: 1/1 passed,
  0.02 s test / 0.06 s total.
- A plane crossing the x-directed patch interface exercises root-only,
  left-refined, fully peer-refined and right-refined layouts. Root-only and
  fully refined peer layouts agree with the existing single-grid weighted
  StateRedist within 2e-14 absolute component discrepancy.
- Partial refinement checks transfers across coarse/fine interfaces. Fine small
  cells can use coarse receivers; coarse small cells requiring several fine
  cells exercise connected fallback. Both signed components conserve physical
  fluid-volume integrals. Relative residuals across the four layouts are
  2.2204460492503131e-16, 6.7709024580337711e-16, 1.2212453270876722e-15 and
  1.3322676295501878e-15.
- Composite outputs remain within the input component minima/maxima within
  2e-14, and constants are preserved within 2e-14. Solid input NaNs are ignored.
- Invalid target, nonfinite fluid input and absent fields reject; rejected
  calls retain the previous output.
- Python syntax, project contract and `git diff --check` passed. No failures,
  skips or interrupted runs were counted.

Artifacts: `build/eb-tree-state-redist-checks/`, containing `provenance.log`,
`build.log`, `ctest.log`, `detail.log` and `contract.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug \
  --target test_eb_tree_state_redistribution_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug \
  -R '^unit_eb_tree_state_redistribution_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The fallback target is raw accumulated fluid volume, not the overlap-adjusted
neighborhood volume; this does not establish a general small-cell CFL bound.
Periodic state-redistribution trajectories, disconnected insufficient-volume
fixtures, MPI ownership, SDC integration, temporal convergence, selected builds,
installed audits, scaling and physical validation were not checked here.
This is not upstream AMR StateRedist parity or higher-order reconstruction.
