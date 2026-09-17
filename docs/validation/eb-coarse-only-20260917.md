# Coarse-only EB update after derefinement — 2026-09-17

Implementation: `ff56e7bf59c456f3baa27f4eabaff7e4960177ee`, based on
`189a059`. The multipatch hydro, transport and full composition APIs accept
zero-length patch arrays and a 3-by-0 lower-bound array. Existing nonempty
patch paths retain their prior operations.

A coarse-only state advances once for the requested physical interval, even
when fine subcycling is enabled. There are no fine fluxes to replace, so the
reflux/EOS recovery pass is omitted. Physical boundary/wall exchanges and
chemical inventory changes remain available. Full-step failure remains
transactional for state, temperature, exchanges and limiter diagnostics.
A positive refinement-ratio argument is still required for API consistency.

## Independently completed checks

GNU Fortran 13.3.0, `/usr/bin/gfortran`, Debug/Ninja, existing
`build/selected-tree-mpi-debug`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_eb_coarse_only_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_eb_coarse_only_3d$' --output-on-failure --output-log build/eb-coarse-only-checks/ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, 12.43 s (total 12.45 s), zero failures/skips. The unit uses
a 4-cubed planar EB mesh, wall z=0.31, density-sheet initialization, interval
1e-8 s and ratio argument 3. Sixteen combinations cover explicit/implicit
chemistry, enabled/disabled reaction and transport, and common/subcycled modes.

An independently composed sequence of chemistry, patch transport and plain
hydro kernels provides the reference: U, T and minimum transport limiter are
exactly equal. Composite component budgets close within
`2e-14 * max(1, abs(initial amount))`; chemistry-enabled cases have nonzero
species exchange. An invalid hydro solver rejects after any requested first
reaction/transport halves while preserving original fields and sentinel
boundary/wall/chemical/limiter diagnostics.

Provenance, contract and whitespace checks passed. No Python files changed.
Logs: `build/eb-coarse-only-checks/build.log` and `ctest.log`.
The only source edit after the test was an explanatory API comment.

This enables the library endpoint needed when temperature tagging returns no
hot regions. Application zero-patch topology, checkpoint/output handling and
temperature-policy dispatch remain open. This is not a full matrix, installed
binary audit, external PeleC parity or physical validation.
