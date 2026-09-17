# Multilevel EB coupled common step — 2026-09-18

Implementation: `03611153f4459d9681b2bef0b0e46e19da7b0dcd`, based on
`05464f5`. API: `advance_reactive_eb_tree_full_3d`.

## Contract

The serial API composes exposed-cell reaction half steps, molecular transport
half steps and a full Euler hydro step in R/2 T/2 H T/2 R/2 order. Chemistry
and transport can be disabled independently; explicit/implicit chemistry uses
the existing optional integrator selector. Each operator retains its tree
halo, synchronization and averaging semantics. Physical hydro faces use
outflow; physical diffusion faces are closed.

The candidate tree and all stage diagnostics remain private through the last
operator. Boundary transfer sums hydro and transport exchanges, wall transfer
comes from hydro, chemical change sums exposed reaction inventories, and the
limiter reports the minimum over both transport halves (one when disabled).
Every output publishes together only on success. A later failure leaves the
entire input tree and all caller diagnostics unchanged.

This composition does not establish second-order temporal accuracy: hydro
and transport stages remain Euler. Common-interval small-cell restrictions,
spatial PCM and repeated snapshot costs remain. Multilevel application,
checkpoint/restart, subcycling, MPI, convergence and physical qualification
are still open.

## Completed focused verification

GNU Fortran 13.3.0, `/usr/bin/gfortran`, existing Debug/Ninja build
`build/selected-tree-mpi-debug`; no toolchain configuration change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_eb_tree_full_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_eb_tree_full_3d$' --output-on-failure --output-log build/eb-tree-full-checks/ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, zero failures/skips, 34.33 s (total 34.36 s). The fixture has
a 4-cubed root, two adjacent children and grandchildren belonging to different
parents, two ratio-2 levels and a planar z=0.31 wall. Node-dependent density
scaling and energy perturbation exercise transport and reaction. Eight cases
combine reaction/transport on/off and explicit/implicit chemistry with HLLC
hydro at dt=1e-8, reaction tolerances 2e-7/1e-12.

Composite change minus boundary/wall/chemical exchanges, normalized by
`max(1,abs(initial_component))`, is at most **2.2753363795462715e-15**
(asserted bound 3e-14). Enabled chemistry has nonzero species exchange and
preserves its first five conserved components within the same bound.
Disabled chemistry returns zero exchange and disabled transport returns
limiter one. The fully enabled explicit case matches separately invoked
public operators exactly in every node state/temperature and all diagnostics.
An invalid hydro solver after successful private reaction/transport halves
preserves all inputs and sentinel diagnostics exactly.

Provenance, project contract and whitespace checks pass. No Python source
changed. No build or test failed. Existing integer-kind comparison warnings
remain in transport. Logs: `build/eb-tree-full-checks/{build,ctest}.log`;
fixture residuals: `build/selected-tree-mpi-debug/Testing/Temporary/LastTest.log`.
