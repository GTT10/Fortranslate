# Multiple-patch EB subcycling — 2026-09-17

Implementation: `fdfd5e4a792d1f3b9aa70226f99fd93319128199`, based on
`5f1a849` on local `main`.

The multiple-patch hydro, transport and full reacting APIs now accept a final
optional `subcycle` argument. Hydro uses `ratio` fine steps per coarse step;
transport uses `ratio**2` (with an integer-overflow guard) within its interval.
The coarse update runs once. Fine patches use linearly interpolated coarse
conserved donors between the initial and unrefluxed coarse endpoint; the halo
routine recovers EOS temperature from those donors. Every peer halo uses the
same accepted fine-time snapshot, independent of patch traversal order.

Each fine substep reconciles shared fluxes and receiver states before the next
snapshot. Accepted reconciled fluxes accumulate over the interval and drive a
single final coarse reflux and average-down. Fine wall inventory accumulates
from each substep's old state; the minimum transport limiter includes all
patches and stages before peer correction. The complete update remains a
transaction across all states and diagnostics. Reaction halves still operate
around the full hydro/transport sequence on composite active cells.

Both application frontends accept multiple-patch `subcycling=.true.` and scale
fine timestep bounds by the applicable substep count. Multiple-patch V3/V4
checkpoint headers bind fixed/automatic subcycling; V1/V2 retain common
stepping. Patch records and topology layout are unchanged.

Scope remains serial, static two levels with a common refinement ratio,
outflow hydro and closed physical transport boundaries. This does not add
MPI, regridding, arbitrary-depth EB, small-cell redistribution, higher-order
time accuracy, or general joined/split transverse-viscous stencil equivalence.

## Independently executed checks

GNU Fortran 13.3.0, Debug/Ninja, existing `build/selected-tree-mpi-debug`;
configured `/usr/bin/mpifort` is OpenMPI 4.1.6. These tests are serial.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target test_eb_multipatch_3d test_eb_multipatch_full_3d \
  pelef_reactive_eb_two_level_3d test_selected_full_reactive_eb_two_level_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^(unit_eb_multipatch_3d|unit_eb_multipatch_full_3d|regression_reactive_eb_two_level_3d|regression_selected_eb_two_level_3d)$' \
  --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/check_selected_eb_two_level_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**4/4 passed**, zero failures/skips, 43.09 s total: fixed application 7.68 s,
hydro/transport unit 12.09 s, reacting unit 12.58 s, selected application 10.70 s.
The first build exposed an overlong Fortran declaration; splitting that line
resolved compilation before these tests. Provenance, syntax, contract and
whitespace checks passed. Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

Hydro/transport cases retain common stepping and add subcycling for ratios 2
and 3, full/partial shared faces and separated boundary patches, over three
coarse intervals. Checks include all-component composite budgets, parent
averages, exact permutation invariance, closed transport boundaries, limiter
range and rollback. The rectangular zero-velocity transport fixture also
compares joined/split results with the existing single-patch subcycled API;
the partial-face case exercises nonzero transverse velocities without a
general joined/split equivalence claim.

Coupled tests cover explicit/implicit reaction, nonzero chemistry, H/O/N
inventory, exact permutation invariance, disabled-operator hydro agreement,
and complete rollback when hydro fails after accepted reaction/transport
halves. The selected application adds two-patch subcycling to its ten-species
fixture: three steps to 3e-8 s, no rejections, normalized composite residual
**6.6425382718082362e-16**, minimum limiter 1. Restart is byte-identical in
state and budget CSVs; changed stepping mode, bundle SHA and composition reject
without overwriting the existing outputs.

Logs: `build/eb-multipatch-subcycle-checks/{build,ctest}.log`.
Application artifacts: `build/selected-tree-mpi-debug/tests/selected_eb_two_level_3d/`.
No broader matrix, installed audit, external solver parity, convergence study
or physical validation was rerun.
