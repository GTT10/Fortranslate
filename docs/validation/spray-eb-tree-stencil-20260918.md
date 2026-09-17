# Composite EB particle stencil — 2026-09-18

Implementation: `17565d8f415788e12662a5feca6ed5dcde4bb1d1`, based on `9f16384`.

The stencil wraps existing composite rectangular-cloud overlap, excludes hidden
parents, multiplies weights and cell volumes by EB fluid fractions, and removes
cells whose fluid centroid is blocked by the original nodal solid. Surviving
weights normalize to one. The particle must lie in fluid; geometry shapes,
spacing and patch bounds must match topology. Physical bounds use a 64-epsilon
roundoff allowance relative to coordinate/cell scale. All output allocations
publish only after the full stencil succeeds.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; existing fixed tests-enabled
and selected H2/O2 tests-disabled configurations reused unchanged.
Compiler `/usr/bin/gfortran`, MPI wrapper `/usr/bin/mpifort`.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_tree_stencil --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_tree_stencil$' --output-on-failure --output-log build/spray-eb-tree-stencil-checks/ctest-final.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/check_project_contract.py tools/check_mechanism_provenance.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **1/1 passed**, **0.01 s**, total **0.04 s**. Selected library
compilation, provenance, syntax, contract and whitespace checks passed.
No failures or skips. An earlier successful check preceded the physical
patch-coordinate consistency guard; it is not an additional test case.

The fixture has a 4-cubed root, a ratio-two child covering its left half and
nodal plane y=1.25 on the finest grid. A particle at (1.9,1.3,1.5) exercises a
cut coarse/fine interface. Assertions check both levels contribute, hidden
parents do not, weights are positive and sum to one, a constant field remains
constant, and returned physical fluid volumes match 0.75 for coarse cut cells
and 0.0625/0.125 for fine cut/regular cells. A solid query preserves previous
outputs; truncated outer-boundary support also normalizes. Arbitrary nonconvex
occlusion is not separately qualified by this planar fixture.

Artifacts: `build/spray-eb-tree-stencil-checks/` contains `build-final.log`,
`build-selected-final.log`, `ctest-final.log` and preceding successful logs.

This is a serial geometry stencil for subsequent distributed field work.
It uses fluid-fraction weighting and centroid visibility, not exact geometric
cloud intersection or linear-exact interpolation. Source deposition, small-cell
redistribution and MPI sampling are not implemented by this API. Geometry and
nodal data must originate from the same gas hierarchy. No convergence,
production scaling, external parity or physical validation was performed.
