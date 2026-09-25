# Periodic composite EB particle stencil — 2026-09-26

Source `94ac088`, parent `b34b0b2`, local main. `eb_tree_spray_stencil` adds
optional periodic axes. Composite cloud construction wraps donors, while
root face metrics and opposing nodal level-set values must match. Positions
are wrapped into the root domain. Visibility to each donor's nearest periodic
centroid image is checked using segment pieces inside the domain, wrapping
each crossing before continuing. Fluid-volume weighting and transactional
output publication remain unchanged.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug, Python 3.12.3; existing tests-enabled build.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_tree_stencil --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_tree_stencil$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **1/1 passed**, 0.02 s test / 0.08 s total, no skips. Checks include
normalized periodic corner support spanning coarse/fine cells, invariance
under domain-length translations, incompatible periodic normal rejection,
and unchanged outputs on mismatched nodal endpoints. A central solid slab
confirms that visibility follows the short path through the periodic seam
instead of the long segment through the solid. Existing nonperiodic cut-cell
volume, constant-preservation and solid-point rollback checks remain active.

Initial compilation exposed initialization inserted into the preceding helper;
it was moved into the intended routine. The initial combined slab/plane test
geometry was rejected by the geometry builder; the final fixture uses a valid
single central slab. Both failures are retained under the artifact directory.
No numerical tolerance or geometry validation was loosened.
Provenance, Python syntax, project contract and whitespace checks passed.

Artifacts under `build/eb-tree-periodic-stencil-checks/`: `build.log`,
`LastTest.log`, `initial-build-failure.log`, `initial-slab-failure.log`.
Mechanism SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

MPI field exchange, particle stepping and periodic MPI spray input remain
unconnected. No full matrix, installed audit, external solver parity or
physical validation was run.
