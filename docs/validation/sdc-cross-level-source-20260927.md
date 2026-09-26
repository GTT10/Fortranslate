# Recursive EB SDC cross-level auxiliary sources — 2026-09-27

Implementation `b5abd928`, based on `7c9a2165`, local `main`, version 0.246.0.
PeleC reference `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
Provenance preflight H2O2 source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
The focused unit fixture uses the existing elementary H2O2 mechanism.

## Source routing

Serial/MPI recursive EB SDC accepts optional `route_auxiliary` (default false).
True requires the auxiliary callback/ledger pair and agrees across MPI ranks.
The callback adapter validates off-level fluid rates, ignores solids and coarse
cells hidden by refinement, and retains separate old and final-new fields.
The selected-level operator continues to consume only its active level.

For each accepted interval, off-level increments are `dt/2*(old+final new)`.
They accumulate in private buffers keyed by the shallower of source and target
levels. After that interval's descendants and flux synchronization complete,
the buffer is added to its uncovered target cells, EOS temperatures are
recovered and parent states restricted. The separate auxiliary inventory adds
the integrated off-level increments exactly once. No internal exchange is
added to the external-forcing ledger. A later failure leaves every caller
field and inventory unchanged. Default false preserves selected-level masking.

This is conservative explicit source synchronization. Off-level increments do
not retroactively rerun target chemistry or transport within the just-completed
interval. No temporal-order or upstream scheduling-parity claim is made.
Recursive particle/film population timing and migration remain unimplemented;
this adds the gas-side cross-level feedback path needed for that work.

## Focused test scope

GNU Fortran 13.3.0, Ninja Debug, fixed-mechanism tests-enabled build
`build/spray-amr-reactive-debug`; OpenMPI 4.1.6 via `/usr/bin/mpifort` and
`/usr/bin/mpiexec`. Existing real-comparison and executable-stack warnings remain.

The three-level/five-patch fixture supplies nonzero source profiles on every
level in each callback, with patch/cell/stage/time dependence. Routing combines
all source-level intervals: the exact integrated leaf inventory is
`7.25*outer_dt*stage_zero_profile`, since the three level contributions are
`2`, `2.5` and `2.75` times the outer interval. Hidden coarse duplicates are
excluded. The default test retains its original per-target-level formula.
Both use two ranks, distributed ownership and all patches on rank zero,
chemistry, molecular transport and gas forcing.

The tests compare serial/MPI gas, temperature and reaction history exactly,
auxiliary inventory to its analytic result within 3e-14 relative, and combined
closure within 3e-14 times max(1,abs(initial inventory)). A final-stage rejection
on the last descendant retains gas, temperature, reaction history, all budgets,
limiter and step counts, including when the rejecting MPI rank owns no patches.
Redistribution combined with routing, nonphysical synchronized-source rejection,
selected-mechanism builds, installed binaries, full matrix and physical validation
are outside this focused check.

## Completed results

CTest **2/2 passed**, 139.46 s total: default recursive auxiliary test
69.42 s; routed auxiliary test 69.99 s. Analytic source inventories and
late collective rollback assertions pass. Both ownership layouts have exact
serial/MPI state, temperature and reaction history. Routed normalized composite
residual is **4.0340877679024619e-15**; the unchanged default gives
4.2613413117553738e-15. Provenance, builds, Python syntax, project contract
and diff whitespace checks passed. No build or test failures occurred.

## Reproduction

Artifacts: `build/sdc-cross-level-source-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_sdc_aux_(route|recursive)_3d_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
