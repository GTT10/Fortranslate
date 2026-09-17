# Multilevel EB molecular transport — 2026-09-18

Implementation: `d29818ddd7fd351315885c9332feb7eec98185a5`, based on
`7ea51b7`. API: `advance_reactive_eb_tree_transport_3d`.

The serial common-interval API connects existing viscosity, heat conduction,
species diffusion and barodiffusion to tree snapshot halos, limited node fluxes,
peer/coarse-fine synchronization and average-down. Integer global origins
identify physical boundaries, which are closed to molecular transport.
Internal parent edges remain donor interfaces even without a same-level peer.

State, temperature, boundary exchange and minimum local limiter value publish
only on success. The returned limiter is the minimum from node transport
updates; synchronized corrections still require EOS acceptance and can reject
the transaction. This does not introduce a new global positivity limiter,
small-cell redistribution, interval estimator, subcycling, combined reaction
splitting, application/restart integration or MPI support. Halo preparation
still copies and averages a tree for each node; scaling is unqualified.

## Completed focused verification

GNU Fortran 13.3.0, `/usr/bin/gfortran`, existing Debug/Ninja build
`build/selected-tree-mpi-debug`; no toolchain configuration change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_eb_tree_transport_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_eb_tree_transport_3d$' --output-on-failure --output-log build/eb-tree-transport-checks/ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, zero failures/skips, 0.64 s (total 0.67 s). Six fixtures
combine x/y/z planar cuts at 0.31 and four/five-node trees, including a
one-sided refined parent edge and grandchildren with different parents.
The root is 4-cubed with two ratio-2 refinements. Node-dependent state scaling
and a 10-percent energy perturbation exercise actual deep transport; EOS
recovers temperatures before stepping. All four transport options are enabled
for a 1e-8 interval.

Composite `after - before - boundary`, normalized by `max(1,abs(before))`,
is at most **5.2293019464126058e-15** (asserted bound 3e-14). External
transport exchange is exactly zero. Deep states change and the limiter lies
in [0,1]. Invalid interval and deep negative density preserve state,
temperatures, exchange and limiter exactly. These checks establish bounded
conservation/transaction behavior, not convergence or physical validation.

Provenance, project contract and whitespace checks pass. No Python source
changed. No failed build/test attempt occurred. Build emitted integer-kind
conversion warnings for global face comparisons; no compilation errors.
Logs: `build/eb-tree-transport-checks/{build,ctest}.log`; fixture output:
`build/selected-tree-mpi-debug/Testing/Temporary/LastTest.log`.
