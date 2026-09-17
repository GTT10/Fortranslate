# Recursive multilevel EB molecular transport — 2026-09-18

Implementation: `7747a1198dd4a998018db3d95e95178c3a8ca922`, based on
`0183f27`. API: `advance_reactive_eb_tree_subcycle_transport_3d`.

Hydro and transport now share a private recursive operator. The hydro public
API retains ratio-based substeps. Molecular transport takes ratio-squared
substeps, using time-interpolated ancestor donors, current peers, accepted
limited fluxes, peer reconciliation, accumulated transfers and interval-local
reflux. Squared substep counts are checked for integer overflow. Physical
domain faces remain closed to diffusion; internal patch edges are enabled.

The transport API returns boundary exchange and the minimum limiter over all
node updates/substeps. It publishes state, temperature and diagnostics only
on complete success. Peer/reflux corrections still require EOS acceptance;
the local limiter diagnostic is not a new global positivity guarantee.
Molecular wall transfer remains zero under the existing closed-wall operator.
No chemistry composition, automatic coarse-interval estimator, application
mode, MPI, convergence or physical validation is added here.

## Completed focused verification

GNU Fortran 13.3.0, `/usr/bin/gfortran`, existing Debug/Ninja
`build/selected-tree-mpi-debug`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_eb_tree_subcycle_transport_3d test_eb_tree_subcycle_hydro_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_eb_tree_subcycle_(hydro|transport)_3d$' --output-on-failure --output-log build/eb-tree-subcycle-transport-checks/ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, zero failures/skips, total 6.19 s (hydro 1.41 s, transport
4.75 s). Transport uses six fixtures combining x/y/z cuts at 0.31 and four/five
nodes with two ratio-2 refinement levels. Node-dependent density scaling and
energy perturbations produce nonzero deep transport. All molecular options
are enabled over a 1e-8 root interval, with four substeps per child interval
and sixteen deepest substeps.

Composite change minus boundary exchange has maximum normalized residual
**5.2293019464126058e-15**, below the unchanged 3e-14 bound. Physical diffusion
exchange is exactly zero, deep states change, limiter lies in [0,1], and invalid
interval/deep density preserve states, temperatures and sentinel diagnostics.
Existing recursive hydro conservation and rollback checks also pass after
sharing the recursive control flow.

Provenance, project contract and whitespace checks pass. No Python source
changed. No failed build/test attempt occurred. Logs:
`build/eb-tree-subcycle-transport-checks/{build,ctest}.log`. Fixture output:
`build/selected-tree-mpi-debug/Testing/Temporary/LastTest.log`.
