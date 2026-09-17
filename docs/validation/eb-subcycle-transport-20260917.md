# Two-level EB subcycled transport — 2026-09-17

Implementation: `26a84ff246fdec8b287d095487908021d2634bd6`, based on
`64e4783` on local `main`.

`advance_reactive_eb_two_level_subcycled_transport_3d` advances coarse
molecular transport once and the fine patch `ratio**2` times over the same
interval. Each fine Euler step uses coarse conserved donors interpolated
between the initial and unrefluxed coarse endpoint, with EOS temperature
recovery during halo construction. Physical transport faces remain closed;
only internal coarse/fine patch faces receive donor fluxes.

Accepted, species-limited flux integrals accumulate across all fine substeps.
Final flux restriction, reflux and parent averaging synchronize the hierarchy.
The returned limiter is the minimum of the coarse step and every accepted
fine substep. State, boundary inventory and limiter publish together only on
success. The shared subcycling loop still handles ratio-based hydro steps;
common-step transport/hydro arithmetic is retained. The integer square of the
refinement ratio is checked before use.

This remains a serial one-patch library path. It does not yet connect the full
reaction/transport/hydro composition, application timestep/restart mode,
multiple patches or MPI. It is not an order-of-accuracy or positivity proof.

## Independently executed checks

GNU Fortran 13.3.0, Debug/Ninja, existing
`build/selected-tree-mpi-debug`. The exercised tests are serial.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target test_eb_two_level_transport_3d test_eb_two_level_3d test_eb_two_level_full_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^(unit_eb_two_level_3d|unit_eb_two_level_transport_3d|unit_eb_two_level_full_3d)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

**3/3 passed**, zero failures/skips, 12.55 s total: hydro 0.97 s,
transport 3.11 s, existing full composition 8.45 s.
Provenance, project contract and whitespace checks passed; no Python changed.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

The expanded transport test exercises ratios 2 and 3 (4 and 9 fine steps),
an interior patch and a physical-boundary patch, with nonuniform velocity
and composition. Three coarse intervals of 1e-5 s check all conserved
components against the boundary ledger, parent/child volume-weighted agreement,
positive temperatures, limiter range, exactly zero physical transport exchange
and preserved covered storage. Fine evolution differs from a common-step
transport update of the same interval. Invalid geometry/ratio and a large
nonphysical interval preserve all states, ledger and limiter on rejection.
Common-step tests and the hydro/full-composition regressions remain included.
This does not inject a failure after a successful fine transport substep or
specifically establish limiter activation during a subcycled run.

Logs: `build/eb-subcycle-transport-checks/{build,ctest}.log`.
No full matrix, installed audit, application subcycling, external solver parity,
convergence study or physical validation was rerun.
