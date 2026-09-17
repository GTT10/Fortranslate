# Redistributed EB exchange inventories — 2026-09-17

Implementation: `be54f14a584328418efc9210c198ed33acf8bdba`, based on
`6343562` on local `main`.

The flux-redistributed and state-redistributed Euler hydro APIs now accept
optional final `boundary_change` and `wall_change` vectors. Boundary exchange
integrates the actual outflow Riemann faces used by the conservative RHS, with
open face areas and timestep. Wall impulse integrates the same initial-state
slip-wall pressure source. Redistribution moves inventory internally and
does not change the global exchange. Vectors publish only after the update
succeeds; failures retain the caller's diagnostic values. Calls without
diagnostics preserve their previous computation path.

`advance_reactive_eb_strang_3d` and `advance_reactive_eb_full_3d` accept optional
boundary, wall and chemistry vectors after the existing integrator argument.
Chemical changes are volume-integrated immediately around each reaction half.
Hydro exchange uses the state at its actual split stage. The existing full
transport operator has closed physical fluxes, so its global boundary exchange
is zero. Diagnostics publish together only after the complete composition is
accepted. Existing state/output failure semantics and limiter behavior are
unchanged; the new diagnostic vectors themselves are transactional.

This connects the numerical update APIs. Application accumulation, persistent
checkpoint ledgers and replacement of the single-level application's
zero-exchange gate remain open. The earlier density-contrast methanol case has
not yet been rerun with application-level exchange accounting.

## Independently executed checks

GNU Fortran 13.3.0, Debug/Ninja, existing `build/selected-tree-mpi-debug`;
the focused unit is serial.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target test_eb_redistributed_exchange_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^unit_eb_redistributed_exchange_3d$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, zero failures/skips, 0.81 s (0.84 s total).
Provenance, contract and whitespace passed; no Python changed.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

Planar fixtures in all three coordinate directions use nonzero outward
velocity, regular/cut density contrast and both redistribution methods.
Boundary mass matches `-dt * rho * velocity * area`; wall impulse matches
`dt * pressure * area`, with zero wall mass/energy/species exchange.
Every conserved component balances boundary plus wall exchange. Coupled
reaction-only/flux-redistribution and reaction/transport/state-redistribution
cases also close after the chemical contribution and have nonzero species
reaction response. Requesting diagnostics leaves state and temperature exactly
equal to calls without diagnostics. Oversized hydro steps retain both ledgers;
an invalid hydro solver after accepted reaction/transport stages retains all
three coupled diagnostic vectors.

Logs: `build/eb-redistributed-exchange-checks/{build,ctest}.log`.
No application restart, broad matrix, install audit, external parity or physical
validation was run for this API change.
