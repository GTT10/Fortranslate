# Multiple-patch EB reacting composition — 2026-09-17

Implementation: `257d7dc0175f8708a03b43db0e80c40d8118aa23`, based on
`5b72fea` on local `main`.

`advance_reactive_eb_multipatch_full_3d` composes reaction, molecular transport
and hydro as `R/2 T/2 H T/2 R/2` over static two-level multiple-patch geometry.
Reaction solves only active, unrefined coarse cells and active cells of each
fine patch. Chemical inventory changes use canonical box order, and fine states
are averaged into parents after each reaction half. Transport and hydro retain
peer synchronization and coarse reflux from their existing multiple-patch APIs.
The minimum patch/stage limiter and boundary/wall/chemical inventories publish
with all coarse/fine state only after the complete composition succeeds.

Chemistry and transport can each be disabled. Disabling both directly uses the
hydro path, with zero chemical inventory and unit transport limiter. Supported
reaction integrators are delegated to the existing cell chemistry interface.
The symmetric composition still contains first-order Euler spatial updates;
no second-order accuracy claim is made.

Scope remains serial, static two levels, common interval and existing physical
boundaries. Multiple-patch subcycling, application/checkpoint integration,
MPI, and general viscous transverse-stencil equivalence remain open.

## Independently executed checks

GNU Fortran 13.3.0, Debug/Ninja, existing
`build/selected-tree-mpi-debug`; the new unit is serial.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target test_eb_multipatch_full_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^unit_eb_multipatch_full_3d$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, zero failures/skips, 5.90 s (5.93 s total).
Provenance, contract and whitespace checks passed; no Python changed.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

Two adjacent fine patches with refinement ratio 2 exercise explicit and
implicit chemistry at 1e-8 s. Refined-parent storage begins deliberately
invalid, demonstrating that it does not own reaction solves and is replaced
by restriction. Checks cover all-component composite boundary/wall/chemical
balance, chemical mass/momentum/energy preservation, H/O/N element inventories,
nonzero species response, limiter range and exact permutation invariance of
states, temperatures and all diagnostics. Disabling chemistry and transport
matches direct hydro exactly. An invalid hydro solver after the first accepted
reaction and transport halves leaves every coarse/fine state and diagnostic
at its original value.

Logs: `build/eb-multipatch-full-checks/{build,ctest}.log`.
No broader matrix, installed audit, application integration, external solver
parity, convergence study or physical validation was rerun.
