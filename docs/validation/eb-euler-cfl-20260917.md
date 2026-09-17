# Plain Euler EB hydro timestep estimate — 2026-09-17

Implementation: `0c76991f638191eb63b5a0fa3073e2f6ad01d9b9`, based on
`2d8f292` on local `main`.

`reactive_eb_cfl_3d_mod` now exposes two additional library APIs:

- `reactive_eb_signal_speeds_3d`: EOS-checked maxima of `abs(u_d)+c` in each
  Cartesian direction and `sum(abs(u))+c` as an orientation-independent wall
  bound. Covered storage is ignored. Empty or invalid active domains fail.
- `compute_reactive_eb_euler_cfl_timestep_3d`: a plain Euler hydro estimate
  using `dt = CFL / max_cell(sum_faces(A * signal)/(kappa * V))`. The sum
  includes both Cartesian faces in each direction and the embedded wall.
  Optional `donor_speeds(4)` raises the local maxima to account for remote
  coarse/fine or rank donors. The caller must provide that bound when donors
  can be faster than the local cells.

The full face sum produces twice the conventional Cartesian acoustic rate in
uniform regular cells. For a thin planar cell, normal open-face and wall terms
scale as `1/kappa`, while tangential open areas scale with `kappa`. The wall
signal bound uses the L1 velocity norm, which bounds its projection along any
unit wall normal. This is a conservative hydro estimate, not a proof of
nonlinear Riemann-solver/reflux positivity or coupled reaction/transport
stability. No automatic application selection or diffusion estimate is added
by this commit. Existing redistributed single-level CFL arithmetic is unchanged.

## Independently executed checks

GNU Fortran 13.3.0, Debug/Ninja, existing
`build/selected-tree-mpi-debug`. The exercised unit is serial.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target test_reactive_eb_cfl_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^unit_reactive_eb_cfl_3d$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, zero failures/skips, 0.04 s (0.07 s total).
Provenance, project contract and whitespace checks passed. No Python changed.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

The expanded unit retains the original redistributed estimate checks and
compares the new rate against analytical regular-cell and planar-slab values
for each Cartesian orientation, with fluid fractions 0.1 and 0.01. Doubling
all remote signal bounds halves the timestep. Covered NaN storage is ignored;
invalid active state, invalid donor bounds, invalid CFL and an entirely
covered domain return failure and zero timestep.

Logs: `build/eb-euler-cfl-checks/{build,ctest}.log`.
No configuration matrix, installed audit, hierarchy stepping, external solver
parity or physical validation was rerun for this increment.
