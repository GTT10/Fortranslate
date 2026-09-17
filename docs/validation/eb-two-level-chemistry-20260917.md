# Two-level EB reaction/hydro composition — 2026-09-17

Implementation `67820bb364071668953dd27b8f73d82e8457c10e`, based on
`e91cd02` on local `main`.

## Behavior and scope

`advance_reactive_eb_two_level_strang_3d` composes a reaction half-step,
the existing synchronized two-level hydro step and another reaction half-step.
Chemical solves use only unrefined fluid cells on the coarse level and active
fine cells. Refined parents are supplied by volume-weighted restriction after
each reaction stage. The existing explicit/implicit chemistry selector is
forwarded to the established cell chemistry implementation.

The wrapper returns separate per-step boundary, wall and chemical inventory
changes. Chemical changes are measured over the composite active-cell inventory;
individual species can change while mass, momentum, total energy and elemental
inventories retain their existing chemistry conservation contract. All four
state arrays and all three inventory vectors commit together. Failure in either
reaction half, hydro, geometry validation or averaging leaves caller data intact.
Disabled chemistry delegates directly to the existing hydro API and returns zero
chemical change; unused reaction controls are ignored in this mode.

The reaction composition is symmetric, but the underlying Euler hydro step
remains first order: this is not an overall second-order accuracy claim.
Support remains one static serial fine patch with common timesteps and outflow
physical boundaries. Molecular transport, subcycling, multiple patches,
regridding, application/checkpoint and MPI EB hierarchy paths remain open.
Small-cell stabilization and timestep selection are unchanged.

## Independently completed checks

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6, Ninja; existing matching
`build/selected-tree-mpi-debug`, MPI enabled and spray disabled. Core and selected
EB runtime compiled. Execution uses fixed elementary chemistry; no MPI or
application matrix was run. No Python source changed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 2 --target test_eb_two_level_chemistry_3d test_eb_two_level_3d pelef_selected_reactive_eb_3d_runtime
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(unit_eb_two_level_3d|unit_eb_two_level_chemistry_3d)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 PASS, 0 failures, 0 skips, total 3.18 s**:

- Existing two-level hydro: 0.31 s.
- Reaction/hydro: 2.85 s. A cut-plane geometry with a 4-cubed coarse domain,
  ratio two and 4-cubed fine patch advances 1e-8 s with each explicit/implicit
  integrator. The reaction response is nonzero, temperatures remain positive,
  and every component satisfies the composite boundary/wall/chemical ledger.
  Chemical mass, momentum and energy changes are near zero, and H/O/N atom
  changes satisfy the element-conservation check. Deliberately invalid refined
  parent states are ignored by chemical solves and replaced by restriction.
  Disabled chemistry is bit-exact with the hydro-only step. Invalid chemistry
  selection and hydro failure after a successful first reaction half preserve
  all state and ledger outputs. Covered coarse storage is unchanged.

Provenance and source contract passed. Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Logs: `build/eb-two-level-chemistry-checks/{build,ctest}.log`.
This focused result is not convergence, external PeleC parity, detailed-fuel
validation or completion of general EB-AMR support.
