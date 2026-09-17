# Connected-small-cell redistribution — 2026-09-17

## Source and scope

Implementation: `1ff853513c6c55a3a8ffeb7eb3a1573926358e5c`, based on
`dc669cac8601ca43f8a1574b1af8db5955faf594` on local `main`.
This record describes independently completed focused checks, not full PeleC
parity, physical validation or a compiler/configuration matrix.

The existing 3D weighted StateRedist path retains its original arithmetic when
all small cells have eligible outward receivers. If any receiver is unavailable,
the entire geometry uses a conservative zeroth-order neighborhood fallback.
For each small cell, breadth-first growth follows positive shared face apertures
in x-, x+, y-, y+, z-, z+ order until the summed volume fractions meet the target.
Large cells own singleton neighborhoods. Covered cells and outer boundaries
cannot be crossed. Insufficient connected volume rejects the operation with
zero redistribution output; callers retain their existing physical-state rollback.

If cell i belongs to N_i neighborhoods, its contribution to each average is
kappa_i / N_i. Each neighborhood average is then shared back with weight 1/N_i.
Summing kappa_i times the output therefore recovers the input weighted sum.
Both passes use positive weights, preserving constants and component bounds.
Queues and membership counters use linear storage; member lists are sparse but
can grow quadratically in a pathological all-small connected component. No
persistent geometry cache is introduced.

The target applies to raw neighborhood volume, not the volume after overlap
sharing. This is not proof of a geometry-independent CFL bound. The fallback
is zeroth order and can smooth resolved gradients. General small-cell stability,
convergence, EB-AMR-MPI integration and scaling remain open.

## Completed checks

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6, Ninja; existing matching build
`build/selected-tree-mpi-debug`, MPI enabled, spray disabled. The fixed and
selected-elementary EB applications were rebuilt. MPI execution was not tested.
Runtime/build commands used `env -u LD_LIBRARY_PATH -u LIBRARY_PATH`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 2 --target test_eb_reactive_redistribution_3d test_reactive_eb_transport_3d pelef_reactive_eb_3d test_selected_elementary_reactive_eb_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 2 --target test_reactive_eb_transport_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(unit_eb_reactive_redistribution_3d|unit_reactive_eb_transport_3d|regression_reactive_eb_level_set_3d)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest: **3/3 PASS, 0 failures, 0 skips, 13.00 s**.

- Redistribution unit: 0.01 s. Existing axis-plane and curved cases plus thin
  sheets in all three orientations with no target-sized receiver; weighted
  conservation, nontrivial averaging, convex bounds and constants. Two fluid
  strips separated by covered cells remain isolated; requesting more volume
  than either connected strip contains rejects with zero output.
- Transport/coupled unit: 0.91 s. Existing tests plus a thin small-cell sheet with
  a constant-pressure density perturbation, hydro alone, hydro/transport and
  hydro/transport/reaction. Positive state, unchanged covered storage, mass,
  energy and element conservation use existing tolerances; separate transport
  and chemistry responses are nonzero.
- Existing nodal level-set application regression: 12.05 s. Fixed/selected
  continuous and checkpoint continuation checks, relocation, changed geometry
  and malformed input rejection. This regression checks compatibility of the
  established sphere case; the new cluster path is exercised by the units.

Provenance preflight and source contract passed; mechanism source SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
No Python source changed; the existing application checker executed successfully.
Build logs retain existing compiler warnings.

Logs: `build/eb-cluster-redist-checks/{build,coupled-build,ctest}.log`.
An initial redistribution-only run also passed (`unit-first.log`); it is not
added to the final test count. Application artifacts:
`build/selected-tree-mpi-debug/tests/reactive_eb_level_set_3d/`.
