# EB tree prescribed inflow stepping — 2026-09-26

## Implementation

Source 3f9ce5e, parent e1f74f1, local main.

The serial EB tree timestep, common and recursive hydro/transport operators,
and full reacting split now forward optional boundary_faces and inflow_state
with inflow_temperature to halo construction. Explicit labels remain limited
to inflow/outflow. CTU receives the same physical face labels for its boundary
stencils. Explicit open transport faces contribute molecular fluxes and the
existing synchronized exterior ledger; absent explicit labels retain legacy
closed physical transport faces, except Mach inlets.

The existing candidate-tree transaction, coarse/fine flux replacement,
average-down and ancestor time interpolation remain in use. Both common and
recursive full stepping forward the controls through each transport half and
the hydro update. Inflow temperatures serve as EOS guesses for the conserved
states. Conflicting Mach and prescribed-state controls reject in the halo.

## Verification commands

GNU Fortran 13.3.0, Ninja Debug, Python 3.12.3, existing tests-enabled
build/spray-amr-reactive-debug. No toolchain/backend changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_subcycle_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_ctu_inflow_(full|common)_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The fixture uses five patches across three levels with CTU and nonzero
artificial viscosity. Each schedule runs explicit/implicit chemistry and all
four chemistry/transport enable combinations. It compares total inventory
against the boundary, embedded-wall and chemistry ledgers at the existing
3e-14 relative tolerance; the coupled split also matches separately called
public operators exactly. Transport energy exchange must be nonzero. Invalid
hydro after private earlier stages and invalid inlet data must roll back.

Artifacts: build/eb-tree-inflow-step-checks/. The initial build failed on a
test line longer than the Fortran limit; the line was split. No numerical
tolerances or golden data changed. Existing compiler warnings remain.
Mechanism SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.

## Results

**2/2 passed, 0 failed, 0 skipped**, 251.56 s total. Recursive stepping took
149.99 s and common stepping 101.51 s. Both passed all eight enable/integrator
combinations, independent split equality and failure rollback checks.
Provenance, Python syntax, project contract and final diff checks passed.
The final source differs from the tested build only in indentation of newly
added continuation lines; numerical statements are unchanged.

## Limits

These are serial library checks. MPI forwarding, tree application input and
restart identity remain pending. Periodic/wall tree dispatch is not added.
The timestep API is wired but has no separate inlet-speed scaling check here.
No full configuration matrix, installed audit, external solver parity or
physical validation was run. Full PeleC parity remains unproven.
