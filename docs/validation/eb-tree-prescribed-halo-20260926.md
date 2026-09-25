# EB tree prescribed inflow halo — 2026-09-26

## Implementation

Source 7ae179f, parent 9827cbc, local main. The serial
fill_reactive_eb_tree_halo_3d API adds optional boundary_faces and inflow_state.
Explicit faces currently accept only inflow/outflow. Inflow requires conserved
states for six faces and the existing inflow_temperature argument as EOS
initial guesses. Explicit faces and Mach controls are mutually exclusive.

Both the face-only one-cell halo and the wide volume-fraction halo apply
prescribed states only outside the root domain. Interior peers and ancestors
retain existing sampling and time interpolation. At a wide corner the first
inflow axis wins, matching the existing Mach corner convention. Prescribed
states undergo EOS validation; failed construction preserves gas, temperature
and volume-fraction output arrays. Covered halo cells stay empty. Root and
patch geometry remain caller-owned.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug, Python 3.12.3. Existing tests-enabled
build/spray-amr-reactive-debug. No compiler/backend changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_halo_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_halo_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final result: **1/1 passed**, 0.16 s test / 0.21 s total, no failures or skips.
The focused test checks two inlet faces on a boundary-touching child,
unchanged peer values in an interior patch, six-layer root halos and corner
precedence, exact target volume fractions, invalid/missing inlet rollback,
legacy Mach inlet success and conflicting-mode rejection. Existing parent,
grandparent, peer, restriction and temporal interpolation cases also pass.

An initial CTest filter omitted the unit_ prefix and selected no tests; it is
not counted. The corrected test first passed in 0.20 s total, then passed again
after adding the Mach compatibility assertions. A trailing space found by
diff --check was removed. No numerical tolerances or golden data changed.
Provenance, syntax, project contract and final whitespace checks passed.
Mechanism SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.

Artifacts: build/eb-tree-prescribed-halo-checks/build.log and LastTest.log.
Existing compiler warnings remain.

## Remaining connections

This is the serial halo API only. Prescribed boundary forwarding through tree
hydro/transport/full updates, MPI exchange, application input and checkpoint
identity remains pending. Periodic and wall mapping are not added here.
No full matrix, installed audit, external solver parity or physical validation
was run; full PeleC parity is not established.
