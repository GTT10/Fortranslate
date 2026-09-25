# SDC hierarchy transport stage — 2026-09-26

Implementation `2bc38854`, parent `e5dc0f6b`, local main, version 0.246.0.

`compute_sdc_amr_patch_tree_transport_3d` returns transport divergence and
dt-integrated face fluxes over the whole tree without updating state or
refluxing. It reuses existing viscous ghost sampling with MC coarse
interpolation, pressure treatment, physical-face detection and molecular
transport fluxes. Wall/inlet data are forwarded. All outputs, including the
minimum transport theta, are committed only after every patch succeeds.

## Completed checks

- GNU Fortran 13.3.0 Debug fixed tests-enabled build
  `build/spray-amr-reactive-debug`, target `test_amr_patch_tree_transport_3d`, passed.
- `unit_amr_patch_tree_transport_3d`: **1/1 passed**, 33.59 s test / 33.64 s total.
  Added nonuniform root-only and branched multilevel periodic cases. The
  root-only rate reconstructs the existing single-grid transport Euler update
  within 1e-12 scaled state error. Applying tree reflux to the returned fluxes
  closes composite component inventories within 1e-12 scaled error. Invalid
  timestep preserves rates, all face fluxes and theta. Existing transport
  boundary/timestep tests remain active.
- Provenance, Python syntax, project contract and whitespace passed. No
  failed/skipped tests; existing compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_amr_patch_tree_transport_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_amr_patch_tree_transport_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/sdc-tree-transport-checks/` contains build, provenance, CTest
and LastTest logs. This exposes the transport stage for old/new SDC evaluation;
the complete hierarchy loop combining hydro, transport, forcing, reactions and
final flux correction remains open. No MPI or physical-validation claim.
