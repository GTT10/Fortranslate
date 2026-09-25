# Common-step AMR SDC integration — 2026-09-26

Implementation `127317f8`, parent `46482763`, local main, version 0.246.0.

`advance_amr_patch_tree_sdc_3d` coordinates hydro, old/new molecular transport,
external force and forced chemistry across the complete serial patch tree.
Each iteration predicts hydro from the same interval-start tree, constructs
all candidate patches before evaluating new sources, and integrates reaction
from the interval-start state. The final hydro flux plus averaged transport
flux is refluxed and restricted; temperatures are then recovered. Separate
boundary and composite external-source inventories are returned. Optional
reaction-source estimates persist across calls. All state and diagnostics
remain private until the complete hierarchy step succeeds.

This is a common timestep across levels, PPM/PeleC hydro and MC coarse
interpolation. It does not implement recursive SDC subcycling or MPI dispatch.

## Completed checks

- GNU Fortran 13.3.0 Debug fixed tests-enabled build
  `build/spray-amr-reactive-debug`, target `test_amr_patch_tree_split_3d`, passed.
- New `unit_amr_patch_tree_sdc_3d`: **1/1 passed**, 9.19 s test / 9.23 s total.
  A three-level branched tree runs two SDC iterations per step, chemistry,
  all molecular transport terms and external force at dt=1e-9 s.
  Uniform reacting node samples match the single-grid SDC result within
  1e-11 scaled state error and all temperatures within 1e-6 K. A subsequent
  nonuniform fine-cell perturbation exercises transport, hydro and reflux.
  Both steps close Euler and element inventories after boundary/source
  accounting within 2e-10 relative bounds and leave covered parents exactly
  synchronized with restriction. Invalid hydro viscosity after source
  evaluation preserves the complete tree and all diagnostic outputs.
- Provenance, Python syntax, project contract and whitespace passed. No
  failed/skipped tests; existing compiler warnings remain. Only the new
  focused SDC test was run; the existing Strang test was not rerun.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_amr_patch_tree_split_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_amr_patch_tree_sdc_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/sdc-tree-loop-checks/` contains build, provenance, CTest and
LastTest logs. AMR application input/diagnostics and reaction-source restart,
SDC subcycling, MPI/EB and external solver/physical validation remain open.
The focused periodic test does not establish all physical-boundary combinations.
