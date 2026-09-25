# SDC reaction stage over an AMR patch tree — 2026-09-26

Implementation `6fb2a1a9`, parent `21472792`, local main, version 0.246.0.

`react_sdc_amr_patch_tree_3d` applies the forced SDC reaction stage to each
patch using interval-start state/temperature and caller-supplied nonreacting
rates. State, temperature and reaction-source trees are published only after
every patch succeeds. Input tree layouts, component counts and scalar controls
are checked. This stage evolves each level's valid cells, including covered
parents; restriction belongs to the hierarchy's later synchronization phase.
Unlike the existing split chemistry routine, it does not average down here.

## Completed checks

- GNU Fortran 13.3.0 Debug fixed tests-enabled build
  `build/spray-amr-reactive-debug`, target `test_amr_patch_tree_chemistry_3d`, passed.
- `unit_amr_patch_tree_chemistry_3d`: **1/1 passed**, 63.58 s test / 63.63 s total.
  Added a three-level branched tree with uniform reacting state and nonzero
  external rates. Node sample states/reaction sources and all temperatures
  match the same standalone reaction kernel on a uniform grid (bounds 1e-10,
  1e-8 and 1e-7 K respectively). A negative-density cell in the last patch
  fails after prior local work without publishing any state, temperature or
  reaction-source tree. Existing split chemistry/Strang checks remain active.
- Provenance, Python syntax, project contract and whitespace passed. No
  failed or skipped tests; existing compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_amr_patch_tree_chemistry_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_amr_patch_tree_chemistry_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/sdc-tree-reaction-checks/` contains build, provenance, CTest
and LastTest logs. This connects only the reaction stage. Hierarchy-wide
hydro/source iteration, coarse/fine flux correction, final restriction and
MPI dispatch remain open; no complete AMR SDC or external solver parity claim.
