# MPI EB SDC selected-level stages — 2026-09-27

Implementation `8496ccb0`, based on `0993dce8` on local `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

MPI SDC hydro and transport stages now accept optional zero-based
`target_depth`, ancestor state/temperature endpoint fields and gas fractions.
Hydro additionally forwards source endpoint fields and independent fractions.
Level choice is agreed collectively. Ancestor inputs require a selected level;
complete groups, shapes and temporal controls are checked by the halo APIs.
Fraction arrays retain the MPI halo convention: root-first, length equal to
the maximum tree level.

Only owned patches on the selected level execute the numerical kernel. Rates
are zero on other owned levels, and inactive flux arrays remain unallocated.
State/temperature remain input-only, and output publication remains collective.
Halo preparation still spans the tree, so valid owner-local donor fields and
geometry must be provided for all levels. This adds no restriction, reflux,
redistribution or recursive time-stepping driver.

## Completed focused checks

GNU Fortran 13.3.0, `/usr/bin/mpifort` and `/usr/bin/mpiexec`
(OpenMPI/OpenRTE 4.1.6), Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance and both touched test targets built successfully.
- CTest: 2/2 passed, 25.57 s total. Two-rank transport: 8.51 s;
  three-rank CTU hydro: 17.02 s.
- New comparisons use three levels, five patches, selected depth 2, round-robin
  ownership and all work on rank 0 (idle ranks included). Gas end states are
  1.01 times their start state; gas fraction is 0.5.
- Hydro predictor source endpoints differ independently and use fraction 0.25.
  Selected rates and all three face-flux arrays match serial exactly.
  Transport's minimum limiter is exact. Inactive rates match serial zeros and
  inactive flux arrays remain unallocated.
- Existing common-stage, temporal-source, rollback and split-step checks remain
  in the same executables and passed. New level-specific failure cases were
  inspected in code but not separately fault-injected.
- Python syntax, project contract and `git diff --check` passed.

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_hydro_3d test_mpi_eb_tree_transport_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_(ctu_hydro_3d_3|transport_3d_2)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/sdc-mpi-eb-level-stages-checks/`.
No selected-mechanism build, full matrix, installed-binary audit, selected-level
coupled SDC driver, recursive MPI EB application or physical validation was run.
