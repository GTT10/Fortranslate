# Independent EB tree inlet mixtures — 2026-09-26

## Change and source

Implementation e37bae9, parent d52f3d8, local main, version 0.246.0.
The shared EB tree application accepts six independent named inlet mixtures,
normalizes their mole fractions and converts them using the active mechanism.
A zero count inherits the initial gas independently on each face. Invalid names,
fractions, unused entries and overrides on non-inflow faces are rejected.
Density and pressure determine inlet temperature through the existing EOS.

The optional PELEF_EB_TREE_INLET_COMPOSITION_V1 envelope stores species count,
six input counts and resolved mass fractions. Changed or removed overrides
reject continuation. Inputs without overrides retain the existing envelope.
The implementation is shared with MPI spray; spray runtime was not rerun here.

Pinned PeleC: bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8.
Pinned PelePhysics: 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.
Mechanism SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3, Ninja Debug.
Existing fixed tests-enabled and selected tests-disabled builds were kept
separate. Selected mechanism: H2O2 full. No backend configuration changed.

Each checker passed **5 successful executions and 10 expected rejections**,
with no failed attempts or skips. Serial, MPI two-rank and continuation on
three ranks have byte-identical final field CSVs across scheduled regridding.
Changing the inlet mixture produces a different final field. Maximum normalized
budget residual: fixed 2.452146e-15; selected 1.102732e-14 (limit 5e-10).
Rejections cover changed hydro method, wrong solver, unknown method, changed
inlet velocity, changed face labels, conflicting Mach input, changed/removed
composition on restart, unknown species and negative mole fractions.
Existing compiler warnings remain; no tolerances or golden data changed.

Commands (all completed successfully):

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_tree_3d_selected pelef_mpi_reactive_eb_tree_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tools/check_eb_tree_ctu_application.py --serial build/spray-amr-reactive-debug/pelef_reactive_eb_tree_3d --mpi build/spray-amr-reactive-debug/pelef_mpi_reactive_eb_tree_3d --mpiexec /usr/bin/mpiexec --input cases/reactive_eb_tree_3d/scheduled.nml --work build/eb-tree-inlet-composition-checks/fixed --inflow --inlet-composition
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tools/check_eb_tree_ctu_application.py --serial build/spray-amr-reactive-selected-debug/pelef_reactive_eb_tree_3d_selected --mpi build/spray-amr-reactive-selected-debug/pelef_mpi_reactive_eb_tree_3d_selected --mpiexec /usr/bin/mpiexec --input cases/reactive_eb_tree_3d/scheduled.nml --work build/eb-tree-inlet-composition-checks/selected --inflow --inlet-composition --selected --uniform-density
python3 -m py_compile tools/check_eb_tree_ctu_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: build/eb-tree-inlet-composition-checks/, including build.log,
build-selected.log, fixed.log, selected.log, fixed/run-sf45i20l and
selected/run-vg_9bfmo. Run folders contain input, stdout, fields and budgets.

## Limits

This is an application input and continuation check, using one overridden
inlet face. Multiple simultaneous mixtures, selected bundles other than H2O2,
MPI spray runtime, full configuration matrix, installed binaries and physical
validation were not checked. Independent single-level inlet mixture input
and domain-wall tree boundaries remain open. Full PeleC parity is unproven.
