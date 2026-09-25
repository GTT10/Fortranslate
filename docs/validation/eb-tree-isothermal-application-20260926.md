# EB tree isothermal application integration — 2026-09-26

Implementation `48041e7`, parent `7fe5cda`, local main, version 0.246.0.
The shared fixed/selected serial/MPI gas entry point now accepts
`boundary_faces='isothermal_wall'` and `wall_temperature(6)` in kelvin.
It forwards temperatures to common/recursive full splits and automatic timestep
selection. Isothermal input requires enabled molecular transport and six finite
positive temperatures. Nonthermal inputs retain zero defaults and no new header.
The optional `PELEF_EB_TREE_ISOTHERMAL_V1` checkpoint prefix binds all six values;
changed, removed or missing thermal context rejects before continuation.
MPI EB spray thermal input still rejects explicitly until coupled forwarding exists.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3, Ninja Debug.
Existing fixed tests-enabled and selected tests-disabled build trees were used;
no compiler/MPI/backend changes. Selected mechanism is full H2/O2.

Both fixed and selected application checkers passed **4 successful runs and
8 expected rejections**, with no unexpected failures. Each runs serial, MPI2,
MPI2 checkpoint stop and MPI3 continuation. The scheduled planar fixture uses
CTU recursive stepping, automatic timestep selection, six 2000 K walls and
collapse/refinement. Fixed chemistry is explicit; selected chemistry is implicit.
Final cell CSVs are byte-identical across serial, MPI2 and restarted MPI3.
Boundary energy is positive and boundary mass/species transfers are exactly zero.
Composite residual maxima: fixed 2.137836e-15, selected 2.697383e-15.

Rejections cover method/solver incompatibilities, changed restart temperature,
missing/negative temperature, disabled transport, and removed thermal boundaries.
A fixed-mechanism CTest entry was registered; the checker was run directly for
both configurations, not separately counted as a CTest execution.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_tree_3d_selected pelef_mpi_reactive_eb_tree_3d_selected --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tools/check_eb_tree_ctu_application.py --serial build/spray-amr-reactive-debug/pelef_reactive_eb_tree_3d --mpi build/spray-amr-reactive-debug/pelef_mpi_reactive_eb_tree_3d --mpiexec /usr/bin/mpiexec --input cases/reactive_eb_tree_3d/scheduled.nml --work build/eb-tree-isothermal-application-checks/fixed --isothermal
env -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tools/check_eb_tree_ctu_application.py --serial build/spray-amr-reactive-selected-debug/pelef_reactive_eb_tree_3d_selected --mpi build/spray-amr-reactive-selected-debug/pelef_mpi_reactive_eb_tree_3d_selected --mpiexec /usr/bin/mpiexec --input cases/reactive_eb_tree_3d/scheduled.nml --work build/eb-tree-isothermal-application-checks/selected --isothermal --selected
python3 -m py_compile tools/check_eb_tree_ctu_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, Python syntax, project contract and whitespace checks passed.
No tolerances or golden values changed. Build warnings remain.
Artifacts under `build/eb-tree-isothermal-application-checks/`:
`build-fixed.log`, `build-selected.log`, `fixed.log`, `selected.log`,
`fixed/run-z8euye44`, `selected/run-_k9ut0jj`.
Mechanism SHA: `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
PeleC: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
PelePhysics: `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.

## Limits

MPI EB spray thermal forwarding is still open. No full matrix, clean install,
external solver comparison or physical wall-flow validation ran. This application
integration is not evidence of complete PeleC parity.
