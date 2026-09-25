# EB tree adiabatic wall application — 2026-09-26

Implementation 50f27ab, parent 6f3f0bc, local main, version 0.246.0.
Fixed/selected serial and MPI gas applications accept slip_wall and no_slip_wall
face labels. Existing configured geometry supplies reflected CTU extensions;
existing boundary checkpoint context binds the labels. MPI tree spray explicitly
rejects domain-wall labels until particle routing is connected. Isothermal
walls remain unsupported in this entry point.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3, Ninja Debug. Existing fixed
tests-enabled and selected tests-disabled H2O2 full builds remain separate.
No backend configuration changes.

Each checker passed **4 successful executions and 6 expected rejections**,
no failures or skips. Serial, MPI two-rank and restart on three ranks produce
byte-identical final field CSVs across scheduled coarsening/refinement. The
fixture encloses the domain with x-slip and y/z-no-slip walls; this includes
walls normal to the planar EB and exercises reflected configured geometry.
Maximum budget residuals: fixed 1.521338e-15; selected 4.429765e-15.
A supplementary read of the six serial/MPI/resumed final budget CSVs verifies
exactly zero external mass, species and energy exchange (components 1 and 5:).
Rejections cover changed hydro method, wrong solver, unknown method, changed
wall type on restart, unsupported isothermal input and conflicting Mach mode.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_tree_3d_selected pelef_mpi_reactive_eb_tree_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tools/check_eb_tree_ctu_application.py --serial build/spray-amr-reactive-debug/pelef_reactive_eb_tree_3d --mpi build/spray-amr-reactive-debug/pelef_mpi_reactive_eb_tree_3d --mpiexec /usr/bin/mpiexec --input cases/reactive_eb_tree_3d/scheduled.nml --work build/eb-tree-wall-application-checks/fixed --wall
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tools/check_eb_tree_ctu_application.py --serial build/spray-amr-reactive-selected-debug/pelef_reactive_eb_tree_3d_selected --mpi build/spray-amr-reactive-selected-debug/pelef_mpi_reactive_eb_tree_3d_selected --mpiexec /usr/bin/mpiexec --input cases/reactive_eb_tree_3d/scheduled.nml --work build/eb-tree-wall-application-checks/selected --wall --selected --uniform-density
python3 -m py_compile tools/check_eb_tree_ctu_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, project contract and whitespace passed. Existing compiler
warnings remain; no tolerance or golden-data changes.
Artifacts: build/eb-tree-wall-application-checks/{build.log,build-selected.log,
fixed.log,selected.log,wall_flux_check.log}, fixed/run-ne1m35l0 and
selected/run-9iuva4di. Each run retains input, stdout, fields and budget CSVs.
Mechanism SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
Pinned PeleC: bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8.
Pinned PelePhysics: 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.

## Limits

This qualifies the planar application fixture and changed-rank continuation.
General level-set wall intersections, other selected mechanisms, isothermal
walls and MPI spray wall routing remain separate work. No full matrix,
installed audit, external solver parity or physical validation ran.
Complete PeleC parity remains unproven.
