# EB tree prescribed inflow application — 2026-09-26

## Source and behavior

Implementation 4c2940d, parent 7be73d2, local main. Serial/MPI fixed and selected
EB tree gas applications accept explicit inflow/outflow faces with density,
pressure and velocity per face. Inlet composition comes from initial gas.
Mach and prescribed-state inputs are mutually exclusive. Tree gas updates
and automatic timestep selection receive the generated conserved states.
MPI inlet generation runs on the root-patch owner, agrees on EOS validity,
and broadcasts six states and temperatures to all ranks. Inlet arrays are
created once and retained when restart restores a scheduled topology.

The optional PELEF_EB_TREE_PRESCRIBED_INLET_V1 envelope binds labels and all
30 density/pressure/velocity values. Missing/changed controls reject restart.
Legacy checkpoint headers remain unchanged. MPI spray explicitly rejects this
new input until its wrapper is connected.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6 via /usr/bin/mpiexec, Ninja Debug,
Python 3.12.3. Existing fixed tests-enabled build/spray-amr-reactive-debug and
selected tests-disabled build/spray-amr-reactive-selected-debug, H2O2 full bundle.
No backend changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_tree_3d_selected pelef_mpi_reactive_eb_tree_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tools/check_eb_tree_ctu_application.py --serial build/spray-amr-reactive-debug/pelef_reactive_eb_tree_3d --mpi build/spray-amr-reactive-debug/pelef_mpi_reactive_eb_tree_3d --mpiexec /usr/bin/mpiexec --input cases/reactive_eb_tree_3d/scheduled.nml --work build/eb-tree-inflow-application-checks/fixed --inflow
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tools/check_eb_tree_ctu_application.py --serial build/spray-amr-reactive-selected-debug/pelef_reactive_eb_tree_3d_selected --mpi build/spray-amr-reactive-selected-debug/pelef_mpi_reactive_eb_tree_3d_selected --mpiexec /usr/bin/mpiexec --input cases/reactive_eb_tree_3d/scheduled.nml --work build/eb-tree-inflow-application-checks/selected --inflow --selected --uniform-density
python3 -m py_compile tools/check_eb_tree_ctu_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final fixed and selected checkers each passed **4 successful executions and
6 expected rejections**, no skips. The executions cover serial, MPI two-rank,
checkpoint stop on two ranks and continuation on three ranks. Scheduled
coarsening and refinement are active. All final field CSVs match byte-for-byte.
Worst budget residuals: fixed 2.220685e-15, selected 1.102559e-14. Rejections
cover changed hydro, wrong solver, unknown method, changed inlet velocity,
changed face labels and conflicting Mach input. No tolerances were loosened.

Two earlier attempts for each mechanism exposed implementation defects:
inlet states were initially generated only in the serial branch; after fixing
MPI generation, scheduled restart attempted to allocate them twice. Both were
fixed before the final complete reruns. MPI launch logs also include X11
authorization messages; these did not prevent the final runs. Compiler
warnings remain. Provenance, syntax, project contract and whitespace passed.
Mechanism SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.

Artifacts under build/eb-tree-inflow-application-checks/:

- Final: fixed/run-berglqxo and selected/run-wjra0ft0.
- Initial MPI failures: fixed/run-zudyhfd_ and selected/run-1lu47zqe.
- Initial restart failures: fixed/run-zfcd3kth and selected/run-f2byfyu1.
- build.log and build-selected.log contain final build output.

## Limits

This qualifies subcycled CTU gas application dispatch and changed-rank restart.
Automatic timestep forwarding is compiled; runtime parity of that API was
checked in the preceding library increment, not separately rerun here.
Independent inlet composition, MPI spray forwarding, and periodic/wall tree
boundaries remain open. No full matrix, installed audit, external solver parity
or physical validation was run. Complete PeleC parity remains unproven.
