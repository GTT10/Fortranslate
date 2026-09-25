# Serial/MPI EB tree CTU applications — 2026-09-25

Implementation `91352d6`, test correction `92b0ee3`, base `70569bd`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

The shared fixed/selected serial/MPI gas application accepts
`hydro_method='plm_ctu'` with `riemann_solver='pelec'`. Common and recursive
full-step calls receive an extended geometry cache. Padded finest geometry is
built once and coarsened for the refinement lattice; patch caches rebuild when
topology descriptors change, including after restore. The existing checkpoint
solver identity binds the method. PCM defaults and checkpoint identity remain
unchanged. EB spray rejects this selector explicitly.

Initial regular/cut density and pressure are now namelist inputs, retaining
0.31/0.93 kg/m3 and 135000 Pa defaults. Finite positive validation and existing
checkpoint identity cover these inputs.

## Evidence

GNU Fortran 13.3.0, OpenMPI `/usr/bin/mpiexec`, Python 3.12.3, Ninja Debug.
Existing builds `build/spray-amr-reactive-debug` and
`build/spray-amr-reactive-selected-debug` (H2O2 bundle). No toolchain change.

The scheduled three-level tree coarsens to one node and refines to four nodes.
Checks use recursive CTU with reaction and transport, serial full execution,
MPI-2 full execution, MPI-2 checkpoint, and MPI-3 continuation. Final CSVs match
serial byte-for-byte. Integrated CSV fields agree with inventory ledgers;
conservation retains the existing `5e-10` criterion and inventory summation
`3e-14` criterion. Method mismatch, incompatible solver and unknown method each
reject as expected.

- Fixed: CTest `regression_eb_tree_ctu_application` **1/1 PASS**, 59.77 s
  (59.81 s total); four successful processes and three expected rejections.
  Maximum residual `1.521338e-15`.
- Selected H2O2, implicit chemistry and default density contrast: all four
  processes completed; final CSV equality and post-run ledger checks passed,
  followed by three expected rejection runs. Maximum residual `1.660887e-15`.
  Timestep retries produced nine steps to 3 ns. Step-2 checkpoint is at 0.5 ns,
  before coarsening, so continuation crosses both scheduled topology changes.
- Supplemental selected implicit case with `cut_density=.31`: four successes
  and three expected rejections; maximum residual `4.429765e-15`.

Initial selected input lacked composition and rejected before execution; fixed
H2/O2/N2 mole fractions were then supplied. The selected explicit trial timed
out at 180 s and is **not** counted as passed. The default-density implicit trial
completed, but the original checker incorrectly required coarsening by step 2.
That assumption was removed for selected cases. The already generated successful
outputs were checked without rerunning simulations by
`build/eb-tree-ctu-app-checks/check_existing_selected.py`; it reran only the three
rejection cases. The final checked-in selected mode retains the density contrast;
`--uniform-density` reproduces the supplemental case. No EOS limits or numerical
tolerances were changed. A GDB snapshot during the long run showed active
transport computation at a shortened interval, not a deadlock.

Commands (from repository root):

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_tree_3d_selected pelef_mpi_reactive_eb_tree_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^regression_eb_tree_ctu_application$' --output-on-failure
# Direct selected checker invocation, with the same environment:
python3 tools/check_eb_tree_ctu_application.py --serial build/spray-amr-reactive-selected-debug/pelef_reactive_eb_tree_3d_selected --mpi build/spray-amr-reactive-selected-debug/pelef_mpi_reactive_eb_tree_3d_selected --mpiexec /usr/bin/mpiexec --input cases/reactive_eb_tree_3d/scheduled.nml --work build/eb-tree-ctu-app-checks/selected --selected
python3 -m py_compile tools/check_eb_tree_ctu_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, contract and whitespace checks passed. Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Build logs and `LastTest-final.log` are in `build/eb-tree-ctu-app-checks/`.
Selected artifacts: `selected/run-w7cww15r` (contrast), `selected/run-e377zir0`
(uniform), `selected/run-wdz6vlow` (explicit timeout), `selected/run-xt6rmd9w`
(missing composition); `selected-stack.log` records the diagnostic snapshot.
Earlier fixed direct check: `fixed/run-18obacn4`; final fixed artifacts are named
in `LastTest-final.log`. Compiler warnings remain. Common-step library checks
are separate evidence; this application check used subcycling. Adaptive tagging,
other EB shapes, spray coupling, full external parity and physical validation
were not newly qualified here.
