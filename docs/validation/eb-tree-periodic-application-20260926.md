# Periodic EB tree gas application — 2026-09-26

Source `a276a7a`, parent `fd155ff`, local main. Fixed/selected serial and MPI
gas executables accept paired periodic faces. Explicit labels reach geometry
construction so extended PLM/CTU geometry wraps donor nodes. Existing boundary
checkpoint identity binds periodic labels without a format change. Unpaired
labels and incompatible opposite EB metrics reject. MPI EB spray input still
rejects periodic faces because its particle coupling is not connected.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Python 3.12.3, Ninja Debug. Existing fixed
tests-enabled and selected tests-disabled builds; H2O2 full selected bundle.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_tree_3d_selected pelef_mpi_reactive_eb_tree_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tools/check_eb_tree_ctu_application.py --serial build/spray-amr-reactive-debug/pelef_reactive_eb_tree_3d --mpi build/spray-amr-reactive-debug/pelef_mpi_reactive_eb_tree_3d --mpiexec /usr/bin/mpiexec --input cases/reactive_eb_tree_3d/scheduled.nml --work build/eb-tree-periodic-application-checks/fixed --periodic
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none python3 tools/check_eb_tree_ctu_application.py --serial build/spray-amr-reactive-selected-debug/pelef_reactive_eb_tree_3d_selected --mpi build/spray-amr-reactive-selected-debug/pelef_mpi_reactive_eb_tree_3d_selected --mpiexec /usr/bin/mpiexec --input cases/reactive_eb_tree_3d/scheduled.nml --work build/eb-tree-periodic-application-checks/selected --periodic --selected --uniform-density
python3 -m py_compile tools/check_eb_tree_ctu_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Both checker runs passed **4 successful executions and 6 expected rejections**,
no failures or skips. Successful runs are serial, MPI two-rank, checkpoint stop
on two ranks, and continuation on three ranks. Scheduled coarsening/refinement
is active; final field CSVs agree byte-for-byte. Maximum budget residuals were
1.521338e-15 fixed and 4.429765e-15 selected. The selected run uses implicit
chemistry and uniform initial density, as requested by its checker flags.
Rejections cover changed hydro method on restart, incompatible solver, unknown
method, changed boundary labels on restart, unpaired periodic faces and
incompatible z-periodic plane geometry.

Provenance, Python syntax, project contract and whitespace checks passed.
Artifacts under `build/eb-tree-periodic-application-checks/`: `build-fixed.log`,
`build-selected.log`, `fixed/run-hm3nu3g6/`, `selected/run-5rszqd34/`.
Mechanism SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

MPI spray periodic coupling remains pending. These bounded gas cases do not
establish arbitrary periodic EB geometry, full configuration coverage,
installed-binary qualification, external solver parity or physical validation.
