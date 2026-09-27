# MPI embedded-boundary MOL rate — 2026-09-27

`compute_mol_mpi_eb_tree_hydro_3d` now exchanges three-layer state/temperature
halos using prepared extended metrics and evaluates the MOL patch operator on
owned cells. It returns owner-local spatial rates and dt-integrated face
fluxes. Rank agreement covers interval, reconstruction order and selected
depth; geometry interiors must match the physical patches. Results are
published only after collective success. Optional ancestor state/temperature
interpolation is forwarded for future level stepping, and requires an explicit
target depth. This API performs no state update, redistribution, averaging or
reflux.

## Completed checks

GNU Fortran 13.3.0 / OpenMPI 4.1.6, existing Debug MPI/spray core build,
tests/SUNDIALS disabled. Final commands:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_support --parallel 2 > build/mpi-eb-mol-build-final.log 2>&1
env -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/mpifort -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -I build/heptane-simple-mpi-debug/mpi_modules -I build/heptane-simple-mpi-debug/modules tests/unit/test_mpi_eb_tree_hydro_3d.F90 build/heptane-simple-mpi-debug/libpelef_mpi_support.a build/heptane-simple-mpi-debug/libpelef_core.a -o build/mpi-eb-mol-unit/check-final
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none /usr/bin/mpiexec -n 2 build/mpi-eb-mol-unit/check-final ctu-wall > build/mpi-eb-mol-unit/run-final.log 2>&1
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

The focused two-rank executable extends the existing nested five-node wall
fixture. It compares MOL cell rates and all three flux arrays against serial
halo preparation plus the patch operator, using a fixed `1e-12` relative gate.
Distributed ownership and all-on-rank-zero layouts both pass, including the
rank with no cells. A rank-local invalid reconstruction order rejects
collectively and preserves previous rates and flux arrays. Existing temporal
source and CTU checks also pass; CTU MPI inventory residual remains
`4.6333897726434849e-15` in both layouts.

An initial build/run passed before adding the explicit target-depth requirement
for ancestor interpolation. The final build and focused executable were rerun
after that guard. No full matrix, selected-runtime rebuild, recursive MOL
time integration or C++ spatial-operator comparison was performed. Linker
executable-stack and MPI X11 messages remain nonfatal. No Python source changed.

MOL old/new time integration, stage transport/reaction/particle coupling,
redistribution/reflux routing and application controls remain open. The separate
jet continuation published step 2300 at `0.0004084 s` and remained live; full
jet completion is still unverified.
