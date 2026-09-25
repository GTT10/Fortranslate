# MPI EB tree adiabatic wall hydro — 2026-09-26

Implementation feb7b87, parent 7a50fd6, local main, version 0.246.0.
MPI flux synchronization now admits slip_wall and no_slip_wall and assigns
distinct collectively checked face codes. Existing common/recursive hydro
already forwards these conditions through reflected halos and physical CTU
fluxes. Accepting the labels completes peer and coarse/fine flux reconciliation
for these wall hydro paths. Periodic seam logic is unchanged.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Python 3.12.3, Ninja Debug. Existing
fixed tests-enabled build/spray-amr-reactive-debug, no backend changes.
Three-rank common/recursive CTU wall tests passed **2/2**, 6.66 and 18.78 s
(CTest total 25.49 s), no failures or skips.

The z-planar EB fixture has five patches on three levels, x-slip walls,
y-no-slip walls and z-outflow. Distributed ownership and all-on-rank-zero
ownership exercise empty ranks. State, temperature and boundary ledger match
the corresponding serial results exactly. Embedded-wall ledger and composite
inventory satisfy existing 3e-14 bounds. Maximum MPI residuals are
4.6333897726434849e-15 common and 4.4679115664776463e-15 recursive.
The serial reference is rerun in these tests, including recursive CTU walls.
Rank-disagreeing labels and the existing invalid state, timestep, geometry,
coefficient and solver cases reject collectively and preserve public outputs.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_hydro_3d test_mpi_eb_tree_subcycle_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_wall_(hydro|subcycle_hydro)_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, project contract and whitespace passed. Compiler warnings
remain. No tolerance or golden-data changes.
Artifacts: build/mpi-eb-wall-hydro-checks/build.log and LastTest.log.
Mechanism SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
Pinned PeleC: bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8.
Pinned PelePhysics: 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.

## Limits

Wall geometry is tangential to the planar EB. General reflected geometry,
viscous/thermal wall transport, full reacting split, application input and
isothermal tree walls remain separate work. No full matrix, installed audit,
external-solver parity or physical wall-flow validation ran. Complete PeleC
parity remains unproven.
