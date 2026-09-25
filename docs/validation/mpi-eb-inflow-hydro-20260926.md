# MPI EB prescribed inflow hydro — 2026-09-26

## Implementation

Source 4a49037, parent 631b12b, local main. Common and recursive MPI EB hydro
APIs now forward boundary_faces/inflow_state/inflow_temperature into their
owner-local halo exchanges. CTU also receives the same explicit face labels
for physical boundary stencils. Existing halo collective validation rejects
inconsistent labels/states before publication. PCM forwarding is present but
the new runtime evidence below uses CTU.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6 via /usr/bin/mpiexec, Ninja Debug and Python
3.12.3; existing tests-enabled build/spray-amr-reactive-debug. No backend changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_hydro_3d test_mpi_eb_tree_subcycle_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_ctu_inflow_(subcycle_)?hydro_3d_3$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed, 0 failed, 0 skipped**, 26.50 s total: common 6.86 s, recursive
19.59 s. Five patches span three levels with a z-normal embedded plane and
prescribed x-low/y-high inlets, using CTU and nonzero artificial viscosity.
Both three-rank fixtures test distributed ownership and all patches on the
last rank, leaving idle ranks. Gas state and temperature match the serial
update exactly. Exterior ledgers match exactly; embedded-wall and total
inventory residuals satisfy the existing 3e-14 tolerance.

Existing collective failure/rollback cases for invalid solver/state, differing
dt/difmag, halo geometry and extended geometry remain in the fixtures. A new
case changes inlet density on only rank zero and verifies collective rejection
without publishing gas, temperature or ledgers. Input density is restored by
exact assignment between ownership cases.

Provenance, syntax, project contract and whitespace checks passed. No failed
build/test attempts; no tolerances or golden data changed. Existing compiler
warnings remain. Mechanism SHA:
0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
Artifacts: build/mpi-eb-inflow-hydro-checks/build.log and LastTest.log.

## Remaining scope

MPI transport/full-split forwarding, timestep controls, application input and
checkpoint context remain to connect. Periodic/wall tree dispatch is not added.
No full matrix, installed audit, external solver parity or physical validation
was run; complete PeleC parity remains unproven.
