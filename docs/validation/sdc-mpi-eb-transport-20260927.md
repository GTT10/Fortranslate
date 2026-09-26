# MPI EB tree SDC transport stage — 2026-09-27

## Source and behavior

Implementation `feb5e331`, based on `92bb9c23` on local `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`compute_sdc_mpi_eb_tree_transport_3d` evaluates molecular transport on owned
patches and returns rates per fluid volume, dt-integrated face flux densities
and a global minimum limiter. It leaves state and temperature unchanged and
performs no restriction, redistribution or reflux. Outputs are private until
all ranks report success. Timestep, transport flags and boundary controls
must agree across ranks.

The MPI EB gas halo adds optional `restrict_snapshot` (default true). SDC
passes false to retain provisional parent states. This path still validates
replicated geometry, topology, owners and state/temperature shapes before
exchanging donors. The option itself is agreed collectively before branching.

## Completed focused checks

GNU Fortran 13.3.0, configured `/usr/bin/mpifort` and `/usr/bin/mpiexec`
(OpenMPI/OpenRTE 4.1.6), Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance preflight passed for `mechanisms/h2o2_full.json`.
- `test_mpi_eb_tree_transport_3d` built successfully.
- CTest `unit_mpi_eb_tree_transport_3d_2`: 1/1 passed, 8.84 s test / 8.88 s total.
- New checks use three levels, five patches and two ranks, with round-robin
  ownership and all patches on rank 0 (rank 1 idle). A parent state deliberately
  differs from the restricted fine state. SDC rates, all three face flux arrays
  and minimum limiter match serial exactly. Input state/temperature are exact,
  and non-owned outputs remain unallocated.
- Invalid last-patch density and, independently, a negative timestep on rank 1
  reject globally and retain prior rate, flux and limiter outputs.
- Rank-dependent halo restriction is rejected before publishing halos.
- Existing split transport and timestep checks run in the same executable.
- Python syntax, project contract and `git diff --check` passed.

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_transport_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_transport_3d_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Logs: `build/sdc-mpi-eb-transport-checks/`.
No isothermal test mode, selected build, full configuration matrix, installed
binary audit, coupled MPI EB SDC application, scaling or physical validation
was run. Selected-level transport and coupled MPI EB SDC hydro remain open.
