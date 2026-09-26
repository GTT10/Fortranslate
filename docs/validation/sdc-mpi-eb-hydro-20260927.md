# MPI EB tree SDC hydro stage — 2026-09-27

## Source and behavior

Implementation `5b35915e`, based on `2f082ae9` on local `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`compute_sdc_mpi_eb_tree_hydro_3d` connects the serial source-aware CTU kernel
to owner-local MPI state and source fields. It validates replicated physical
and six-layer extended geometry, exchanges gas without restricting provisional
parents, computes flattening, and returns rates and dt-integrated fluxes only
after collective success. State and temperature are input-only. No restriction,
redistribution, reflux or coupled timestep is performed by this stage.

`fill_sdc_mpi_eb_tree_source_halo_3d` exchanges signed source values using
same-level donors followed by ancestor PCM donors. Physical boundaries use
scalar extrapolation or periodic wrapping. There is no EOS conversion or
momentum reflection. Covered targets remain zero; unused covered source NaNs
are ignored. Open targets require a fluid donor. Owner-local target volume
fractions are checked for shape, range and agreement at physical cells.
Width, component count and boundary behavior must agree across ranks.

## Completed focused checks

GNU Fortran 13.3.0, configured `/usr/bin/mpifort` and `/usr/bin/mpiexec`
(OpenMPI/OpenRTE 4.1.6), Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance preflight passed for `mechanisms/h2o2_full.json`.
- `test_mpi_eb_tree_hydro_3d` built successfully.
- CTest `unit_mpi_eb_tree_ctu_hydro_3d_3`: 1/1 passed, 12.21 s test / 12.26 s total.
- New checks use three levels, five patches and three ranks, with round-robin
  ownership and all patches on rank 0 (two idle ranks). A provisional root
  differs from restricted fine state. Predictor sources include positive and
  negative momentum rates, energy rates and covered-cell NaNs.
- All rate and face-flux arrays match serial exactly. Input state and temperature
  remain exact; non-owned output arrays remain unallocated.
- An active source NaN on the last patch's owner rejects globally and retains
  prior rates and fluxes. A separate invalid timestep on the last rank, including
  an idle rank, also rejects globally with unchanged outputs.
- Existing CTU split-step regression checks remain in the same executable.
- Python syntax, project contract and `git diff --check` passed.

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_ctu_hydro_3d_3$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Logs: `build/sdc-mpi-eb-hydro-checks/`.
The executed test uses outflow boundaries. Separate inflow, periodic and wall
modes, selected builds, full matrices, installed-binary audit, scaling and
physical validation were not run. Temporal source interpolation, selected-level
MPI EB hydro and the coupled MPI EB SDC application remain open.
