# Sparse MPI coarse PPM hydro — 2026-09-25

Source: base `6219726`, implementation `9dd4b09`. GNU Fortran 13.3.0,
Open MPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug, Python 3.12.3.
Reused fixed tests-enabled `build/spray-amr-reactive-debug` and selected H2/O2
tests-disabled `build/spray-amr-reactive-selected-debug` separately.
Artifacts: `build/ppm-mpi-two-level-coarse-checks/`.

Sparse two-level coarse halo exchange accepts optional width two (default) or
four, with collective width agreement. Layer-by-layer exchange forwards cells
across slabs narrower than the halo. Physical boundary filling uses a width-sized
edge stencil and preserves candidate publication on success. No global field
replication was added. The supplied-halo CTU kernel now accepts one owned cell,
which has all reconstruction neighbors in its supplied halo.

`advance_mpi_amr_sparse_boundary_hydro_3d` supports a single conservative PPM
CTU step with collective flux/EOS acceptance and the existing surface ledger.
The periodic API forwards PPM through that path and preserves its upper-only
transverse flux interface. Default PCM/PLM branches are unchanged. Fine PPM
interpolation, subcycling, hierarchy reflux and application selection still need
connection; this increment does not enable a full two-level MPI PPM run.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_reactive_3d --parallel 2`
  and `cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_mpi_amr_reactive_3d_runtime --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^mpi_amr_ppm_coarse_3d_np8$' --output-on-failure`:
  **1/1 passed**, 4.51 s test / 4.55 s total; no failures or skips.
  Eight ranks own one x cell each of an 8x8x8 coarse grid. Four-layer physical
  and pressure halos match serial faces/corners exactly. Tests include moving
  viscous walls, inflow, outflow and periodic transverse directions, plus width,
  face and rank-local metadata rejection with preserved halo outputs.
  PPM periodic/physical/pressure coarse state, temperature and fluxes match the
  serial reference exactly. Physical surface ledger reduction uses its existing
  tolerance. CFL and rank-dependent reconstruction rejection preserve state and
  zero diagnostics. Periodic upper-only transverse flux dispatch also matches.
- In the same environment:
  `/usr/bin/mpiexec -n 8 build/spray-amr-reactive-debug/test_mpi_amr_reactive_3d --halos-only`:
  passed, retaining default two-layer behavior. The log contains nonfatal display
  authorization messages and floating-point underflow/denormal flags; exit was 0.
- `python3 -m py_compile tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

Selected runtime execution, full matrix, complete two-level MPI PPM, convergence,
reference-solver parity and physical validation were not run. Compiler warnings
remain. Reference primitive-pressure tracing fidelity remains open.
