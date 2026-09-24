# Sparse two-level MPI PPM hierarchy — 2026-09-25

Source: base `219f9f3`, implementation `172f29a`. GNU Fortran 13.3.0,
Open MPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug, Python 3.12.3. Reused fixed
MPI tests-enabled `build/spray-amr-reactive-debug` and selected H2/O2 MPI
tests-disabled `build/spray-amr-reactive-selected-debug` separately.
Artifacts: `build/ppm-mpi-two-level-hierarchy-checks/`.

Fine PPM steps now exchange four fine layers, interpolate coarse donors from
four-layer coarse halos at each substep start, trace the CTU fluxes and perform
one conservative update. Interpolation detects the supplied coarse halo width
for periodic and physical slabs. Empty fine ranks retain their existing
collective participation. Existing interface flux accumulation, reflux,
restriction, EOS recovery and hierarchy candidate publication are retained.
Reaction and transport splits accept PPM. Coarse input temperatures follow the
serial PPM path without redundant pre-recovery.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_reactive_3d --parallel 2`
  and `cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_mpi_amr_reactive_3d_runtime --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^mpi_amr_ppm_hierarchy_3d_np8$' --output-on-failure`:
  **1/1 passed**, 34.02 s test / 34.06 s total; no failures or skips.
- Eight ranks own one coarse x cell each. The ratio-two tests cover a periodic
  interior patch, physical interior/lower-corner/upper-corner patches and complete
  refinement, including ranks with no fine cells. Coarse/fine conserved state,
  temperature and maximum reflux correction equal the serial result exactly.
  Physical surface inventories retain the existing tolerance. The upper corner
  uses a prescribed x outlet pressure. A lower-corner full reaction/transport
  split exercises moving isothermal walls, elemental conservation and chemistry
  activity. Invalid rank-local boundary metadata and a periodic CFL failure
  preserve the hierarchy and zero the failure diagnostics.
- `python3 -m py_compile tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

Selected execution, other refinement ratios, full matrix, convergence, external
solver parity and physical validation were not run. Compiler warnings remain.
Application/regrid/restart connection is recorded separately. Reference
primitive-pressure tracing fidelity remains open.
