# Two-level AMR PPM inflow replacement — 2026-09-25

Source: base `f404ba3`, face-composition correction `3a37061`, AMR wiring
`c75bfe9`. GNU Fortran 13.3.0, OpenMPI 4.1.6, Ninja Debug, Python 3.12.3;
reused fixed tests-enabled `build/spray-amr-reactive-debug`.
Artifacts: `build/ppm-two-level-inflow-checks/`.

Serial fine-patch CTU now passes its existing physical-touch mask intersected
with inflow faces. MPI coarse CTU enables X inflow replacement only on the
rank owning the corresponding global domain end. Sparse fine CTU additionally
requires the fine patch to touch that physical boundary and the rank to own its
first/last fine cell. Y/Z masks use patch physical extent. Empty fine ranks
retain their existing empty-state handling. Periodic and inter-rank faces keep
CTU fluxes. Patch-tree and regular MPI spray wiring remain separate unfinished work.

## Failure and correction

The first two-test run passed coarse MPI comparison but failed the serial
reference of the reacting/transporting hierarchy. A one-rank diagnostic located
the failure at lower-Y inflow Riemann evaluation: independently limited species
slopes produced a face composition sum `0.999999999989751`, violating the
existing EOS composition contract. The initial and diagnostic logs are retained.

The PLM predictor now verifies nonnegative finite face fractions and normalizes
each face composition before EOS use. Negative/invalid predictions still fail
transactionally. No EOS tolerance was relaxed. Cell inventories remain updated
through the existing conservative flux and species-closure path. This face
normalization is a PeleF admissibility policy in addition to the pinned PLM
arithmetic; it must not be described as bitwise reference PLM parity.

## Completed checks

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`: passed.
- With `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset:
  `cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_reactive_3d test_reactive_ppm_halo_3d --parallel 2`: passed.
- With `DISPLAY`, `LD_LIBRARY_PATH`, `LIBRARY_PATH` unset,
  `HWLOC_COMPONENTS=-x86`, `OMPI_MCA_hwloc_base_binding_policy=none`:
  `ctest --test-dir build/spray-amr-reactive-debug -R '^(mpi_amr_ppm_(coarse|hierarchy)_3d_np8|unit_reactive_ppm_halo_3d)$' --output-on-failure`:
  **3/3 passed**, 27.78 s total; no final failures/skips.
  - Halo unit, 5.76 s: added independently limited species normalization;
    existing PLM, six-face replacement, inventory and rollback checks remain.
  - Eight-rank hierarchy, 19.75 s: exact serial state/temperature/reflux parity
    across periodic, interior, physical-edge and full-refinement geometries;
    lower-X/lower-Y inflow, pressure outlet, reacting transport split, empty
    fine ranks and rollback. Existing boundary-ledger tolerances retained.
  - Eight-rank coarse, 2.23 s: one coarse cell per rank, four-layer halos,
    both X inflow domain ends plus lower-Y inflow, exact serial flux/state
    agreement, pressure/periodic cases and rollback.
- `python3 -m py_compile tools/check_project_contract.py`,
  `python3 tools/check_project_contract.py --source .`, `git diff --check`: passed.

Selected build/execution, application restart, full matrix, external CTU oracle,
convergence and physical validation were not run. Compiler warnings and
nonfatal display-authorization messages remain in logs. The initial failed run
also reported floating-point status flags; it is not counted as a pass.
