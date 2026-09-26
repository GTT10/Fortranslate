# MPI EB SDC ancestor-source time interpolation — 2026-09-27

Implementation `f545b663`, based on `3661295c` on local `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`fill_sdc_mpi_eb_tree_source_halo_3d` now accepts optional `ancestor_start`,
`ancestor_end` and `ancestor_fraction`. The group must be complete and present
on every rank. Fractions are finite, in [0,1], root-first and sized to the
maximum tree level. Values must agree across ranks before donor exchange.
Same-level and interior cells still use the current source field. Coarser
ancestor donors use their depth's interpolated endpoint value.

Current/start/end requests are exchanged separately. Endpoints with zero
weight are never fetched, so unused endpoint NaNs do not reject a valid stage.
Signed source values are neither clipped nor converted through an EOS. The
complete output stays private until all exchanges and finite-result checks
succeed collectively. No time-stepping or level-selection driver is added here.

## Completed focused checks

GNU Fortran 13.3.0, `/usr/bin/mpifort` and `/usr/bin/mpiexec`
(OpenMPI/OpenRTE 4.1.6), Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance and `test_mpi_eb_tree_hydro_3d` build passed.
- CTest `unit_mpi_eb_tree_ctu_hydro_3d_3`: 1/1 passed,
  12.37 s test / 12.41 s total.
- New three-level, five-patch source-halo checks use three ranks with both
  round-robin owners and all work on rank 0 (two idle ranks).
- Full width-3 source halos match serial exactly at fractions 0, 0.5 and 1.
  Endpoint values differ; unused whole endpoint arrays are NaN at 0 and 1.
  Covered source cells also retain NaNs that must not be fetched.
- A fraction mismatch on the last rank rejects globally, including when idle.
  A used root-ancestor endpoint NaN also rejects globally. Previously successful
  halo outputs remain exact in both failure cases.
- Existing source-aware MPI SDC hydro and split CTU checks remain in the same
  executable and passed.
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

Artifacts: `build/sdc-mpi-eb-source-time-checks/`.
The executed fixture uses outflow faces. No separate physical-boundary mode,
selected build, full matrix, installed-binary audit, level-selective coupled
SDC, recursive MPI EB application or physical validation was run.
