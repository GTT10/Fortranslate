# MPI AMR LES Mach inlet integration — 2026-09-25

Implementation `1f3f618`, based on `2dbe5f2`.

## Change

The distributed LES timestep selector, common-step and per-level steppers now
accept Mach and inlet-pressure controls, forwarding them and inlet temperature
to each viscous ghost fill, including Dynamic's wide filter halo. Collective
validation and value agreement precede reservoir handling; only ordinary inlet
faces require reservoir tables. Recursive stepping forwards the same controls
through all intervals. Coupled application and input integration remain pending.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug. Fixed MPI
test targets and selected H2O2 spray-tree runtime built after the line-wrap
correction below. Selected execution was not run. **2/2 tests passed**, no
failures or skips:

- `spray_les_mpi_amr_subcycle_2`: 158.53 s (CTest total 158.57 s). Added Dynamic
  Mach 0.2 with extrapolated pressure and moving isothermal walls. State,
  temperature, boundary/reflux diagnostics, maximum diffusivity and 1/4/36 step
  counts match serial exactly. Existing four-model cases and rollback remain.
- `spray_les_mpi_amr_tree_4`: 148.13 s (CTest total 148.17 s). Added fixed-pressure
  Mach 0.2/120000 Pa inlet for all four models. Common-step and timestep parity
  pass across distributed and empty-rank ownership arrangements at existing
  tolerances (state 2e-12 scaled, temperature 1e-8 K, boundary budget 1e-14).
  Existing mismatch and rollback checks remain included.

Provenance, Python syntax, project contract and whitespace passed. No tolerance
or timeout changed. This is focused numerical parity, not physical validation.
The initial test regex included nonexistent `tree_3`; it ran only the recursive
test. The separate `tree_4` execution above supplies the common-step evidence.

## Commands and artifacts

Logs: `build/mach-inlet-mpi-amr-les-checks/`. The initial `build-fixed.log`
records a line-length compile failure, corrected by wrapping continued arguments.
Final builds are in `build-fixed-final.log` and `build-selected.log`.
Copied test logs are `ctest-subcycle.log` and `ctest-tree.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_amr_tree test_les_mpi_amr_subcycle --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_mpi_amr_(tree_3|subcycle_2)$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_mpi_amr_tree_4$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

MPI test environment overrides are host-local workarounds, not defaults.
