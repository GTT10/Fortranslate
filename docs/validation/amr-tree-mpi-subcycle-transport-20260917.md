# MPI recursive transport subcycling — 2026-09-17

## Source and scope

Base: `7a4353c6d608d21c63f82262f79c256dd9fa794b`, local `main`.
Adds owner-local recursive transport intervals, ratio-squared substeps,
coarse-time endpoints, child face accumulation, distributed interval reflux,
EOS recovery and scaled diffusion/isothermal-wall timestep bounds. The new API
is included in fixed and selected runtime libraries. Serial algorithms,
mechanism data, tolerances and golden hashes are unchanged.

Verified SHA-256:

- `src/parallel/mpi_amr_patch_tree_transport_subcycle_3d_mod.F90`:
  `0f69ccee27dfb5fc0b0d17e718a2492a71b2c7050f5eec2fc4bdcf9d90235c6d`
- `tests/unit/test_mpi_amr_patch_tree_transport_subcycle_3d.F90`:
  `bfec432bbe8b1c78c951945dda780b270f02e7c8cb9ed1ff0382bd4f948fc62c`

## Verification commands

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, existing system-MPI configuration
`build/selected-tree-mpi-debug`; compiler/backend settings unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_mpi_amr_patch_tree_transport_subcycle_3d pelef_selected_mpi_amr_tree_reactive_3d_runtime -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_mpi_amr_patch_tree_transport_subcycle_3d_[23]$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Final fixed/selected runtime builds passed. CTest **2/2 passed**, 101.20 s total
(2 ranks 50.30 s, 3 ranks 50.86 s), no failures or skips. Provenance, project
contract and whitespace checks passed. Existing compiler warnings remain.

The nonuniform three-level H2/O2 fixture has refinement ratios 2 then 3.
Four modes cover PCM/MC coarse donors, periodic boundaries, outflow,
reservoir inflow and tangentially moving isothermal walls. Ownership includes
split payloads and all payload on one rank with empty ranks. The interval is
1e-9 seconds, with viscosity, conduction, species diffusion and barodiffusion.
Assertions require exact serial equality of scaled dt, maximum diffusivity,
owned U/T, boundary budget, maximum reflux correction, minimum theta and
1/4/36 per-depth step counts. Invalid owner transport data and mismatched dt
must reject collectively and preserve state, temperatures, diagnostics and counts.
This test does not specifically inject failure during a later child interval.

Initial build failed on an overlong test line, which was wrapped. Initial tests
rejected the serial oracle because the fixture supplied nonzero wall velocity
on periodic faces. An intermediate positive-trace-species diagnostic did not
change that failure. The fixture now sets wall velocity only on physical walls
and retains its original zero-species composition; no solver guard was loosened.

Logs: `build/amr-tree-mpi-transport-subcycle-checks/` (`build.log`,
`test-build.log`, `initial-ctest.log`, `ctest.log`).
Binary: `build/selected-tree-mpi-debug/tests/test_mpi_amr_patch_tree_transport_subcycle_3d`.

## Remaining scope

Recursive reaction/transport/hydro composition and application/checkpoint mode
selection remain open. The advance API does not select/enforce dt or retry.
Coefficient identity remains the caller's contract. No full-field gather is
introduced; no production scaling study, full matrix, installed audit,
external-solver parity or physical validation was performed. No Python sources changed.
