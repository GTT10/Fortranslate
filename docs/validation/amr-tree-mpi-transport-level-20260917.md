# MPI single-depth molecular-transport SSPRK2 — 2026-09-17

## Source and scope

Base: `c66df19050c4e7d57481af853b73fdb2e33a0c7b`, local `main`.
Adds an MPI transport level-step API, fixed/selected runtime wiring and focused
serial parity tests. Existing synchronous transport and serial numerical kernels
remain unchanged. Recursive level intervals/reflux remain the caller's work.

Verified SHA-256:

- `src/parallel/mpi_amr_patch_tree_transport_level_3d_mod.F90`:
  `676cd64458df93cf591e2c45d9c8da577c6f2468eadc1a149845454b5d0a9014`
- `tests/unit/test_mpi_amr_patch_tree_transport_level_3d.F90`:
  `db32a8dfb6a15752390cbf897cb91db6b03d166c046e55aae6c9dcc87b9368f2`

## Completed checks

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, existing clean system-MPI configuration
`build/selected-tree-mpi-debug`; compiler/backend options unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 4 --target test_mpi_amr_patch_tree_transport_level_3d pelef_selected_mpi_amr_tree_reactive_3d_runtime
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_mpi_amr_patch_tree_transport_level_3d_[23]$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Initial compilation rejected overlong API/test argument lines; line wrapping
corrected them. Final fixed/selected builds passed. CTest **2/2 passed**, 11.57 s
total: 2 ranks 5.76 s, 3 ranks 5.78 s; no failures/skips. Provenance, project
contract and whitespace checks passed. No Python sources changed. Existing
compiler warnings remain.

The three-level H2/O2 fixture has ratios 2 then 3. Depths 0/1/2 are tested with
viscosity, conduction, species diffusion and barodiffusion enabled, dt=1e-9 s,
PCM/MC donors, periodic boundaries, reservoir inflow and outflow with moving
isothermal walls. Independent stage fractions [0.1,0.25] and [0.6,0.75] sample
coarse endpoints whose conserved state changes by 2%; unused endpoint fields are
absent. Owner maps include empty ranks and all payload on one rank.

Every owned state and temperature field matches serial exactly, with inactive
fields unchanged. Active integrated x/y/z fluxes and the global minimum positivity
factor match serial exactly. Inactive flux slots are unallocated. Invalid diameter
on the active owner and valid but rank-dependent second-stage fractions trigger
collective rejection, preserving state, temperature, flux arrays and minimum
factor. No tolerance changes or golden hashes were used.

Logs: `build/amr-tree-mpi-transport-level-checks/{runtime-build,build,ctest}.log`.
The runtime-build log retains the initial line-length failure; build/ctest logs
contain the final successful incremental build/test.
Binary: `build/selected-tree-mpi-debug/tests/test_mpi_amr_patch_tree_transport_level_3d`.

## Limits

No averaging down, reflux, timestep enforcement or retry occurs in this API.
Species/transport coefficient identity remains a caller contract. MPI recursive
subcycling and coupled/application mode dispatch are not connected. No full
matrix, installed-prefix audit, scaling study, external solver parity or physical
validation was run.
