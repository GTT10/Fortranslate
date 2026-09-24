# Regular MPI Dynamic LES physical fluxes and stepping — 2026-09-25

Implementation `e74daca`, based on `94d5ffd`.
`compute_mpi_dynamic_les_boundary_rhs_3d` connects distributed physical donors
to the ghosted Dynamic flux kernel. It returns local RHS, a global maximum
diffusivity and an optional global physical-boundary inventory rate. Only ranks
owning physical x faces contribute those faces; all ranks contribute their y/z
face areas. Fluxes at internal slab faces remain part of local divergence.
Global rate reduction may change floating-point summation order relative to serial.

The diffusion maximum includes each padded workspace and is conservative;
it need not equal the serial full-domain bound and may vary with partition.
Physical-face classification prevents wall heat flux from being applied at
internal slab edges. The all-periodic case retains the existing distributed
coefficient/filter operator. Model/filter, spacing, wall temperatures and faces
agree collectively; the donor API agrees inlet/wall-motion/pressure inputs.
RHS, maximum diffusivity and rate publish only on success.

`advance_mpi_dynamic_les_ssprk2_3d` appends optional boundary inputs and an
integrated boundary-change output. Both stages enforce the global diffusion
bound, recover temperatures collectively and publish only after complete
success. Regular MPI spray coupling/application dispatch remains to be connected.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug.
Fixed unit targets and selected H2O2 spray runtime built. Selected execution
was not run. Final CTest **3/3 passed**, total 11.04 seconds:

- Existing periodic `spray_les_mpi_dynamic_2`: 0.67 seconds.
- Extended `spray_les_mpi_boundary_2`: 5.25 seconds.
- Extended `spray_les_mpi_boundary_4`: 5.08 seconds.

The 5x3x2 test includes uneven/one-cell slabs and six periodic/physical boundary
combinations. All local RHS values equal serial exactly. Global boundary rates
match at 1e-12 relative maximum-rate scale. Reservoir inflow, moving isothermal
walls and pressure-outlet/moving-wall combinations also compare SSPRK2 states
and temperatures exactly, with integrated budgets at the same 1e-12 scale.
Excessive dt preserves accepted state, temperature and budget. Rank-disagreeing
wall-temperature controls preserve RHS, diffusion and rate sentinels. Existing
halo equality and width/wall/pressure disagreement checks remain.

An initial run also passed 3/3 (10.91 seconds). Inspection then added explicit
pressure-target validation before the periodic/physical dispatch. Both builds
and the focused checks were rerun for that final source. No failures, skips or
tolerance changes occurred. Mechanism provenance, Python syntax of the contract
checker, project contract and whitespace passed. This is distributed operator
evidence, not regular MPI spray application execution, physical validation,
external solver parity or a full configuration matrix.

## Commands and artifacts

Logs: `build/mpi-dynamic-physical-checks/` (initial and final build/test logs).
Launcher overrides below are host-local.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_boundary test_les_mpi_dynamic --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_mpi_(boundary_(2|4)|dynamic_2)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
