# MPI two-level PCM/PLM viscosity — 2026-09-25

## Source and scope

Code `b7fd2ba`, parent `eebc787`. MPI two-level AMR accepts finite nonnegative
PCM/PLM coefficients under its existing collective control agreement. Periodic
coarse hydro with a positive coefficient uses the full boundary halo path;
physical coarse fluxes receive the correction before their update. Fine face
fluxes receive it from the existing coarse-time-interpolated extended state at
each SSPRK2 stage, before flux averaging and reflux. No-slip masks include both
domain geometry and rank-local ownership. Empty fine owners retain the existing
zero-size path. Omitted/zero PCM/PLM and PPM policies remain unchanged.

Application input/default/restart policy and EB remain outstanding. This
increment exercises hydro; coupled chemistry/transport propagation is not
independently rerun. No full matrix or external PeleC/physical validation claim.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Ninja Debug, Python 3.12.3; reused
`build/spray-amr-reactive-debug`. Focused eight-rank CTest **1/1 passed**,
11.23 s (11.27 s elapsed). No failures or skips. Compiler warnings remain.

The fixture uses coefficient 0.037 with PCM and PLM. Periodic, interior physical
patch and wall-touching physical patch cases compare coarse/fine states,
temperatures and maximum reflux exactly against serial. Physical cases combine
slip/no-slip walls, inflow and pressure outlets; boundary-ledger tolerance
remains relative 2e-11 with the existing 1e-20 floor. Periodic large-step and
rank coefficient mismatch checks retain the hierarchy. Existing physical
rollback checks run too. No tolerances were changed.

After the passing run, the test selector was initialized to empty for calls
without a command argument and rebuilt successfully; the tested command always
overwrites it via `get_command_argument`, so that initialization does not alter
the tested path. No second runtime was needed for this test-only initialization.

Mechanism provenance, Python syntax, project contract and whitespace checks
passed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_reactive_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^mpi_amr_viscosity_hierarchy_3d_np8$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/mpi-two-level-viscosity-checks/` contains `build.log`,
`rebuild.log` and copied `LastTest.log`.
