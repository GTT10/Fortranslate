# Regular MPI spray Dynamic physical boundaries — 2026-09-25

Implementation `7425c28`, based on `9599f3c`.
The regular MPI spray coupled step forwards physical faces, reservoir inputs,
wall temperature/motion and pressure targets into both Dynamic LES half steps.
Their integrated exchanges enter the existing gas/particle boundary ledger.
The coupled and application periodic-only guards are removed; the underlying
physical halo, flux and SSPRK2 APIs now provide those operations. Existing
candidate publication and whole-step timestep retries are retained. Positive
les_cs enables Dynamic SGS; zero still disables it.

The application checker adds `--dynamic-walls` for a mixed adiabatic moving wall
and isothermal wall fixture. Periodic axes remain periodic. A changed boundary
on restart now tests incompatible input context instead of an unsupported
Dynamic configuration. Existing physical inputs already belong to checkpoint
identity; no checkpoint format change is introduced.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug.
Fixed tests-enabled and selected H2O2 tests-disabled MPI spray applications
built. Selected execution was not run. CTest `spray_mpi_dynamic_wall_restart`:
**1/1 passed**, 20.01 seconds, total 20.05. No failures, skips or tolerance changes.

An 8x4x4 case uses Gaussian Dynamic filtering (test ratio 4, grid ratio 2),
four mixed-liquid parcels, an adiabatic no-slip wall moving tangentially at 2 m/s
and a 600 K isothermal wall. Net boundary energy is nonzero. Existing relative
mass/energy residual bounds of 2e-11 pass. Same-rank restart retains exact gas
and particle CSV values; 2-to-4-rank restart agrees at rel_tol=2e-11 and
abs_tol=1e-25. Model/filter changes and changed physical boundaries reject
checkpoint continuation, while repeated invocation preserves the saved history.
Other physical application combinations and selected execution were not rerun;
the previous distributed operator tests cover the narrower boundary operators.

Mechanism provenance, Python syntax for the application/contract checkers,
project contract and whitespace passed. This is focused coupled application and
rank-neutral restart evidence, not physical validation, external solver parity
or a full configuration matrix.

## Commands and artifacts

Logs: `build/mpi-dynamic-spray-wall-checks/`.
Artifacts:
`build/spray-amr-reactive-debug/tests/spray_mpi_dynamic_wall_restart/run-vhb90xtk/`.
Launcher environment overrides below are host-local.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_dynamic_wall_restart$' --output-on-failure
python3 -m py_compile tests/spray/test_mpi_spray_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
