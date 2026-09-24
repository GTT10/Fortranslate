# MPI Mach inlet gas and mixture stepping — 2026-09-25

Implementation `3a15248`, based on `59a6f8f`.
`advance_mpi_spray_hydro`, `advance_mpi_spray_transport`,
`advance_mpi_spray_les` and `advance_mpi_mixture_reaction` append optional
`inflow_mach` and `inflow_pressure`. All physical ghost calls forward these
controls, as do both Dynamic/static LES and molecular half steps in the coupled
mixture routine. The common private validation helper collectively rejects
invalid or rank-disagreeing controls before periodic/no-op dispatch. Active
Mach requires an inlet face; inlet pressure requires a Mach array, and positive
pressure requires an active face. The mixture entry additionally requires hydro
when Mach is active. Existing candidate transactions, species limiting, flux
budgets and diagnostic contracts are preserved.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6 via cached `/usr/bin/mpiexec`, Ninja Debug.
Fixed test/runtime and selected H2O2 spray runtime built; selected execution was
not run. Final `spray_mpi_migration_2`: **1/1 passed**, 9.43 s, total 9.47 s,
no failures/skips. The earlier operator-only extension also passed (9.36 s,
total 9.40 s); the final run adds coupled Mach coverage.

The existing boundary comparison now adds an extrapolated x-low Mach 0.2 case
with moving/isothermal walls. Hydro, all molecular transport terms and static
SGS states, temperatures and global ledgers match serial operators within their
existing bounds (state 2e-12 scaled, temperature 1e-5 K hydro/1e-8 K transport
and SGS, surface exchange 1e-10 scaled). The new Mach case uses 1e-7 s; existing
cases retain 1e-5 s. Integrated changes match boundary fluxes. The coupled
mixture test now uses Mach 0.2 and 800 K at x-low: an escaping parcel is removed
and accounted, while gas evolution and the complete system budget pass.
Invalid wall temperature on one rank rolls back gas, temperature, parcel,
component masses, ID watermark and ledger after the attempted coupled work.
Dynamic Mach operator parity has separate preceding evidence; this coupled
execution uses the default static Smagorinsky model.

The first fixed test compilation failed on a newly added overlong Fortran call
line. Wrapping that line fixed compilation; no numerical change was needed.
Mechanism provenance, Python syntax of the contract checker, project contract
and whitespace passed. No existing tolerances changed. This is two-rank
operator/coupled API evidence, not MPI application input/restart, selected
execution, AMR Mach integration, nozzle parity or physical validation.

## Commands and artifacts

Logs: `build/mach-inlet-mpi-step-checks/`, including the fixed build failure and
successful rebuild, selected build, operator/final CTest logs and final
`LastTest.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_migration_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

MPI-related environment overrides are host-local workarounds, not defaults.
