# Distributed Mach inlet donors and MPI Dynamic LES — 2026-09-25

Implementation `95dedf2`, based on `9241732`.
`fill_mpi_les_boundary_ghosts_3d` appends optional `inflow_mach` and
`inflow_pressure`. It agrees Mach/pressure controls collectively, requires a
reservoir table only for ordinary inlet faces, and constructs Mach states from
owner-routed raw donor planes. The shared face mapper retains corner precedence
and moving-wall transforms. No global field gather is introduced. Candidate
ghosts and temperatures publish only after collective success.

`compute_mpi_dynamic_les_boundary_rhs_3d` and
`advance_mpi_dynamic_les_ssprk2_3d` forward the controls through both stages.
Active Mach on non-inflow faces and malformed pressure controls are rejected
before periodic dispatch. Physical controls, reservoir data and temperatures
are collectively agreed by the halo constructor. Existing global surface budgets,
atomic outputs and diffusion gates remain unchanged. The padded-workspace
maximum diffusivity remains a conservative bound and can depend on partition;
this check does not assert serial equality of that bound.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6 via cached `/usr/bin/mpiexec`, Ninja Debug.
Fixed MPI test/runtime and selected H2O2 spray runtime built; selected execution
was not run. Final CTests: **2/2 passed**, no failures/skips:
`spray_les_mpi_boundary_2` 9.08 s and `_4` 8.54 s, total 17.66 s.

The global 5x3x2 grid uses seven donor layers, including local one-cell x slabs.
Added cases cover an extrapolated x-low Mach inlet crossing moving y walls and
a prescribed z-high Mach inlet coexisting with an ordinary x-low reservoir.
All ghost states/temperatures, active RHS values and SSPRK2 states/temperatures
match full-domain serial results exactly. Global boundary rates/integrated
budgets meet the existing 1e-12 relative check. Excessive timesteps preserve
state and diagnostics. Mach-only input works without a reservoir table;
rank-disagreeing Mach, inlet pressure or temperature rejects collectively with
outputs unchanged. Existing periodic, reservoir, wall and pressure-outlet cases
also pass.

The initial run failed both rank-count cases at the serial pressure-donor check:
the new mixed-inlet fixture also imposed a 101000 Pa x-high outlet next to its
120000 Pa Mach inlet. Across seven layers the pressure correction can produce
inadmissible corner states (the kernel uses p + 2*layer*(target-p)). The new
Mach fixture now uses ordinary outflow; the pre-existing pressure-outlet case
remains separately exercised. No production acceptance checks or tolerances
were changed to accommodate the fixture.

Mechanism provenance, Python syntax of the contract checker, project contract
and whitespace passed. This is distributed donor/LES operator evidence, not MPI
hydro/transport/spray application or AMR Mach integration, nozzle parity or
physical validation.

## Commands and artifacts

Logs: `build/mach-inlet-mpi-les-checks/`, including initial/final CTest logs and
copied final `LastTest.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_boundary --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_mpi_boundary_(2|4)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

MPI-related environment overrides are host-local workarounds, not defaults.
