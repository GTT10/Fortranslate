# MPI spray Mach ghost entry — 2026-09-25

Implementation `4f290f0`, based on `02e5ea3`.
`fill_mpi_spray_ghosts` appends optional `inflow_mach` and `inflow_pressure`.
Active Mach faces require species thermodynamics. After collective face,
reflection-mode and Mach/pressure agreement, the entry routes active Mach
construction through the owner-routed donor-plane helper with width two.
Inactive Mach retains the existing nearest-neighbor ghost path. Candidate
outputs publish only on success.

`fill_mpi_les_boundary_ghosts_3d` now accepts optional `viscous_wall`, collectively
agreed and defaulting true for existing LES callers. The spray entry forwards
its existing default-false choice. Thus hydro retains normal-only wall
reflection; viscous boundaries retain full velocity reflection and moving-wall
energy correction. No whole-domain gather is added. The active Mach path can
fill one-cell x slabs, though the higher-level MPI steppers still retain their
existing extent requirements. This change only connects the ghost entry;
MPI hydro/transport/static LES and coupled application controls remain to be
forwarded.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6 via cached `/usr/bin/mpiexec`, Ninja Debug.
Fixed test/runtime and selected H2O2 spray runtime builds passed. Selected
execution was not run. CTests **2/2 passed**, no failures/skips:
`spray_les_mpi_boundary_2` 9.05 s, `_4` 8.65 s, total 17.74 s.

The added two-layer comparisons exercise both inviscid and viscous modes on a
5x3x2 global domain: an extrapolated x-low Mach inlet crossing moving y walls,
and a prescribed z-high Mach inlet coexisting with an ordinary x-low reservoir.
MPI conserved ghost states and temperatures match serial construction exactly,
including corners and local one-cell slabs. Existing seven-layer Dynamic RHS,
SSPRK2, boundary budget and disagreement/rollback checks also pass. No numerical
tolerances changed. Mechanism provenance, Python syntax of the contract checker,
project contract and whitespace passed.

This is boundary-construction evidence, not MPI coupled application execution,
AMR integration, external nozzle parity or physical validation.

## Commands and artifacts

Logs: `build/mach-inlet-mpi-spray-ghost-checks/`, including fixed/selected builds,
CTest and copied `LastTest.log`.

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
