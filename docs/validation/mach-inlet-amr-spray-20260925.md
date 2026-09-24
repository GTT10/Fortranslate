# MPI AMR reacting-spray Mach inlet forwarding — 2026-09-25

Implementation `ac60266`, based on `6f1a78a`.

## Change

The AMR spray composition accepts optional Mach/inlet-pressure controls and
forwards them to the gas reaction/transport/hydro composition and both LES half
steps, including recursive LES. Inlet temperature uses its existing argument.
Validation remains in the collective gas/LES operators, whose failure returns
before private gas/parcel candidates are published. Application input, timestep
selection and restart-context integration remain pending.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug. Fixed MPI
test and selected H2O2 spray-tree runtime built. Selected execution was not run.
`spray_mpi_amr_les_coupled_4` passed: **1/1 tests**, no failures or skips,
85.90 s (CTest total 85.94 s).

New cases use Mach 0.2/800 K, extrapolated pressure in common stepping and
120000 Pa in recursive stepping, with Vreman LES and molecular transport.
Gas-plus-parcel inventories close against the boundary ledger at the existing
3e-12 scaled bound. Invalid inlet temperature preserves gas, temperature,
parcel position/velocity/mass/temperature, component masses, IDs and the boundary
ledger. New Mach cases disable chemistry; the preceding reacting continuation
and existing periodic/physical-boundary cases remain included. Provenance,
Python syntax, project contract and whitespace passed. No tolerance or timeout
changed. This is focused coupled conservation evidence, not nozzle validation.

## Commands and artifacts

Logs: `build/mach-inlet-amr-spray-checks/`, containing `build-fixed.log`,
`build-selected.log` and copied `ctest.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_les_coupled_4$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Environment overrides are host-local workarounds, not application defaults.
