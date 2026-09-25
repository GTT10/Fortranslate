# MPI recursive/coupled EB viscosity — 2026-09-25

## Source and scope

Code `6f8f8ea`, parent `f510a7b`. MPI recursive EB hydro accepts finite
nonnegative explicit coefficients and agrees on them collectively before
stepping. Each recursive level requests complete state/volume-fraction halos
with the existing ancestor time interpolation, then passes the coefficient
into Euler flux correction before accumulation and reflux.

The full reaction/transport wrapper validates and collectively agrees on the
coefficient before private chemistry work, then forwards it to common or
recursive hydro. Molecular transport controls are unchanged. Omitted/zero
retains previous policy. EB application/default/restart policy and reference
geometry fidelity remain outstanding. The PCM ancestor geometry and inlet
corner policies retain the limitations in earlier EB halo evidence.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Ninja Debug, Python 3.12.3; reused
`build/spray-amr-reactive-debug`. Focused three-rank CTest **1/1 passed**,
43.16 s (43.21 s elapsed). No failures or skips. Compiler warnings remain.

The five-node z-plane fixture uses coefficient 0.037, Mach inflow, chemistry
and molecular transport together, recursive stepping, and explicit/implicit
reaction integrators. It covers distributed owners and all nodes on the last
rank. States, temperatures and boundary ledger match serial exactly; existing
wall/chemical/inventory comparisons pass without changed tolerances. A
rank-zero coefficient of 0.05 rejects with state/temperature/diagnostic rollback.
Existing invalid solver, inlet and mode-agreement failure checks also run.

Common-step full composition was compiled but not independently run here.
This is bounded coupled MPI/serial evidence, not a full configuration matrix,
external PeleC geometry parity or physical validation. Mechanism provenance,
Python syntax, project contract and whitespace checks passed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_subcycle_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_viscosity_subcycle_full_3d_3$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/mpi-eb-viscosity-split-checks/` contains `build.log`,
`rebuild.log` and copied `LastTest.log`.
