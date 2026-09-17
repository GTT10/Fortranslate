# Synchronized two-level EB molecular transport — 2026-09-17

Implementation `61a6732b3a75f647662f4d49114017c553a18f62`, based on
`734d832` on local `main`.

## Implementation

`advance_reactive_eb_transport_patch_euler_3d` builds limited transport fluxes,
applies their conservative divergence to active cells and recovers temperature
with the existing EOS. It returns dt times the actual limited flux per open
face area. State, temperature, flux object and minimum species-limiter factor
publish together only after success. Covered storage is retained. This plain
Euler patch API does not apply state redistribution.

`advance_reactive_eb_two_level_transport_3d` connects initial averaging, PCM donor
halos, both patch updates, area-weighted flux restriction, reflux and final
averaging through the same private synchronization sequence used by two-level
hydro. The existing public hydro signature remains a forwarding wrapper.
Only internal fine-patch faces enable donor transport. Physical domain faces
retain zero molecular flux; EB wall transport remains adiabatic/nonpenetrating
under the established kernel. The accepted composite physical-boundary inventory
and minimum patch species-limiter factor are returned with the state transaction.

The minimum factor describes local patch species limiting, not a stability proof
for reflux. Nonphysical synchronized states reject without changing either grid,
the boundary inventory or the limiter output. Callers must provide a suitable
small-cell diffusion timestep. The path remains serial, one static patch, common
Euler interval and PCM donors. Combined reaction/transport/hydro composition,
subcycling, multi-patch, hierarchy application/checkpoint and MPI remain open.
The donor-gradient limitations in the transport-halo record still apply.

## Completed focused verification

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6, Ninja; existing matching
`build/selected-tree-mpi-debug`, MPI on, spray off. Fixed core and selected EB
runtime compiled. No MPI or application matrix was run. No Python source changed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 2 --target test_eb_two_level_3d test_eb_two_level_chemistry_3d pelef_selected_reactive_eb_3d_runtime
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 2 --target test_eb_two_level_transport_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(unit_eb_two_level_3d|unit_eb_two_level_chemistry_3d|unit_eb_two_level_transport_3d)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

**3/3 PASS, 0 failures, 0 skips, total 3.74 s**:

- Existing two-level hydro: 0.30 s.
- Existing reaction/hydro composition: 2.90 s, including disabled-chemistry parity.
- New transport: 0.50 s. A cut-plane 4-cubed coarse domain with internal or
  boundary-touching fine patches, ratios two and three, three 1e-8 s steps each.
  Density/temperature, velocity and composition variation drives nonzero
  transport with viscosity, conduction, species diffusion and barodiffusion
  enabled. Composite inventory balances hold for every conserved component;
  physical-boundary transport is exactly zero. Refined parents match fine-cell
  inventories, temperatures stay positive and covered coarse storage is unchanged.
  Excessive timestep, geometry mismatch and invalid ratio preserve both grid
  states, the boundary inventory and limiter factor.

Provenance and source contract passed. Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Logs: `build/eb-two-level-transport-checks/{build,test-build,ctest}.log`.
This is focused conservation/transaction evidence, not convergence, external
PeleC parity, calibrated physical validation or general EB-AMR completion.
