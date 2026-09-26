# Coupled recursive MPI EB spray core — 2026-09-27

Implementation `a73a56c6`, based on `e5eb6167`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`advance_mpi_eb_spray_sdc_subcycle` joins the recursive gas SDC driver to the
physical/virtual/ghost population coordinator and complete particle predictor
halo composition. It constructs three-cell particle volume/centroid halos from
extended geometry, keeps private gas/temperature/reaction fields and population
snapshots, and provides the stage-time source and predictor callbacks.

Old particle sources/remainders feed every hydro predictor. Intermediate SDC
iterations retain those sources; final iterations evaluate the new half kick.
Per-level particle endpoint fields remove parent particle interpolation before
the actual peer/deposited particle halo is inserted. Completed physical
wall/outflow ledgers are accumulated once per level interval. After all clocks
reach the outer-step end, only physical populations are collected, with current
geometric owners and globally checked IDs. Gas, particles, components, source
inventories and optional reaction/level-count outputs publish together.

`particle_steps` counts accepted half-kick stages across level intervals.
The ID watermark is validated and retained; this core does not create children.
It supports frozen-gas drag/heating/evaporation with nonperiodic outflow or
stationary slip/no-slip boundaries and optional gas acceleration, chemistry,
transport and gas redistribution arguments. Invalid gas shapes, bounds and
extended metrics reject before publication.

Breakup models (including computational splitting), films, particle acceleration,
inflow/isothermal/periodic boundaries and application dispatch are not connected
to this new entry. Unsupported breakup/boundary input is rejected. End-of-level
subtree migration/cleanup and virtual aggregation still need reference-complete
integration; current physical level populations are collected at outer-step
synchronization. These remain work, not exclusions from the broader goal.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. Final `spray_mpi_eb_sdc_recursive_2` passed:
1/1, 7.29 s; CTest total 7.34 s. No build failures or test skips occurred.
Provenance, Python syntax, project contract and whitespace checks passed.

The two-level distributed planar EB case advances an interior fine-level mixed
water/peroxide particle with chemistry and all molecular transport flags enabled,
two SDC iterations, one parent interval and two child intervals. It checks six
accepted half-kick stages, level counts [1,2], final physical ownership/ID,
nonzero gas feedback, heating and evaporation. Composite gas-plus-particle
inventories close after boundary, wall, chemical and external changes with a
maximum relative residual of 1.5189353956769058e-15 (existing 3e-12 tolerance).
A negative hydro viscosity coefficient rejects after the old particle stage;
original gas, temperature, particle components, ID and sentinel ledgers/counts
are retained.

The initial run and one diagnostic repeat failed because the fixture requested
a positive computational splitting threshold with breakup_model='none', which
existing options validation rejects. The fixture now uses the supported model
and a hydro rejection for rollback; the unreachable splitting branch was removed.
No validation contract or tolerance was loosened. The corrected focused run and
the final run after input-validation fixes passed. Failure/debug logs are retained;
temporary trace statements were removed.

This evidence covers one interior-particle, two-level fixture. It does not prove
interface-crossing conservation, deep hierarchy lifecycle parity, physical spray
validation or complete application support. No full matrix, installed audit,
external solver numerical parity, restart or regrid run was performed.

## Reproduction

Artifacts: `build/spray-recursive-core-checks/`, including build variants,
initial/debug/retry/final CTest logs, provenance, contract and final `LastTest.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_recursive_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
