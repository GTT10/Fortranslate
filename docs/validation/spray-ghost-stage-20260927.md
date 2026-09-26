# Integrated ghost halo stage — 2026-09-27

Implementation `765bd9b4`, based on `83b1bfd9`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`advance_mpi_amr_ghost_halo` composes the existing halo kick/drift/reaction
transfer with redistribution of the evolved ghost population. It uses private
particle, component, node, source and ledger copies until final redistribution
succeeds. Coarse reaction routing is enabled. The caller supplies independent
kick/drift durations, stage-time gas halos, a target level, grow width and ghost
capacity. Optional EB metrics and exchange capacity reach the underlying paths.

The result is one stage for an existing ghost population. It does not establish
recursive timing, creation/removal lifecycle, breakup/children or the accounting
needed to combine physical, virtual and ghost feedback without duplication.
The recursive spray application remains unconnected.

## Completed verification

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. `spray_mpi_amr_step_2` and `_4` passed
(2/2, 0.87 s and 0.77 s; total 1.69 s).

The added fixture starts a ghost just outside a fine patch, moving inward.
The stage uses a 1e-8 s kick followed by a 1e-5 s drift. With final ghost
capacity zero, redistribution rejects after the preceding private operations;
checked positions/velocities, components/count, node indices, source/remainder
fields and sentinel boundary ledgers remain unchanged. With sufficient capacity,
the ghost survives inside the fine valid region, retains heating/evaporation,
and moves by the kicked velocity times the drift interval within 1e-15.
The interior stage has zero boundary exchange. Existing halo source conservation,
drift, outflow, capped and composite coupling assertions also execute.

Provenance, target build, Python syntax, project contract and whitespace checks
passed. No build/test failures or skips occurred. This new wrapper's EB path was
compiled but not independently exercised. No full matrix, installed audit,
recursive application, external solver parity or physical validation was run.

## Reproduction

Artifacts: `build/spray-ghost-stage-checks/`, including copied CTest details.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_step --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_step_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
