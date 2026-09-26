# Level population SDC stages — 2026-09-27

Implementation `178018ce`, based on `33bb64ff`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`advance_mpi_eb_population_sdc_stage` and `spray_sdc_population_stage` provide
one level's old/new particle staging. Stage zero performs a half-interval kick
and full-interval drift on private physical/virtual/ghost copies, retaining the
midpoint populations. Stages 1 through the configured iteration count each
restart from that midpoint and perform one new half kick without drift. Only
the final iteration updates the caller's particle populations.

The controller validates stage sequence, iteration count, level, interval and
interval-start time. MPI agrees on numerical stage controls. Valid sources and
halo remainders are divided by the actual half-kick duration and returned as
rates, with overflow rejection. Stage-zero wall/outflow ledgers are stored;
later stages return old plus current-new physical boundary changes. These
ledgers are not sums over iterations. On success of the final iteration the
stage cursor returns to zero for reuse on another interval.

Failed population updates or rate conversion preserve the controller,
populations, rates and ledgers. The caller must keep this controller private
until its encompassing gas step succeeds. The caller also supplies correct
stage-time gas/ancestor fractions, refreshes virtual/ghost lifecycles and owns
breakup, films and particle acceleration. This does not yet connect the complete
recursive spray application or assemble its complete gas predictor halo.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. `spray_mpi_eb_sample_2` and `_4` passed:
2/2, 0.69/0.67 s; CTest total 1.41 s. No build/test failures or skips occurred.
Provenance, Python syntax, project contract and whitespace checks passed.

The nonperiodic planar EB fixture exercises all three populations. The old
stage uses ancestor fraction 0.25, a 1e-8 s half kick and a 2e-8 s drift directed
inward from a fine ghost-cell center. New stages use fraction 1 and no drift.
The midpoint position follows the kicked velocity; old and intermediate stages
retain the original published component inventories. Rates exactly match a
separate half-kick result divided by 1e-8 s. Two new iterations produce the same
final component masses, velocities and positions as one new half kick.

An invalid node in the cached third population forces failure during the final
iteration. The controller remains at that iteration and the prior rate fields
and sentinel boundary ledgers remain unchanged. Correcting the cached node and
retrying then succeeds. Earlier combined-stage rollback and independent
frozen-gas/temporal-halo checks still execute.

The check covers one level interval under two ownership sizes. Clock rejection,
multiple consecutive intervals, nonzero boundary-event ledger accumulation and
the full recursive gas/particle callback were not independently tested here.
No full matrix, installed audit, external solver numerical parity or physical
validation was performed.

## Reproduction

Artifacts: `build/spray-population-sdc-checks/`, including initial/final build,
provenance, contract, CTest summary and copied detailed `LastTest.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_sample --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sample_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
