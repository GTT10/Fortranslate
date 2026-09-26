# Recursive population stage coordination — 2026-09-27

Implementation `e3a934d4`, based on `0a019c0f`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`advance_mpi_eb_recursive_population_stage` joins logical-level copy refresh
with the existing final-only SDC particle stage controller. The caller owns
physical, virtual and incoming-ghost arrays plus one stage controller per level.
On stage zero, the selected level and descendants must have completed particle
clocks at the requested start. For a child, its immediate parent must have
completed its provisional interval and that interval must contain the child
interval. Fresh controllers describe time zero. Clock comparisons allow 32
machine epsilons scaled by interval/start magnitude for floating-point scheduling.
MPI agrees time, interval, stage, iteration count and level before branching.

Before motion, the routine refreshes this level's virtual particles and creates
its next-level ghost population. Incoming ghosts are retained between child
substeps. It then runs the old/final-new particle stage against caller-supplied
stage-time gas/ancestor snapshots. Publication of the selected populations,
stage controller, outgoing ghosts, source rates and ledgers follows successful
completion of the whole call. Callers must keep the whole hierarchy private
until the encompassing gas step succeeds.

This is a callback-level coordinator, not the complete recursive gas/spray
wrapper. The caller still schedules intervals, supplies matching temporal gas
data, handles end-of-level cleanup/synchronization, constructs complete predictor
halos and connects application dispatch, breakup, films and forcing.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. `spray_mpi_eb_sample_2` and `_4` passed:
2/2, 0.74/0.74 s; CTest total 1.53 s. No test failures or skips occurred.
The initial build rejected an overlong new Fortran declaration line; splitting
that line fixed it. Both build logs are retained. Final build, provenance,
Python syntax, project contract and whitespace checks passed.

The nonperiodic two-level EB fixture initializes one physical particle, runs
one 2e-8 s parent interval and two 1e-8 s child intervals with two SDC iterations
each, and checks child clock advancement, retention of the incoming ghost
across the two child intervals and unchanged parent physical component inventory
during child updates. A child call before the parent finishes and a child
interval extending beyond the parent window both reject with sentinel ledgers
unchanged. Existing independent frozen-gas/temporal sampling, final-only source
scheduling and late population rollback checks still execute.

The schedule is invoked directly by the test with supplied gas snapshots;
actual gas SDC callback wiring, deeper hierarchies, full recursive gas/particle
conservation and end-of-level synchronization were not run here. No full matrix,
installed audit, external solver numerical parity or physical validation was
performed.

## Reproduction

Artifacts: `build/spray-recursive-population-checks/`, including initial/final
build, provenance, contract, CTest summary and copied detailed `LastTest.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_sample --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sample_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
